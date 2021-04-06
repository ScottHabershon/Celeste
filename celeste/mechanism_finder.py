#
############################################################################
#
# Celeste - chemical reaction discovery code.
# Copyright (C) 2021, Scott Habershon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
############################################################################


"""Defines the Mechanism object and related functions.

The Mechanism object is the main container for performing either double-ended mechanism-finding or
single-ended mechanism generation.
"""

import sys
import celeste.cxs as cx
import numpy as np
import constants as const
import json
import copy
import chemical_data as chemdata

class Mechanism:
    """Defines the ``Mechanism`` object as a container for a double-ended reaction mechanism.

    This module contains the ``Mechanism`` object, as well as related mechanism-finding routines. These routines
    implement the double-ended mechanism-finding approach, as well as a single-ended
    mechanism generation scheme, both as developed in the Habershon group.

    Args:
        reactantfile (str) : Structure containing the initial reactant file.
        productfile (str,optional) : Structure containing the target product file (used in the case of
            double-ended mechanisms search).
        maxreactions (int) : Maximum number of allowed reactions in the mechanism (including
            non-progressing null reactions).
        valencerange (list of characters and ints): Defines the allowed valence ranges of atoms in the system,
            used to judge whether graphs generated along a proposed reaction-mechanism are sensible or not.
            For example, ``valencerange = [['C',1,4],['O',1,2]`` defines the allowed carbon valence to be in the
            range 1-4 and the allowed oxygen valence to be in the range 1-2.
        reactiveatoms (list): A list of atomic labels indicating those atoms which are allowed to participate in
            reactions. If an atom label is not included in this list, it can't be used in reactions.
        errortype (string): Defines the graph-error function to be used, either ``FullStructure`` or ``SingleMolecule``.

    Attributes:
        maxreactions (int): Maximum number of reactions in reaction-mechanism. Note that this includes ``null``
           reactions whcih don't progress the mechanism.
        ReactiveAtoms (list): This is the local storage of the input ``reactiveatoms`` list.
        ns (int): This is the total number of structures along a given reaction mechanism,
        nrxn (int): Number of reactions in mechanism (same as maxreactions).
        cxr (CXS object): Chemical structure object defining the reactants.
        cxp (CXS object): Chemical structure object defining the target products (for double-ended search).
        ValenceRange (list): Local copy of ``valencerange`` input list.
        natoms (int): Number of atoms in the structures considered, as taken from the reactant structure.
        atomchange (logical list): List of True/False defining whether atoms are allowed to react. This list is
            ``natom`` atoms long for an ``natom``-size system. The ``atomchange`` list is based on the input
            ``reactiveatoms`` list.
        rxindex (int list): List of the indices of reactive atoms in the system, as defined by ``atomchange``.
        ngmoves (int): Number of defined moves.
        namove (list of ints): Number of modified atoms in each move-class.
        move (3-d list): For each move, ``move`` is the graph defining the change.
        moveid (list of strings): The ``id`` string for each move read from the ``MoveFile``.
        ErrorType (string): Defines the graph-error function, either ``FullStructure`` or ``SingleMolecule``.
        structures (CXS objects): The set of ``CXS`` objects created by application of the reaction-moves
            in self.movenum to the reactive atoms in self.moveatoms. Including the reactants, the total number of entries in structures[] is self.nrxn - we
            need to make sure this is accounted for in later indexing.

    """

    def __init__(self, reactantfile = None, productfile = None, maxreactions = 1,
                 valencerange=None, reactiveatoms = None, errortype = None):
        """Initializes a ``Mechanism`` object.
        """

        # Make local copy of errortype:
        self.ErrorType = errortype

        # Convert the string of reactive atoms into a list:
        if reactiveatoms is None:
            print("* ERROR in mechanism_finder.__init__: reactiveatoms is None.")
        else:
            self.ReactiveAtoms = json.loads(reactiveatoms)

        # Total maximum number of structures in the mechanism, including reactants plus all
        # intermediates generated by maxreactions reactions:
        self.ns = 1 + maxreactions
        self.maxreactions = maxreactions

        # Store maximum number of reaction steps:
        self.nrxn = maxreactions

        # First, create reactant structure.
        self.cxr = cx.CXS(FileName=reactantfile)
        if productfile:
            self.cxp = cx.CXS( FileName = productfile )

        # Check that number of atoms and elements in reactants and products match.
        if productfile:
            Flag = self.cxr.CompareAtoms( self.cxp )
            if Flag == False:
                print("* ERROR in mechanism_finder.py - reactant and product structures have different chemical formulae")
                sys.exit(1)

        # Set valence ranges and check all atoms have identified valence ranges:
        self.ValenceRange = valencerange
        for label in self.cxr.atomlabels:
            ok = False
            for j in range(len(self.ValenceRange)):
                if label in self.ValenceRange[j]:
                    ok = True
            if ok != True:
                print("* ERROR in mechanism_finder.py (__init__) - atomlabe does not have a defined ValenceRange")
                sys.exit(1)

        # Finally, set the atomchange list, which is a logical flag list of length self.cxr.natoms,
        # defining whether or not the atom cna react (at least as judged by its atomlabel and the
        # ReactiveAtoms input flag):
        self.atomchange = []
        self.natoms = self.cxr.natoms
        for i in range(self.natoms):
            if self.cxr.atomlabels[i] in self.ReactiveAtoms:
                self.atomchange.append(True)
            else:
                self.atomchange.append(False)

        # Set up a list of atoms which can react, to help with random selection later:
        self.rxindex = []
        for i in range(self.natoms):
            if self.atomchange[i] == True:
                self.rxindex.append(i)

        return

    def ReadMoveSet(self, movefile = None):
        """
        Reads the set of allowed reaction moves for a mechanism search simulation.

        Each move must have the following format:

        .. code-block::

            *move
            id: A string describing the reaction
            namove: Number of atoms in move
            graph: number of lines following must be same as na
            0 1
            1 0
            *end

        Each move has to begin with ``*move``, has an ``id``, and then has the move-graph as a list.

        Args:
            movefile: A file containing the moves.

        """

        # Check that the move-file exists.
        try:
            f = open(movefile)
        except IOError:
            print("Move-file named " + movefile + " does not exist...")
            sys.exit(1)


        # Read the move-file.
        data = f.readlines()
        block = []
        blocks = []
        found = False
        for line in data:

            # Are we reading a move-block?
            if found and not line.startswith('*end'):
                block.append(line.rstrip('\n'))

            # If line begins with #, then ignore it as a comment line.
            if line.startswith('#'):
                continue
            elif line.startswith('\n'):
                continue
            elif line.startswith('*move'):
                found = True
            elif line.startswith('*end'):
                found = False
                blocks.append(block)
                block = []


        # Setup data arrays containing moves.
        self.ngmoves = len(blocks)
        self.namove = np.zeros(self.ngmoves, dtype='int')
        self.move = np.zeros([self.ngmoves, const.NAMOVEMAX, const.NAMOVEMAX ])
        self.moveid = []
        ic = 0
        for block in blocks:

            # First line of block contains id:
            key, value = block[0].split(":",maxsplit=1)
            if key != "id":
                print("* ERROR: First line in move block should contain id.")
            else:
                self.moveid.append(value)

            # Second line is number of atoms:
            key, value = block[1].split(":", maxsplit=1)
            if key != "na":
                print("* ERROR: Second line in move block should contain na.")
            else:
                self.namove[ic] = int(value)

            # Second line is number of atoms:
            key, value = block[2].split(":", maxsplit=1)
            if key != "graph":
                print("* ERROR: Third line in move block should contain graph.")
            else:
                for k in range(self.namove[ic]):
                    a = block[3+k].split()
                    for m in range(self.namove[ic]):
                        self.move[ic,k,m] = a[m]
            ic += 1

        # Convert the move array into an numpy array:
        self.move = np.array(self.move)

        return

    def ReportMoves(self, File=None):
        """Creates a report of the reaction moves

        Args:
            File: An output file (open already).
        """

        if File is None:
            out = sys.stdout
        else:
            out = File

        out.write("\n*** Available graph-move report *** \n")
        n = self.ngmoves
        for i in range(n):
            out.write("\n===========\n")
            out.write("\n- Move number: "+str(i)+"\n")
            out.write("- ID:"+self.moveid[i]+"\n")
            out.write("- Number of atoms: "+str(self.namove[i])+"\n")
            out.write("- Move graph: \n")
            for j in range(self.namove[i]):
                out.write(" "+str(self.move[i,j,0:self.namove[i]]) + "\n")
            out.write("\n===========")

        return

    def RunMechanismFinder(self, iterations = 100, InitialTemperature = 100.0, progressfile = None,
                           File = None, OutputFrequency = 500):
        """Performs a double-ended reaction-mechanism search.

        Args:
            iterations (int): Maximum number of SA/MC iterations
            InitialTemperature (float): Initial temperature for SA simulation (in K).
            progressfile (str): A filename where reaction-mechanism search progress will be output.
            File (str): A filename where general output is placed.
            OutputFrequency (int): Frequency with which output is generated during simulated annealing run.

        Returns:
            converged (logical): This is true if a mechanism was correctly found within the
                                 required number of iterations.

            At the end of the simulation, if ``converged == True``, then:

                - self.movenum[:] contains the list of integer graph-moves at each reaction-step.
                - self.moveatoms contains the list of atoms that each graph-move applies to.

        """

        # Setup the output file-stream:
        if File is None:
            out = sys.stdout
        else:
            out = File

        # Write a header to the output file:
        out.write("\n\n* Double-ended reaction-mechanism search starting...\n")


        # Open the progress output file and write the header:
        fprogress = open( progressfile, "w")
        fprogress.write(
            "\n-----------------------------------------------------------------------------------\n")
        fprogress.write("# Iteration     |    Graph error function (au)    |   SA Temperature (K) \n")
        fprogress.write(
            "\n-----------------------------------------------------------------------------------\n")


        # Set up the arrays for the move numbers and reactive atoms.
        #
        nrxn = self.nrxn
        self.movenum = np.zeros(nrxn, dtype = 'int')
        self.moveatoms = np.zeros( [nrxn, const.NAMOVEMAX], dtype = 'int')


        # Now, the identity movenum = -1 corresponds to a "null" reaction whcih does
        # not change the input graph - set all initial reactions to the "null" reaction:
        #
        self.movenum[:] = -1

        # Evaluate the initial path error:
        Error = self.GetPathError()

        # Print information on initial path-error:
        fprogress.write("  0   " + str(Error))
        out.write("\n- Initial graph error = "+str(Error)+"\n")

        # Loop over SA iterations:
        converged = False
        for it in range(iterations):

            # Check convergence:
            if Error < const.GCONV:
                out.write("\n* Mechanism search converged...\n")
                break

            # Calculation the temperature for this step - assume linear scaling:
            Temperature = InitialTemperature * (1.0 - (float(it)/float(iterations)))
            beta = 1.0 / (const.kboltz * Temperature)

            # Store the old error before graph update:
            Error_old = Error

            # Store the mechanism:
            movenum_store = copy.deepcopy( self.movenum )
            moveatoms_store = copy.deepcopy( self.moveatoms )

            # Make an update to the mechanism.
            self.UpdateMechanism()

            # IF something has gone horribly wrong, reject this update:

            # Evaluate the new error:
            Error = self.GetPathError()

            # Decide whether to accept new mechanism or not - the switching here is to avoid
            # exponential overflow (np.exp can handle up to about exp(709)).
            if Error - Error_old > 700.0:
                P = 0.0
            elif Error - Error_old < 0.0:
                P = 1.0
            else:
                P = min(1.0, np.exp(-beta * (Error - Error_old)))
            if np.random.uniform(low=0.0, high=1.0) > P:
                self.movenum = copy.deepcopy( movenum_store )
                self.moveatoms = copy.deepcopy( moveatoms_store )

            # Print out some useful information:
            if it%OutputFrequency == 0:
                fprogress.write(str(it) + "     " + str(Error)+"\n")
                out.write(str(it) + "     " + str(Error)+"\n")

        # If we get here and the error is too large, alert the user.
        if Error > const.GCONV:
            converged = False
            out.write("\n* Mechanism search completed without locating mechanism...\n")
        else:
            converged = True

        # End of routine - return flag defining whether calculation converged or not.
        return converged


    def UpdateMechanism(self):
        """Performs a trial move by modifying a current mechanism string by either changing
         a move and its atoms, or just changing the atoms of an existing move.

        .. note::

            When updating the atoms, we need to be careful to select atoms in a bonding pattern
            which is consistent with the requirements of the bonding-graph to be performed for the
            relevant move.

            For example, if the chosen move requires a bond formation between atoms *i* and *j*, then
            we can't select atoms *i* and *j* which are **already** bonded!

        """

        # Set flag indicating success in searching for atoms.
        found = False

        # First select a reaction to change:
        irx = np.random.randint(low=0,high=self.nrxn)

        # Decide which move-type to perform. If the movenum is -1, it's a null reaction so we can
        # only perform a move which selects a new reaction-class and a new set of atoms:
        if self.movenum[irx] == -1:
            movetype = 2
        else:
            rx = np.random.uniform(low=0.0,high=1.0)
            if rx < 0.5:
                movetype = 1
            else:
                movetype = 2

        # Try to change the movenum and moveatoms for reaction irx.
        #
        # movetype = 1: Update the atoms of the move only.
        #
        if movetype == 1:

            imove = self.movenum[irx]
            nx = self.namove[imove]
            nrx = len(self.rxindex)

            # We'll try up to const.NTRYMAX times to sample new reactive atoms.
            found = False
            for itry in range(const.NTRYMAX):

                # Select random atoms - ix is a list of positions in the rxindex list of atoms:
                ix = []
                for j in range(nx):
                    ix.append(np.random.randint(low=0,high=nrx))

                # Check non of the selected atoms in ix are the same:
                duplicates = any(ix.count(element) > 1 for element in ix)

                # If there are no duplicates.....
                if duplicates == False:

                    # Check that the selected atoms are viable reactive atoms
                    # based on the reaction graph:
                    found = True
                    for i in range(nx):
                        ii = ix[i]
                        for j in range(nx):
                            jj = ix[j]
                            gx = self.move[imove, i, j]
                            if gx == 1 and self.structures[irx-1].graph[ii, jj] == 1:
                                found = False
                            if gx == -1 and self.structures[irx-1].graph[ii, jj] == 0:
                                found = False

                # If we get here, and found == True, we're fine!
                if found == True:
                    for k in range(nx):
                        self.moveatoms[irx, k] = int(ix[k])
                    break

        # If movetype == 2, then we change both the move number and the reactive atoms:
        elif movetype == 2:

            # Select a new move number for reaction irx:
            imove = self.movenum[irx]
            while imove == self.movenum[irx]:
                imove = np.random.randint(low=-1,high=self.ngmoves)

            # Select new reactive atoms:
            if imove == -1:
                self.moveatoms[irx, :] = 0
                self.movenum[irx] = imove
                found = True
            else:

                nx = self.namove[imove]
                nrx = len(self.rxindex)

                # We'll try up to const.NTRYMAX times to sample new reactive atoms.
                found = False
                for itry in range(const.NTRYMAX):

                    # Select random atoms - ix is a list of positions in the rxindex list of atoms:
                    ix = []
                    for j in range(nx):
                        ix.append(np.random.randint(low=0, high=nrx))

                    # Check non of the selected atoms in ix are the same:
                    duplicates = any(ix.count(element) > 1 for element in ix)

                    # If there are no duplicates.....
                    if duplicates == False:

                        # Check that the selected atoms are viable reactive atoms
                        # based on the reaction graph:
                        found = True
                        for i in range(nx):
                            ii = ix[i]
                            for j in range(nx):
                                jj = ix[j]
                                gx = self.move[imove, i, j]
                                if gx == 1 and self.structures[irx - 1].graph[ii, jj] == 1:
                                    found = False
                                if gx == -1 and self.structures[irx - 1].graph[ii, jj] == 0:
                                    found = False

                    # If we get here, and found == True, we're fine!
                    if found == True:
                        for k in range(nx):
                            self.moveatoms[irx,k] = int(ix[k])
                        self.movenum[irx] = imove
                        break

        return

    def GetPathError(self):
        """"Calculates the permutationally-invariant graph-error function for the current
        mechanism.

        The current mechanism is defined through:
            - movenum[:]: containing the integer list of reaction-moves at each step;
            - moveatoms[:,:]: containing the integer list of atom indices
        """

        # First propagate the graph to get the set of structures from this mechanism.
        errflag = self.PropagateGraph()

        # If the error flag errflag is True, it means that an error was detected with the
        # proposed reaction mechanism, so we assign an arbitrarily large error:
        if errflag == True:

            Error = const.LARGE_ERROR

        else:

            # This ErrorType compares the mass-weighted eigenvalues of the full structure:
            if self.ErrorType == 'FullStructure':

                # Create mass-weighted graph:
                n = self.natoms
                mgraph1 = np.zeros([n,n])
                mgraph2 = np.zeros([n,n])
                mgraph1 = self.structures[self.maxreactions].graph
                mgraph2 = self.cxp.graph
                for i in range(n):
                    # Get atomic number.
                    mm = chemdata.ElementSymbol[self.cxp.atomlabels[i]]
                    mgraph1[i,i] = mm
                    mgraph2[i,i] = mm

                # Diagonalize mass-weighted graphs:
                Eig1,C1 = np.linalg.eigh(mgraph1)
                Eig2,C2 = np.linalg.eigh(mgraph2)

                # Compare eigenvectors:
                Error = np.linalg.norm(Eig1-Eig2)
                #Error = 0.0
                #for i in range(n):
                #    Error += (Eig1[i]-Eig2[i])**2

            # This ErrorType
            elif self.ErrorType == 'SingleMolecule':
                sys.exit(1)


        return Error


    def PropagateGraph(self):
        """For the current mechanism (as encapsulated in self.moveatoms and self.movenum), propagate
        for self.nrxn reaction-steps and return the final graph.

        Returns:
            errflag: A logical flag which is True if there is an error detected.
        """



        # Set the current structure:
        current = copy.deepcopy( self.cxr )

        # Create an array containing the returned structures.
        self.structures = []
        self.structures.append(current)

        # Set the errorflag, indicating a forbidden move:
        errflag = False

        # Loop over reactions...
        for irxn in range(self.nrxn):

            # Copy the current structure.
            new = copy.deepcopy( current )

            # Perform the reaction - a movenumber of -1 is a "null" reaction which does nothing,
            # so we ignore those:
            ir = self.movenum[irxn]
            if ir != -1:

                # Create a list containing the atoms reacting in this move:
                ix = []
                nx = self.namove[ir]
                for i in range(nx):
                    ix.append(self.moveatoms[irxn,i])


                # Check allowed reactivity of atoms, as defined in ReactiveAtoms from input file:
                for i in range(nx):
                    ii = ix[i]
                    if self.atomchange[ii] == False:
                        errflag = True
                        break

                # TODO: Add bondchange check.

                # Check that the selected atoms in self.moveatoms for this reaction actually
                # match the intended graph.
                for i in range(nx):
                    ii = ix[i]
                    for j in range(nx):
                        jj = ix[j]
                        gx = self.move[ir,i,j]
                        if gx == 1 and new.graph[ii,jj] == 1:
                            errflag = True
                        if gx == -1 and new.graph[ii,jj] == 0:
                            errflag = True

                # If we get here and errflag == False, then we can go ahead and apply the reaction:
                if errflag == False:
                    for i in range(nx):
                        ii = ix[i]
                        for j in range(nx):
                            jj = ix[j]
                            new.graph[ii,jj] += self.move[ir,i,j]

                # Check the valencies of the atoms in the new structure:
                for i in range(self.natoms):
                    CurrentValence = sum(new.graph[i,:])
                    found = False
                    for j in range(len(self.ValenceRange)):
                        if new.atomlabels[i] == self.ValenceRange[j][0]:
                            found = True
                            if CurrentValence < self.ValenceRange[j][1] or CurrentValence > self.ValenceRange[j][2]:
                                errflag = True


            # Append the current updated structure:
            self.structures.append( new )

            # Finally, copy the new graph to be the start point for the next reaction:
            current = copy.deepcopy( new )

        return errflag

