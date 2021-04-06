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

"""Defines the Chemical Structure (CXS) class object.

The CXS class is the main container for loading, writing, and manipulating chemical structures;
The chemical structure objects might represent single molecules or collections of molecules.
"""
import sys
import os
import celeste.chemical_data as chemdata
import celeste.constants as const
import numpy as np
import copy
import scipy

class CXS:
    """Defines the chemical structure (CXS) object.

    Args:
        natoms (int) : Number of atoms
        positions (float, list) : Atomic positions, assumed to be in atomic units (``Bohr``) by default.
        atomlabels (str, list) : Atomic labels of all atoms.
        comment (str) : A text comment describing the structure.
        FileName (str, optional) : File from which structure may be read, if ``FileName is not None``.
        PositionUnits (str, optional) : Units for input atomic positions. Atomic units are used internally, so
            coordinates are converted after reading. Options for ``PositionUnits`` are ``Bohr`` or ``Angstrom``.
        Empty (bool): Flag indicating whether to create an empty ``CXS`` object.

    Attributes:
        natoms (int) : Number of atoms in structure.
        positions (float, list) : Atomic positions.
        comment (str) : Comment string describing the structure.
        mass (float) : Atomic mass list.
        nmol (int) : Number of molecules
        namol (int, list) : Number of atoms in each molecule. For example, ``namol[0]`` is the number of
                atoms in the first identified molecule.
        molid (int, list) : Atomic indices of atoms in each identified molecule. For example,
                ``molid[0][2]`` is the atomic index of the third atom in the first molecule.
        mol (Molecule) : A created list of ``Molecule`` objects, generated based on ``nmol``, ``namol``,
                and ``molid``.
        dist (numpy, float) : Interatomic distance matrix, where element (i,j) gives the interatomic distances
                between atoms i and j.
        graph (numpy, int) : Adjacency matrix for the structure, with entries 1 for bonded atoms and 0 for non-bonded.
                Note that the type of bonding does not matter here - we only report whether or not two atoms are bonded.
        TotalCharge (int) : Formal charge on structure...
        MolCharge (int, list) : Formal charge on each molecule in the structure.
        dvdr (float,list) : Derivatives of energy with respect to Cartesian coordinates.
        success (bool): Flag indicating whether last energy calculation performed for this structure
                was successful or not.

    Example:
        A CXS object can be created by either directly reading a structure from file:

    .. code-block:: python3

            a = CXS(FileName='test.xyz')

    or by defining the number of atoms ``natoms``, atomic ``positions`` and ``atomlabels`` parameters:

    .. code-block:: python3

        a = CXS(natoms = 2, positions =[(0.0,0.0,0.0), (1.5,0.0,0.0)], atomlabels =['C','O'],
                comment = "Carbon monoxide")

    The example above creates a carbon monoxide molecule, with the carbon atom at (0,0,0) and the
    oxygen atom at (1.5,0,0) (atomic units used by default).
    """

    def __init__(self, natoms=None, positions=None, atomlabels=None, comment=None, FileName=None,
                 PositionUnits='Bohr', Empty=False):
        """Initialization of a CXS object, either from file FileName, or from direct definition.
        """

        # If requested, create and empty CXS object....
        if Empty:
            self.natoms = 0
            self.positions = np.array([])
            self.dvdr = np.array([])
            self.velocities = np.array([])
            self.NEBForces = np.array([])
            self.PerpendicularForces = np.array([])
            self.ParallelForces = np.array([])
            self.atomlabels = []
            self.comment = "Empty CXS object"

        else:

            # If FileName is defined, read from file - otherwise, use passed values.
            if FileName is not None:
                # self.natoms,self.comment,self.atomlabels,self.positions = \
                self.ReadStructureFromFile(FileName)
            else:

                # Check number of atoms is at least sensible.
                if natoms is None:
                    print("\nError in cxs.__init__: Number of atoms must be defined.\n")

                try:
                    assert (natoms > 0)
                except AssertionError:
                    print("Error in cxs.__init__: Number of atoms must be greater than zero.")

                if atomlabels is None:
                    print("\nError in cxs.__init__: Atom labels must be defined.\n")

                self.natoms = natoms
                self.positions = np.array( positions )
                self.dvdr = np.zeros([natoms, 3])
                self.comment = comment
                self.atomlabels = atomlabels
                self.dvdr = np.zeros([self.natoms, 3])
                self.velocities = np.zeros([self.natoms, 3])
                self.PerpendicularForces = np.zeros([self.natoms, 3])
                self.ParallelForces = np.zeros([self.natoms, 3])
                self.NEBForces = np.zeros([self.natoms, 3])

        # Set the atomic masses.
        self.mass = self.SetMasses()

        # Calculate the adjacency matrix.
        self.GetGraph()

        # Calculate the number and details of molecules - note that this results in a new set of molecule
        # objects being created.
        self.GetMols()

        # Zero-out the energy:
        self.energy = 0.0000


    def SetMasses(self):
        """Allocates the atomic masses of each atom.

        The masses (in atomic units, where mass of electron is 1) are stored in a dictionary in ``chemical_data.py``.

        Returns:
             mass (float, list) : A list of atomic masses.
        """
        mass = np.zeros(self.natoms)
        for i in range(self.natoms):
            try:
                m = chemdata.AtomicMass[self.atomlabels[i]]
            except:
                print("Error in cxs.SetMasses - element not in chemical_data.py AtomicMass dictionary: "+self.atomlabels[i])
            mass[i] = m
        return mass


    def GetGraph(self):
        """Calculates the adjacency matrix ``self.graph`` for the structure.

        Notes:
            The determination of whether atoms are `bonded' or not is based on the covalent radii in the ``chemical_data``
            module, as well as the ``BondingScaleFactor`` in ``constants.py``. The cut-off distance :math:`R_{cut}`
            for each pair of atoms (i,j) is:

            .. math::
                R_{cut} = \\alpha (R_{i} + R_{j})

            where :math:`\\alpha` is ``BondingScaleFactor`` and :math:`R_{i,j}` are covalent radii.


        Example:
            After creating a CXS object, the adjacency matrix ``graph`` and the distance matrix ``dist``
            can be calculated as follows:

            .. code-block:: python3

                a = CXS(FileName='test.xyz')
                a.GetGraph()
        """
        n = self.natoms
        a = np.zeros(3)
        b = np.zeros(3)
        self.dist = np.zeros([n,n],dtype='float')
        self.graph = np.zeros([n,n],dtype='int')
        for i in range(n):
            self.dist[i,i] = 0.0
            radius_i = chemdata.CovalentRadii[self.atomlabels[i]]
            for j in range(i+1,n):
                a[:]= self.positions[:][i]
                b[:]= self.positions[:][j]
                self.dist[i,j] = np.linalg.norm(a - b)
                self.dist[j,i] = self.dist[i,j]
                radius_j = chemdata.CovalentRadii[self.atomlabels[j]]

                # Calculate the cutoff distance for this atom pair based on covalent radii
                # and the BondingScaleFactor (in constants.py).
                rcut = (radius_i + radius_j) * const.BondingScaleFactor

                if self.dist[i,j] <= rcut:
                    self.graph[i,j] = 1
                    self.graph[j,i] = 1


    def PrintEnergy(self,file=None):
        """Simple routine to print energy in pretty format.

        Args:
            file (str): File to print energy. If ``File is None``, print to stdout.

        """
        if file is None:
            out = sys.stdout
        else:
            out = file
        out.write("\n* Calculated energy = "+self.energy+" Eh\n")


    def Printdvdr(self,file=None):
        """Simple routine to print derivatives in pretty format.

        Args:
            file (str): File to print energy. If ``File is None``, print to stdout.

        """
        if file is None:
            out = sys.stdout
        else:
            out = file
        out.write("\n* Calculated derivatives:\n")
        for j in range(self.natoms):
            out.write('{} {} {:f} {:f} {:f} \n'.format(j,self.atomlabels[j],self.dvdr[j][0],self.dvdr[j][1],self.dvdr[j][2]))


    def PrintGraph(self,file=None):
        """Prints the adjacency matrix to either ``file`` or standard output.

        Args:
            file (str) : File to print adjacency matrix stored in ``self.graph``.

        Notes:
            If ``file is None``, then prints to standard output.
        """

        if file is None:
            out = sys.stdout
        else:
            out = file
        out.write("\n* NON-ZERO adjacency matrix elements for structure: "+self.comment+"\n")
        for i in range(self.natoms-1):
            for j in range(i+1,self.natoms):
                if self.graph[i,j] == 1:
                    out.write("Atom "+self.atomlabels[i]+" "+str(i)+" - Atom "+self.atomlabels[j] + " " +str(j)+"\n")
        out.write("\n")


    def GetMols(self):
        """Uses the adjacency matrix in ``self.graph`` to determine the molecules defined by the chemical
        structure object.

        This is achieved using the Floyd-Warshall shortest-path algorithm to split the adjacency matrix into
        disconnected sub-graphs.

        """

        # Get shortest-paths between atoms based on connectivity matrix.
        dsp = self.GetShortestPaths()

        # Get the molecules based on the shortest-path lengths.
        self.nmol = 0
        self.namol = [0 for i in range(const.NMOLMAX)]
        self.molid = [ [0 for i in range(const.NMOLMAX)] for j in range(const.NAMOLMAX)]
        self.MolCharge = [0 for i in range(const.NMOLMAX)]
        na = self.natoms

        ifound = [False for i in range(na)]

        while not all(ifound):

            for i in range(na):
                if ifound[i] is False:
                    # Started a molecule....
                    ifound[i] = True
                    j = self.nmol
                    self.namol[j] += 1
                    self.molid[j][0] = i
                    self.nmol += 1
                    break


            for j in range(na):
                if ifound[j] is False:
                    jfound = False
                    for k in range(self.nmol):

                        for l in range(self.namol[k]):
                            kl = self.molid[k][l]
                            if dsp[j,kl] < const.BIG_BONDING - 1:
                                self.molid[k][self.namol[k]] = j
                                self.namol[k] += 1
                                ifound[j] = True
                                jfound = True
                                break
                        if jfound is True:
                            break

        # Finally, assign molecule objects.
        nmol = self.nmol
        self.mol = []
        for i in range(nmol):
            na = self.namol[i]
            pos = []
            lab = []
            for j in range(na):
                ii = self.molid[i][j]
                x,y,z = self.positions[:][ii]
                pos.append([x,y,z])
                lab.append(self.atomlabels[ii])
            self.mol.append( Molecule(natoms=na,positions=pos,atomlabels=lab) )
            self.mol[i].formula = self.mol[i].GetFormula()


    def GetShortestPaths(self):
        """Determine the shortest-path between atoms using the Floyd-Warshall algorithm.

        Notes:
            If atoms i and j are not bonded, the distance between them is assigned some very
            large value ``BIG_BONDING`` which is stored in ``constants.py``.
        """
        n = self.natoms
        dist = np.zeros([n, n], dtype='float')
        for i in range(0,n-1):
            for j in range(i+1,n):
                if self.graph[i,j] == 1:
                    dist[i,j] = 1.0
                else:
                    dist[i,j] = const.BIG_BONDING
                dist[j,i] = dist[i,j]

        # Floyd-Warshall algorithm - note that dist is destroyed.
        for k in range(n):
            for  i in range(n):
                for j in range(n):
                    dist[i,j] = min(dist[i,j], dist[i,k]+dist[j,k])
        return dist


    def PrintAllMols(self,file=None):
        """Prints information on all molecules in structure.

        Args:
           file (str) : File to which molecule information is printed. If ``file is None``, then we print to stdout.
        """

        if file is None:
            out = sys.stdout
        else:
            out = file
        out.write("\n* Molecule information for structure: "+self.comment)
        for i in range(self.nmol):
            out.write("\n"+"+ Molecule: "+str(i)+"\n")
            out.write("- Number of atoms = "+str(self.namol[i])+"\n")
            out.write("- Molecular formula = "+str(self.mol[i].formula)+"\n")


    def ReadStructureFromFile(self,FileName):
        """Reads a ``CXS`` chemical structure object from file. Currently, only .xyz file formats are read.

        Args:
            FileName (int) : Input filename. For now, this should be a standard .xyz file.

        Returns:
            (tuple) : tuple containing

                - natoms (int): Number of atoms.
                - comment (str): A text comment, read from second line of xyz file.
                - atomlabels (str): List of atom labels.
                - positions (float, numpy array): x,y,z positions of all atoms.
        """

        # First check that FileName exists.
        try:
            assert(os.path.isfile(FileName))
        except AssertionError:
            print("Error: In ReadStructureFromFile, file '"+FileName+"' does not exist")


        # Check that this is identified as an xyz file.
        try:
            assert(FileName.lower().endswith('xyz'))
        except AssertionError:
            print("Error: In ReadStructureFromFile, FileName must end with .xyz")

        # Read the file.
        self.atomlabels = []
        self.positions = []
        with open(FileName) as f:
            self.natoms = int(f.readline())
            self.comment = f.readline()
            for line in f:
                line = line.strip()
                if len(line) != 0:
                    try:
                        a,b,c,d = line.split()
                    except:
                        print("Error: In ReadStructureFromFile, problem reading line '"+line)


                    # Check that b,c,d can be converted to floats.
                    #
                    try:
                        float(b)
                        float(c)
                        float(d)
                    except:
                        print("Error: In ReadStructureFromFile, (x,y,z) coordinates not floats in line: "+line)

                    self.atomlabels.append(a)
                    self.positions.append([float(b)*const.ang_to_bohr,float(c)*const.ang_to_bohr,float(d)*const.ang_to_bohr])

        # Check that the number of coordinate lines matches the number of atoms.
        try:
            assert( self.natoms == len(self.atomlabels) )
        except AssertionError:
            print("Error: In ReadStructureFromFile, number of atoms and number of coordinate lines don't match.")

        # Return values.
        self.positions = np.array(self.positions)
        return \
            # natoms,comment,atomlabels,np.array(positions)


    def PrinttoFile(self, file, format='xyz'):
        """Prints structure to file.

        Args:
            file (str): Filename to write the structure to.
            format (str): Output format. Allowed values: *xyz* only....
        """
        if format != 'xyz':
            print("* Only xyz supported in cxs.PrinttoFile right now...")
            sys.exit(1)
        elif format == 'xyz':
            f = open(file,"w")
            f.write( str(self.natoms)+"\n" )
            for i in range(self.natoms):
                x = self.positions[i][0] * const.bohr_to_ang
                y = self.positions[i][1] * const.bohr_to_ang
                z = self.positions[i][2] * const.bohr_to_ang
                line = "{}  {:14f}{:14f}{:14f}\n".format(self.atomlabels[i],x,y,z)
                f.write(line)
            f.close()

    def GetSpinandCharge(self):
        """Calculates the total charge and the spin multiplicity on the structure.
        """

        # First count electrons assuming neutral....
        self.nel = 0
        self.TotalCharge = 0
        for j in range(self.natoms):
            a = self.atomlabels[j]
            self.nel += chemdata.ElementSymbol[a]
        for j in range(self.nmol):
            self.nel -= self.MolCharge[j]
            self.TotalCharge += self.MolCharge[j]

        # Spin multiplicity.....based on odd or even electrons.
        if self.nel % 2 == 0:
            self.multiplicity = 1
        else:
            self.multiplicity = 2


    def CopyCXS(self):
        """Makes a copy of the self.CXS object in a safe manner.

        Returns:
            cxnew (``CXS``): A new ``CXS`` object which is a copy of ``self``.
        """
        # cxnew = CXS(Empty=True)
        # pos = self.positions
        # cxnew.atomlabels = self.atomlabels
        # cxnew.positions = pos
        # cxnew.mass = self.mass
        # cxnew.comment = self.comment
        # cxnew.natoms = self.natoms
        # cxnew.MolCharge = self.MolCharge
        # cxnew.graph = self.graph
        # cxnew.mol = self.mol
        # cxnew.molid = self.molid
        # cxnew.namol = self.namol

        cxnew = copy.deepcopy(self)

        return cxnew


    def NavigationFunction(self, cxp, knav = 1.0):
        """Evaluates the navigation function for MEP calculation, and calculates the derivatives.

        In the current implementation, the navigation function is similar to that used in IDPP, namely:

        .. math::
            V_{nav}(\mathbf{r}) = \sum_{i,j} k_{nav} ( r_{ij} - r_{ij}^{P})^2,

        Here, :math:`k_{nav}` is a spring-constant for the navigation function, :math:`r_{ij}^{P}`
        is the interatomic distance between :math:`i` and :math:`j` in the products (``cxp``), and :math:`r_{ij}` is
        the corresponding distance in the current structure (``self``).

        Notes:
            The navigation function is returned by this routine, and the derivatives with respect to the
            Cartesian coordinates are placed into ``self.dvdr_nav``.

        Returns:
            Enav (float): Value of the navigation function.
            cxp (``CXS`` object): Chemical structure object containing the target product structure.
        """

        self.Enav = 0.00000
        d = np.zeros(3)
        dp = np.zeros(3)
        dxr = np.zeros(3)
        n = self.natoms
        self.dvdr_nav = np.zeros([n,3])
        for i in range(n):
            for j in range(i+1,n):
                d[:] = self.positions[i,:] - self.positions[j,:]
                dp[:] = cxp.positions[i,:] - cxp.positions[j,:]
                dr = np.linalg.norm(d)
                dr_P = np.linalg.norm(dp)
                onr = 1.000/dr
                dxr[:] = d[:] * onr
                t1 = dr-dr_P
                self.Enav += knav * t1 * t1
                dE = 2.0 * knav * t1
                self.dvdr_nav[i,:] += dE * dxr[:]
                self.dvdr_nav[j,:] -= dE * dxr[:]

        return

    def NavigationDot(self):
        """Calculates the dot product between the derivatives of the PES and those given by the
        navigation function.

        Returns:
            The calculated dot product is stored in ``self.NavigationDotForces``.
        """
        mag1 = np.linalg.norm( self.dvdr )
        mag2 = np.linalg.norm( self.dvdr_nav )
        sum = 0.0000
        for i in range(self.natoms):
            for j in range(3):
                sum += self.dvdr[i,j] * self.dvdr_nav[i,j]
        self.NavigationDotForces = sum / (mag1 * mag2)

    def GetAllPairDistances(self):
        """Calculates ALL interatomic distances in the current structure, and returns them in
        a numpy array ``distances``. Also stores the distances in ``self.distances``.

        Returns:
            distances (np, float): n-by-n matrix containing all atomic pair-distance (in atomic units)
        """
        n = self.natoms
        distances = np.zeros([n,n],dtype=float)
        self.distances = np.zeros([n,n],dtype=float)
        for i in range(n):
            for j in range(i+1,n):
                distances[i,j] = np.linalg.norm(self.positions[i,:] - self.positions[j,:])
                distances[j,i] = distances[i,j]
        self.distance = distances
        return distances


    def KabschAlign(self, cxt):
        """Uses the Kabsch algorithm to align structure to a target structure ``cxt``. First, the centroids of the
        two structures are aligned (so the centroid of ``self`` is shifted to match that of ``cxt``), then the
        orientation of ``self`` is changed to match that of ``cxt`` using the Kabsch algorithm. This results
        in ``self`` having the minimized RMS displacement from cxt.

        Args:
            cxt (a ``CXS`` structure): A target structure, against which ``self`` is optimally aligned.

        """

        # Check that the number of atoms in both structures is the same.
        if self.natoms != cxt.natoms:
            print("* Error in cxs.KabschAlign: Number of atoms in two structures is NOT the same!")
            sys.exit(1)

        # Calculate centroid.
        n = self.natoms
        x = np.zeros(3)
        for i in range(n):
            x[:] += cxt.positions[i,:]
        x = x / float(n)

        # Shift the current structure to the centroid.
        for i in range(n):
            self.positions[i,:] -= x[:]

        # Finally, rotate self.positions to match cxt.positions using the Kabsch algorithm, as implemented
        # in scipy:
        a = cxt.positions
        b = self.positions
        Rotation, rmsd = scipy.spatial.transform.Rotation.align_vectors(a,b)

        # Finally, apply the rotation:
        self.positions = Rotation.apply(self.positions)



    def IDPPFunction(self,TargetDistances):
        """Calculates the IDPP potential function and its derivatives for the current structure.

        The ``TargetDistances`` need to be passed into this function; in IDPP, these are linearly
        interpolated between the reactants and products according to:

        .. math::
            r_{ij}^{target} = r_{ij}^{R} + \frac{k}{n} (r_{ij}^{R} - r_{ij}^{P})

        where :math:``r_{ij}^{R}`` is the inter-atomic distance between two atoms in the reactant structure,
        :math"``r_{ij}^{P} is the corresponding distance in the product.

        Args:
            TargetDistances (float): An [n,n] array of target distances, where n is the number of atoms.

        Returns:
            Eidpp (float): Calculated IDPP energy value.
            derivs (float): Derivatives of IDPP cost function, size [N,3].
        """

        n = self.natoms
        self.energy = 0.000
        self.dvdr = np.zeros([n,3])
        d = np.zeros(3)
        dxr = np.zeros(3)

        # Calculate function and derivs.
        for i in range(n):
            for j in range(i+1,n):
                d[:] = self.positions[i,:] - self.positions[j,:]
                dr = np.linalg.norm( d )
                onr = 1.000 / dr
                onr4 = onr**4
                dxr[:] = d[:] * onr
                t1 = dr - TargetDistances[i,j]
                self.energy += onr4 * t1 * t1
                dEdr = 2.0 * t1 * onr4 - 4.0 * t1 * t1 * (onr4 * onr)
                self.dvdr[i,:] += dEdr * dxr[:]
                self.dvdr[j,:] -= dEdr * dxr[:]
        return


    def CompareAtoms(self, cx):
        """
        Compares the number of atoms and the elements in the atom lists in ``self`` and ``cx``.

        Args:
            cx: A ``CXS`` object to compare against.

        Returns:
            Flag: A logical flag which is ``True`` if the total chemical formulae (i.e. number of
            atoms and element-types) are identical.
        """

        # First check number of atoms - return false if they are different.
        Flag = (self.natoms == cx.natoms)
        if Flag == False:
            return Flag

        # Now compare chemical formulae:
        Flag = ( self.GetFormula() == cx.GetFormula() )

        return Flag


    def GetFormula(self):
        """Creates a string representation of the molecular chemical formula.

        Returns:
            formula (string) : A string containing the chemical formula of the molecule.
        """
        nelem = [0 for i in range(chemdata.nelements)]
        for j in range(self.natoms):
            for key in chemdata.ElementSymbol:
                if self.atomlabels[j] == key:
                    nelem[chemdata.ElementSymbol[key]] += 1

        formula = ""
        for key in chemdata.ElementSymbol:
            if nelem[chemdata.ElementSymbol[key]] > 0:
                formula += str(key)
                formula += str(nelem[chemdata.ElementSymbol[key]])
        return formula


########################

class Molecule(CXS):
    """Defines the "molecule" sub-class. This is an inherited class (base-class is ``CXS``) defining
    individual molecules in a ``CXS`` object.

    Attributes:
        formula (string) : A string containing the chemical formula of the molecule.

    """

    def GetMols(self):
        """When GetMols is called for a molecule object, it should not return anything new - this function ensures
        that happens.....
        """
        self.nmol = 1
        self.GetFormula()
        return


    # def GetFormula(self):
    #     """Creates a string representation of the molecular chemical formula.
    #
    #     Returns:
    #         formula (string) : A string containing the chemical formula of the molecule.
    #
    #     """
    #     nelem = [0 for i in range(chemdata.nelements)]
    #     for j in range(self.natoms):
    #         for key in chemdata.ElementSymbol:
    #             if self.atomlabels[j] == key:
    #                 nelem[chemdata.ElementSymbol[key]] += 1
    #
    #     formula = ""
    #     for key in chemdata.ElementSymbol:
    #         if nelem[chemdata.ElementSymbol[key]] > 0:
    #             formula += str(key)
    #             formula += str(nelem[chemdata.ElementSymbol[key]])
    #     return formula





