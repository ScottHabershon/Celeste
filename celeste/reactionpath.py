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


"""Defines the Reaction-Path class.

The ``RP`` class is the main container for creating, reading, writing and
manipulating reaction paths, comprising collections of chemical structure (``CXS``)
objects.
"""
import sys
import celeste.cxs as cx
import numpy as np
import celeste.pes as pescalc
import celeste.constants as const
import scipy.optimize as opt
from celeste.pes import PES
import celeste.input_output as io


class RP:
    """
    Defines the reaction-path class object.

    Initializes a ``RP`` object, comprising a set of chemical structures along a reaction-path.

    If all of ``reactantfile``, ``productfile`` and ``pathfile`` are defined simulataneous,
    we're note sure which to use, so we throw an exception.

    Args:
        reactantfile (str): File containing reactant structure.
        productfile (str): File containing product structure.
        pathfile (str): File containing the full reaction-path.
        nimages (int): Number of images in the reaction path (including reactants and products).
        OptimizeEnds (bool) : Flag indicating whether or not the end-points of a reaction-path should be
                              geometry-optimized.
        EnergyCalc (PES object list): PES object list defining a set of energy calculation objects. These are used
                             if ``OptimizeEnds == True``.
        EmptyPath (bool): Flag indicating whether or not to create an empty path, containing just  the reactant
                        and product structures.
        AlignEndPoints (bool): Flag indicating whether or not to align the product structure to the reactant structure.
                         If ``True``, this alignment is achieved by translatin centroids, then minimizing the RMS structural
                         difference using the Kabsch algorithm. Note that ``AlignEndPoints`` only works if we are
                         initializing from ``reactantfile`` and ``productfile`` (not ``pathfile``).

    Examples:
        This class provides two methods of initializing ``RP`` objects. First, if ``pathfile`` is defined, we start from
        the structurs in the given ``pathfile``.

            .. code-block:: python3

                a = RP(pathfile='path.xyz')

        The file ``path.xyz`` (in this case) should contain a sequence of xyz-format frames, one for each
        image in the reaction-path. In this case, the structure and number of images ``nimages`` are read from file.

        Alternatively, if the reactant structure and product structure are given (in ``reactantfile`` and ``productfile``,
        respectively), an empty path can be generated (then used for later interpolation. Here, the TOTAL number
        of images in the path is then ``nimages``.

            .. code-block:: python3

                a = RP(reactantfile='reactants.xyz',productfile='products.xyz', nimages = 10)

      """

    def __init__(self,reactantfile=None,productfile=None, pathfile=None, nimages=0,OptimizeEnds=False,
                 EnergyCalc=None, EmptyPath = False, AlignEndPoints = False):
        """
        Initializes a reaction-path object ``RP``.
        """

        # If EmptyPath == True, we just create an empty path with only the reactant and produc structure.
        # This empty path is useful in some kinds of calculations, such as MEP approximation using
        # navigation functions or the growing-string method.

        if EmptyPath:

            self.nimages = 1
            self.start = 0
            self.end = 1
            try:
                assert(reactantfile is not None)
            except:
                print("\nError in reactionpath.__init__: pathfile is None, but no reactantfile is given\n")
            try:
                assert(productfile is not None)
            except:
                print("\nError in reactionpath.__init__: pathfile is None, but no productfile is given\n")

            # Read the reactant structure.
            self.rpcxs = [cx.CXS(FileName=reactantfile)]

            # Read the product structure.
            self.rpcxs.append(cx.CXS(FileName=productfile))

        else:

            # Set the number of images.
            self.nimages = nimages
            self.start = 0
            self.end = self.nimages - 1

            # Check which files are given as input - either the pathfile or the reactantfile and productfile:
            if pathfile is None:
                try:
                    assert(reactantfile is not None)
                except:
                   print("\nError in reactionpath.__init__: pathfile is None, but no reactantfile is given\n")
                try:
                    assert(productfile is not None)
                except:
                    print("\nError in reactionpath.__init__: pathfile is None, but no productfile is given\n")

                # Read the reactant structures.
                self.rpcxs = [cx.CXS(FileName=reactantfile)]

                # Fill in the intermediate images from 1 to nimage-2, by copying the reactant structure.
                # This is really just place-holding to make sure that the cxs objects have space so that the
                # interpolation scheme can fill them in...
                for i in range(1,nimages):
                    self.rpcxs.append(cx.CXS(FileName=reactantfile))

                # Read the product structure.
                self.rpcxs.append(cx.CXS(FileName=productfile))

                # At this point, align the structures using Kabsch algorithm, if requested.....
                if AlignEndPoints:
                    self.rpcxs[nimages].KabschAlign(self.rpcxs[0])

            # If pathfile is defined, read the full path from file.
            else:

                # Read the full path file.
                self.nimages,self.rpcxs = io.ReadPathFromFile(pathfile)

        # If requested, perform geometry optimization of the reactant and product.
        if OptimizeEnds:
            success = pescalc.Driver(self.rpcxs[0], EnergyCalc, type='Optimize')
            e1 = self.rpcxs[0].energy
            if not success:
                print("\n** Harsh Realm...optimization of REACTANTS failed in reactionpath.__init__\n")
                sys.exit(1)
            success = pescalc.Driver(self.rpcxs[self.nimages], EnergyCalc, type='Optimize')
            e2 = self.rpcxs[self.nimages].energy
            if not success:
                print("\n** Harsh Realm...optimization of PRODUCTS failed in reactionpath.__init__\n")
                sys.exit(1)


    def LinearInterpolation(self):
        """
        Generates interpolated atomic coordinates for the intermediate images in a reaction-path object.
        The end-points are not changed by this routine.
        """

        # Simple interpolation in Cartesian coordinates (not great for curvilinear reactions!):
        x0 = self.rpcxs[0].positions
        xn = self.rpcxs[self.nimages].positions
        i = 1
        for c in self.rpcxs:
            lam = float(i-1)/float(self.nimages)
            i += 1
            c.positions = self.rpcxs[0].positions + lam * (xn-x0)

        return


    def GetFullPathEnergy(self,EnergyCalcs):
        """
        Calculates the potential energy at every image along the reaction-path object ``RP``.

        Args:
              EnergyCalcs (PES, list): List of ``PES`` calculation objects.
        """
        for c in self.rpcxs:
            c.success = pescalc.Driver(c, EnergyCalcs, type='Energy')


    def ReportFullPathEnergy(self,File=None,comment=None):
        """
        Reports the energy at each point along a reaction-path. This is output to either File or stdout.

        Args:
            File (str): A file in which to output the path energy.
            comment (str): An optional comment describing the current path.
        """
        # Setup output destination.
        if File is None:
            out = sys.stdout
        else:
            out = File

        if comment is not None:
            out.write("\n# "+comment+" #\n")
        out.write("--------------------------------------------------------------------------------------------------\n")
        out.write("# Index  |  Arc length (au)  |  Energy (au)  |  Relative energy (au)  |  Relative energy (kJ/mol)\n")
        out.write("--------------------------------------------------------------------------------------------------\n")
        for i in range(self.nimages+1):
            if i == 0:
                ArcLength = 0.0000
                E0 = self.rpcxs[i].energy
            else:
                ArcLength += np.linalg.norm(self.rpcxs[i].positions - self.rpcxs[i-1].positions)
            dE = self.rpcxs[i].energy-E0
            out.write('{:6d} {:18f} {:16f} {:20f} {:23f} \n'.format(i,ArcLength,self.rpcxs[i].energy,
                                                                    dE, dE * const.hartrees_to_kjmol))
        out.write("--------------------------------------------------------------------------------------------------\n")


    def NEB(self, PotentialType=None, spring = 0.05, Iterations = 10, ForceTolerance = 1E-3,
            MaxForce = 5E-3, ClimbingImage = False, Fileroot = None, EnergyCalc = None, Alpha = 0.1,
            TimeStep = 20.0, Nmin = 5, File = None):
        """Performs a nudged-elastic band calculation on the current reaction path.

        This routine implements the climbing-image NEB method, using the fast inertial relaxation
        method (FIRE) as the optimizer.

        Args:
            PotentialType (str): Either 'IDPP' for an image-dependent pair potential interpolation, or 'PES' to
                                 perform full NEB optimization on the PES.
            spring (float): NEB spring constant (in atomic units)
            Iterations (int): Maximum number of NEB optimization steps.
            ForceTolerance (float): RMS force-tolerance to signify convergence (in atomic units)
            MaxForce (float): Maximum force tolerance to signify convergence (in atomic units). Set to a negative
                              value to just use ForceTolerance to judge convergence.
            ClimbingImage (bool): Flag to indicate whether to perform climbing-image variant of NEB.
            Fileroot (str): A file-root name, which is prefixed onto the output files from the NEB calculation. If
                            ``Fileroot is None``, then a set of locally-named files is created instead.
            EnergyCalc (``PES`` objects): Set of PES calculation objects, used if ``PotentialType = PES``.
            Alpha (float): Initial alpha parameter used in the FIRE scheme.
            TimeStep (float): Time-step (in atomic units) for the FIRE scheme.
            Nmin (int): Nmin parameter from the FIRE scheme.
            File (file object): An open file object, or stdout, where output is placed.

        """

        # Setup output destination.
        if File is None:
            out = sys.stdout
        else:
            out = File

        # Define the number of internal images, assuming the end-points are fixed.
        ninternal = self.nimages - 2

        # Define number of atoms.
        n = self.rpcxs[0].natoms

        # Open output files.
        if PotentialType == 'PES':
            suffix = "neb_pes"
        elif PotentialType == 'IDPP':
            suffix = "neb_idpp"
        fconv = open(Fileroot + "_"+suffix+".conv", "w")
        finitial = open(Fileroot+"_"+suffix+"_initial.xyz","w")

        out.write("\n----- NEB output files -----\n")
        out.write("\n* Convergence analysis file: "+str(fconv.name)+"\n")
        out.write("* Initial MEP file: " + str(finitial.name) + "\n")

        # Pretty-print some initial information..
        out.write("\n----- Initial NEB parameters -----\n")
        out.write("\n* Potential Type = "+PotentialType+"\n")
        out.write("* NEB spring constant = "+str(spring)+" Eh / a0**2 \n")
        out.write("* RMS force tolerance = " + str(ForceTolerance) + " Eh / a0 \n")
        out.write("* Maximum image force tolerance = " + str(MaxForce) + " Eh / a0 \n")
        out.write(
            "\n-----------------------------------------------------------------------------------------\n")
        out.write("# Iteration  |  <NEB Forces> (au)  |  <Perp. Forces> (au)  |  Max. Perp. Force (au)  \n")
        out.write(
            "-----------------------------------------------------------------------------------------\n")

        # Calculate initial IDPP or PES forces on all internal images.
        if PotentialType == 'IDPP':
            self.EvaluateIDPP()
        else:
            self.GetFullPathEnergy(EnergyCalc)

        # Calculate total forces, accounting for spring and projection.
        self.GetNEBForces(spring, ClimbingImage)

        # Zero initial velocities of all images.
        for c in self.rpcxs:
            c.velocities = np.zeros([n,3])

        # Output initial information.
        rmsforce, rmsperpforce, maxforce, maxperpforce, converged = self.CheckNEBConvergence(ForceTolerance, MaxForce)
        out.write('{} {:20f} {:20f} {:24f} \n'.format("  initial",rmsforce,rmsperpforce,maxforce))


        # Loop over iterations.
        for ineb in range(Iterations):

            # Quickmin projection.
            for i in range(1, self.nimages):

                TotalForceNorm = 0.0
                for j in range(self.rpcxs[i].natoms):
                    for k in range(3):
                        TotalForceNorm += self.rpcxs[i].NEBforces[j,k]**2
                bconv = np.sqrt(TotalForceNorm / (3 * self.rpcxs[i].natoms))
                TotalForceNorm = np.sqrt(TotalForceNorm)

                dot = 0.0
            #for i in range(1, self.nimages):
                for j in range(self.rpcxs[i].natoms):
                    for k in range(3):
                        dot += self.rpcxs[i].velocities[j,k]*(self.rpcxs[i].NEBforces[j,k]/TotalForceNorm)

                if dot < 0.0:
                #for i in range(1, self.nimages):
                    self.rpcxs[i].velocities[:,:] = 0.0
                else:
                #for i in range(1, self.nimages):
                    for j in range(self.rpcxs[i].natoms):
                        for k in range(3):
                            self.rpcxs[i].velocities[j, k] = dot * (self.rpcxs[i].NEBforces[j, k] / TotalForceNorm)


            # Update positions and velocities using Euler step.
            #for i in range(1, self.nimages):
                for j in range(self.rpcxs[i].natoms):
                    for k in range(3):
                        self.rpcxs[i].positions[j,k] += TimeStep * self.rpcxs[i].velocities[j,k]
                        self.rpcxs[i].velocities[j,k] += TimeStep * (self.rpcxs[i].NEBforces[j,k] / self.rpcxs[i].mass[j])



            # Evaluate new forces.
            if PotentialType == 'IDPP':
                self.EvaluateIDPP()
            else:
                self.GetFullPathEnergy(EnergyCalc)

            # Calculate total forces, accounting for spring and projection.
            self.GetNEBForces(spring, ClimbingImage)

            # Check convergence.
            rmsforce, rmsperpforce, maxforce, maxperpforce, converged = self.CheckNEBConvergence(ForceTolerance, MaxForce)
            out.write('{:3} {:26f} {:20f} {:24f} \n'.format(ineb+1,rmsforce,rmsperpforce,maxperpforce))

            if converged:
                out.write(
                    "-----------------------------------------------------------------------------------------\n")
                out.write("\n*** CONVERGENCE ACHIEVED AFTER "+str(ineb+1)+" ITERATIONS ***\n")
                return

        # End of iterations - notify if failed to converge.
        if converged != True:
            out.write(
                "-----------------------------------------------------------------------------------------\n")
            out.write("\n*** FAILED TO CONVERGE AFTER " + str(ineb + 1) + " ITERATIONS ***\n")
            out.write("\n*** PLEASE CHECK RESULTS CAREFULLY TO IDENTIFY THE PROBLEM! ***\n")
        return


    def PrintToFile(self, File = None, format="xyz" ):
        """Prints an entire reaction-path object ``rp`` as a series of xyz-format frames.

        Args:
            format (str): Define the file output extension, currently only "xyz" is allowed.
            File (str): File to print reaction path in.
        """

        fout = open(File, "w")
        for c in self.rpcxs:
            n = c.natoms
            fout.write(str(n)+"\n")
            fout.write("\n")
            for i in range(n):
                x = c.positions[i,0]*const.bohr_to_ang
                y = c.positions[i,1]*const.bohr_to_ang
                z = c.positions[i,2]*const.bohr_to_ang
                fout.write(c.atomlabels[i]+"  "+str(x)+"  "+str(y)+"  "+str(z)+"\n")
        fout.close()


    def GetNEBForces(self, spring = 0.05, ClimbingImage = False):
        """Calculates the NEB forces on all images.

        The forces from the PES or from IDPP are assumed to be present in self.rpcxs[:].dvdr, and this routine
        calculates the perpendicular forces plus the parallel spring contribution.

        The total forces [n,3] at each ``CXS`` object are returned in ``self.rpcxs[:].NEBforces[n,3]``.

        If ``ClimbingImage == True``, the climbing-image force is calculated on the highest-energy image.

        Notes:
            Climbing-image should not be used for IDPP, as the underlying effective potential is not related
            to the actual potential energy surface of the system.

        Args:
            spring (float): NEB spring constant.
            ClimbingImage (bool): A flag indicating whether to use climbing-image forces.

        """
        nim = self.nimages
        n = self.rpcxs[0].natoms
        for c in self.rpcxs:
            c.tangent = np.zeros([n,3])


        # First, calculate the tangent vectors at each image.
        self.GetTangents()

        # If using climbing-image, identify the highest-energy image.
        if ClimbingImage:
            En = self.rpcxs[:].energy
            max_value = max(En)
            CI_image = En.index(max_value)

        # Now we have the tangents, calculate the projected forces at each image.
        for i in range(1,nim):

            dot = 0.000
            for j in range(n):
                for k in range(3):
                    dot += self.rpcxs[i].dvdr[j,k] * self.rpcxs[i].tangent[j,k]
            self.rpcxs[i].PerpendicularForces = -(self.rpcxs[i].dvdr - dot * self.rpcxs[i].tangent)

            d1 = np.linalg.norm(self.rpcxs[i+1].positions - self.rpcxs[i].positions)
            d2 = np.linalg.norm(self.rpcxs[i].positions - self.rpcxs[i-1].positions)
            self.rpcxs[i].ParallelForces = spring*(d1-d2)*self.rpcxs[i].tangent
            #
            # dot = 0.000
            # for j in range(n):
            #     for k in range(3):
            #         a = (self.rpcxs[i+1].positions[j,k] - self.rpcxs[i].positions[j,k])
            #         b = (self.rpcxs[i].positions[j,k] - self.rpcxs[i-1].positions[j,k])
            #         dot += spring * (a-b) * self.rpcxs[i].tangent[j,k]
            # self.rpcxs[i].ParallelForces = dot * self.rpcxs[i].tangent

            self.rpcxs[i].NEBforces = self.rpcxs[i].PerpendicularForces + self.rpcxs[i].ParallelForces

        return


    def CheckNEBConvergence(self, ForceTolerance = 1E-3, MaxForce = 5E-3):
        """Calculates indicators of NEB convergence, and decides whether or not convergence
        has been achieved.

        Returns:
            rmsforce (float): RMS total force on the NEB string. This includes the perpendicular component
                              from the PES plus the parallel spring component.
            rmsperp (float): RMS of the perpendicular forces from the PES (in atomic units).
            maxforce (float): Maximum single force component for the total forces (in atomic units).
            maxperpforce (float): Maximum perpendicular force component (in atomic units)
            converged (bool): Flag indicating whether the string is converged. This is based on comparing
                        ``rmsperp`` and ``maxperpforce`` to the input thresholds ``ForceTolerance`` and ``MaxForce``.

        """

        rmsforce = 0.000
        rmsperpforce = 0.000
        maxforce = -1E6
        maxperpforce = -1E6
        n = self.rpcxs[0].natoms
        icount = 0
        for i in range(1,self.nimages):
            for j in range(n):
                for k in range(3):
                    icount += 1
                    rmsforce += self.rpcxs[i].dvdr[j,k]**2
                    if abs(self.rpcxs[i].dvdr[j,k]) > maxforce:
                        maxforce = abs(self.rpcxs[i].dvdr[j,k])
                    rmsperpforce += self.rpcxs[i].PerpendicularForces[j,k]**2
                    if abs(self.rpcxs[i].PerpendicularForces[j,k]) > maxperpforce:
                        maxperpforce = abs( self.rpcxs[i].PerpendicularForces[j,k] )
        rmsforce = np.sqrt( rmsforce / float(icount) )
        rmsperpforce = np.sqrt(rmsperpforce / float(icount) )

        converged = ( (rmsperpforce < ForceTolerance ) and (maxperpforce < MaxForce) )

        return rmsforce, rmsperpforce, maxforce, maxperpforce, converged


    def GetTangents(self):
        """Calculates the tangents to the NEB string at each image.

        Notes:
            The tangents are estimated using the improved tangent method described in:

            Improved tangent estimate in the nudged elastic band method for finding minimum energy paths
            and saddle points, G. Henkelman and H. Jónsson, J. Chem. Phys., 113, 9978 (2000).
            https://doi.org/10.1063/1.1323224

        """
        n = self.rpcxs[0].natoms

        for i in range(1,self.nimages):
            self.rpcxs[i].tangent = np.zeros([n,3])

            #self.rpcxs[i].tangent = self.rpcxs[i + 1].positions - self.rpcxs[i-1].positions

            E = self.rpcxs[i].energy
            Ep = self.rpcxs[i+1].energy
            Em = self.rpcxs[i-1].energy

            # Case 1: V(i+1) > V(i) > V(i-1).
            if (Ep > E) and (E > Em):
                self.rpcxs[i].tangent = self.rpcxs[i+1].positions - self.rpcxs[i].positions

            # Case 2: V(i+1) < V(i) < V(i-1).
            elif (Ep < E) and (E < Em):
                self.rpcxs[i].tangent = self.rpcxs[i].positions - self.rpcxs[i-1].positions

            # Case 3: image i is a minimum or maximum.
            else:
                isMin = ( (Ep > E) and (Em > E) )
                isMax = ( (Ep < E) and (Em < E) )
                if (isMin == True) or (isMax == True):
                    dvmax = max( abs(Ep-E), abs(Em-E) )
                    dvmin = min( abs(Ep-E), abs(Em-E) )
                    dp = self.rpcxs[i + 1].positions - self.rpcxs[i].positions
                    dm = self.rpcxs[i].positions - self.rpcxs[i-1].positions
                    if Ep > Em:
                        self.rpcxs[i].tangent = dp * dvmax + dm * dvmin
                    else:
                        self.rpcxs[i].tangent = dp * dvmin + dm * dvmax

            # Final job is to normalize the tangent
            norm = np.linalg.norm( self.rpcxs[i].tangent )
            self.rpcxs[i].tangent = self.rpcxs[i].tangent / norm

        return


    def EvaluateIDPP(self):
        """Calculates the IDPP function and derivatives at all images along the ``rp`` object.
        """

        # Define space for targets.
        n = self.rpcxs[0].natoms
        targets = np.zeros([n,n])

        # Calculate the pair distance matrix for reactants and products.
        Rdist = self.rpcxs[0].GetAllPairDistances()
        Pdist = self.rpcxs[self.nimages].GetAllPairDistances()

        # Loop over internal images.
        for i in range(1,self.nimages):
            lam = float(i) / float(self.nimages)
            targets = Rdist + lam * (Pdist-Rdist)
            self.rpcxs[i].IDPPFunction(targets)

        return