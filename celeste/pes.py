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

"""Contains routines relating to potential energy surface calculations for
``CXS`` objects.
"""

import sys
import os
import celeste.chemical_data as chemdata
import celeste.constants as const
import numpy as np
import subprocess
import ast
import glob

class PES:
    """Defines the PES object, the main object used in interfacing with energy calculations.

    Args:
        calctype (str) : Calculation type. Allowed values:

                - *Energy*: Performs an energy and gradient evaluation.
                - *Optimize*: Performs geometry optimization.
                - *Frequency*: Performs frequency/thermodynamics calculation.

        code (str) : Identifies which code to run to perform calculation. Allowed values are **orca**
        executable (str): Identifies the external executable for ab initio calculation. *This should be an available
            executable installed on your computer.*
        calctemplate (str, optional) : Identifies template file for ab initio calculation.
        arglist (str,optional) : Gives directives to be used by external ab initio codes,
                for example relating to DFT calculation type or basis set.)

    Notes:
        PES objects can be initialized by either providing the ``arglist`` containing the
        set of arguments for the energy calculation, or by providing the name of a template file ``calctemplate``.

        If both ``arglist`` and ``calctemplate`` are provided to a ``PES`` object, it's not clear
        what the user intends (*i.e.* which one do we use?), so the code will flag an error.


    Example:
        To use the ``arglist`` to create a ``PES`` object to perform ORCA energy calcuations
        using DFT B3LYP and a def2-SVP basis:

            .. code-block:: python3

               d={ 'calctype' : "Energy",
                    'code' : "orca",
                    'executable' : "~/orca",
                    'arglist' : "! DFT B3LYP def2-SVP engrad" }
                ecalc = PES(d)

        Now suppose we have a ``CXS`` object created as follows:

            .. code-block:: python3

                a = CXS(natoms = 2, positions =[(0.0,0.0,0.0), (1.5,0.0,0.0)], atomlabels =['C','O'],
                comment = "Carbon monoxide")

        To perform PES calculation ``ecalc`` on the ``CXS`` object ``a``, we then do:

            .. code-block :: python3

                E = ecalc.RunPES( a )

        What happens when we run this command? First, the ``PES`` object writes a file which is in the correct
        format for the external *ab initio* code. Next, the external code is executed with the input file. Finally,
        the results of the calculation are read back into the ``CXS`` object. So, in the example above, the
        ORCA input file which was created looks like this:


             .. code-block :: console

                    ! DFT B3LYP def2-SVP engrad
                    *xyzfile 0 1 temp.xyz


        Here, the first directive line is the same as the one we had in the dictionary which initialized the
        ``PES`` object. The second line identifies the charge (0) and spin-multiplicity (1), as well as the
        temporary xyz file which stores the coordinates of the ``CXS`` structure; this temporary file (``temp.xyz``)
        is generated when ``ecalc.RunPES`` is executed.

        Note that the arglist entry in the input list can also be a list; in this case, each list member is
        interpreted as a new line in the input file required for the ab initio code. For example,

              .. code-block:: python3


                d={ 'calctype' : "Energy",
                    'code' : "orca",
                    'executable' : "~/orca",
                    'arglist' : ["! DFT B3LYP def2-SVP", "! engrad"] }
                ecalc = PES(d)

        The resulting ORCA input file in this case will look as follows:

               .. code-block :: console

                    ! DFT B3LYP def2-SVP
                    ! engrad
                    *xyzfile 0 1 temp.xyz

        Further examples of PES calculation set-ups can be found in the Examples section.


    """

    def __init__(self, calctype = "Energy", code = "orca", executable = "~/orca", arglist = None,
                 calctemplate = None, comment=None):
        """Creates a ``PES`` object for energy calculations.
        """
        self.calctype = calctype
        self.code = code
        self.executable = executable
        self.arglist = arglist
        self.calctemplate = calctemplate
        self.comment = comment
        self.fromtemplate = False
        self.fromarglist = False

        # Check the input - either from template file or arglist. If both are defined, we need to flag
        # an error....
        if self.calctemplate is None:   # No template file defined...

            self.fromarglist = True

            # Check that arglist is defined...if not, we stop.
            try:
                assert( self.arglist is not None )
            except:
                print("\nError in PES.__init__: Both arglist and calctemplate are empty!\n")
                sys.exit(1)

        else: # Template file found....

            self.fromtemplate = True

            # Check that arglist is empty...if not, we stop.
            try:
                assert( self.arglist is None )
            except:
                print("\nError in PES.__init__: both calctemplate are defined (not sure which to use)!\n")
                sys.exit(1)


    def PrintPESDetails(self, file = None):
        """Prints the details of the PES object to either file or stdout (if File is None)
        """
        if file is None:
            out = sys.stdout
        else:
            out = file
        out.write("\n*** PES Calculation details ***\n\n")
        out.write("- Calculation comment: "+str(self.comment)+"\n")
        out.write("- Calculation type: "+str(self.calctype)+"\n")
        out.write("- Calculation code: "+str(self.code)+"\n")
        out.write("- Executable: "+str(self.executable)+"\n")
        out.write("- Argument list: "+str(self.arglist)+"\n")
        out.write("- Template file: "+str(self.calctemplate)+"\n")
        out.write("\n*******************************\n")


    def RunPES(self, cx):
        """Performs a PES calculation on the chemical structure object ``cx``.

        Args:
            cx (``CXS`` object) : A chemical structure object provided by ``cxs.py``. This should be a valid ``CXS``
                                object containing atomic labels and coordinates.

        Returns:
            success (bool): Flag indicating whether calculation was successful.

        Notes:
            As a result of the electronic structure calculation, cx.energy, cx.dvdr and cx.positions will be
            updated, depending on calculation options:

                - *calctype = Energy*: ``cx.energy`` and ``cx.dvdr`` are read from calculation output, but
                        ``cx.positions`` remains untouched.
                - *calctype = Optimize": ``cx.energy``, ``cx.dvdr`` and ``cx.positions`` are all read from
                        calculation output.

        Notes:
            The return flag does not indicate anything about the reliability or quality of the results of the
            energy calculation; instead, it reports on simply whether the calculation ran as expected or not.
            It's up to you to make sure that the results are sensible!
        """

        if self.code == 'orca':
            success = self.RunOrcaCalc( cx )
        else:
            print("non-orca not yet implemented")
            sys.exit(1)
        cx.success = success
        return success


    def RunOrcaCalc(self, cx, FileRoot="orca_temp"):
        """Runs an ORCA calculation.

        This function writes an ORCA input file using either the argument list or template file.

        Args:
            cx (CXS object): Chemical structure object to run calculation on.
            FileRoot (str) : Root to use for temporary input/output files.

        Returns:
            success (bool): Flag indicating whether the calulation completed successfully.

        Notes:
            The ``success`` flag does not indicate anything about the energy calculation itself; it
            simply indicates whether or not the calculation finished, as judged by the appearance of
            the expected output files.

        """

        # Set the file root for the orca file.
        fr = FileRoot

        # Remove any previous orca data...
        ofiles = glob.glob('orca_temp*')
        subprocess.call(['rm'] + ofiles)

        # First, print the structure to xyz file.
        cx.PrinttoFile( file='temp.xyz',format = 'xyz' )

        # Calculate spin multiplicity...
        cx.GetSpinandCharge()

        # Option 1: Run from argument list.
        if self.fromarglist:

            f = open(fr+'.in',"w")
            for item in self.arglist:
                f.write(item+"\n")
            f.write("* xyzfile "+str(cx.TotalCharge)+" "+str(cx.multiplicity)+" "+"temp.xyz\n")
            f.close()

        # Option 2: Run from template file.
        else:
            print("* Not yet done from file in pes.RunOrcaCalc")

        # Execute the calculation.
        fout = open(fr+'.output',"w")
        subprocess.call( [str(self.executable),fr+".in"],stdout=fout )

        # Check that the engrad file exists, signifying calculation complete.
        if not os.path.exists(fr+".engrad"):
            success = False
            return success

        # Read in the lines from the engrad file.
        fen = fr + ".engrad"
        pattern1 = "# The current total energy in Eh"
        pattern2 = "# The current gradient in Eh/bohr"
        pattern3 = "# The atomic numbers and current coordinates in Bohr"
        with open(fen) as f:
            while True:
                line = f.readline().strip()
                if not line:
                    break

                # For the energy line, read the 2 lines after pattern1 - the energy is the second line...
                if line == pattern1:
                    data = [f.readline().strip() for i in range(2)]
                    cx.energy = float(data[1])

                # Read in derivative data in atomic units.
                if line == pattern2:
                    cx.dvdr = []
                    data = [f.readline().strip() for i in range(1+3*cx.natoms)]
                    for j in range(cx.natoms):
                        cx.dvdr.append([float(data[3*(j)+1]),float(data[3*(j)+2]),float(data[3*(j)+3])])
                    cx.dvdr = np.array(cx.dvdr)

                # If the calctype is optimization, then the structure itself needs to be updated by
                # reading in the new coordinates.
                if line == pattern3:
                    data = [f.readline().strip().split() for i in range(1 + cx.natoms)]
                    for j in range(cx.natoms):
                        cx.positions[j] = [float(data[(j)+1][1]),float(data[(j)+1][2]),float(data[(j)+1][3])]
                    cx.positions = np.array(cx.positions)

        # If we get here, everything (at least in reading) was fine...
        success = True
        return success


def InitializePESCalcs(npes = 0, PESdetails = None):
    """
    Initializes a set of PES calculation objects using PESdetails input dictionaries.
     
    Args:
        npes (int): Number of PES objects to be initialized (limited to 0, 1 or 2). 
        PESdetails (dict): Dictionary of keywords for PES set-up. These will be PES calculation dependent!

    Returns:
        EnergyCalcs (PES): List of PES calculation objects.
    """

    EnergyCalcs=[]
    for i in range(npes):
        d = ast.literal_eval(PESdetails[i])
        EnergyCalcs.append(PES(**d))

    # Check the details of the PES calculation types - we need to make sure there is only 1 Energy calculation
    # type and 1 optimization type.
    nopt = 0
    nen = 0
    for en in EnergyCalcs:
        s = en.calctype.strip()
        if s.lower() == 'optimize':
            nopt += 1
        if s.lower() == 'energy':
            nen += 1
    if nopt > 1:
        print("\n** Harsh Realm...Number of defined OPTIMIZE calculation objects must be less than 2!\n")
        sys.exit(1)
    if nen > 1:
        print("\n** Harsh Realm...Number of defined ENERGY calculation objects must be less than 2!\n")
        sys.exit(1)

    return EnergyCalcs


def Driver(cx,EnergyCalc,type='Energy',):
    """
    This is a PES evaluation driver. Given a set of energy calculation types and a ``CXS`` object, this
    routine identifies the correct energy calculation according to ``type``, then runs the calculation.

    This driver provides an alternative route to performing energy calculations which is useful in
    calculations where several mixtures of calculation types are requested (e.g. NEB calculations).

    Args:
        type: Type of calculation to be performed. Either 'Energy', 'Optimize' or 'Hessian'.
        cx: Chemical structure object to perform calculation on.
        EnergyCalc: List of energy calculation types.

    Returns:
        success (bool): Flag indicating problems (or not) with energy calculation.
    """
    OK = False
    for en in EnergyCalc:
        if en.calctype == type:
            OK = True
            success = en.RunPES(cx)

    if not OK:
        print("\n** Harsh Realm... " + type + " is not an available PES calculation type in pes.Driver\n")
        sys.exit(1)

    return success







