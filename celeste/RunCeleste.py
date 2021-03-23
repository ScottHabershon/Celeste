#!/Users/scott/opt/anaconda3/bin/python
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

""" Main driver routine to run the Celeste reaction discovery program.

    Celeste can be used to perform the following simulations:
        (1) Single-ended graph-driven discovery (GDS),
        (2) Double-ended mechanism search (DEGDS),
        (3) Reaction path analysis, particularly (CI)NEB.

    Example:
        To execute Celeste as a stand-alone code, you'll need an input file.
        You can find information and example input files in :ref:`input`.

        Once you have an input file, you can run Celeste code using:

            .. code-block:: bash

                    RunCeleste.py -i [input file] -o [output file]

"""

import sys
import celeste.input_output as io
import celeste.params as params
import celeste.reactionpath as rpath
import celeste.mechanism_finder as mfind
import celeste.pes as pes
import celeste.cxs as cxs

def main(argv):
    """Main driver routine for RunCeleste.py.

        Performs a variety of reaction-discovery and reaction-path finding simulations, as controlled through the
        ``input`` file.

        Args:
            argv (str) : List of input arguments from command line.

    """

    # Get the inputfile and output file from command line.
    inputfile, outputfile = io.ParseCommandLine(argv)

    # Read and check the input file.
    io.ReadInputFile(inputfile)

    # Startup the output file.
    io.StartupOutput(outputfile, inputfile)

    # Setup and check PES calculations using data from input file:
    EnergyCalcs = pes.InitializePESCalcs(npes = len(params.IOpesinput), PESdetails = params.IOpesinput )
    for en in EnergyCalcs:
        en.PrintPESDetails(params.OutputFile)


    # Decide which calculation type to run....
    # - MEP finder....
    if params.IOvars['CalculationType'] == "NEB":

        params.OutputFile.write("\n*** Starting NEB Calculation ***\n")

        # Initialize the reaction path; either read the whole thing, or interpolate.
        rp = rpath.RP( reactantfile = params.IOvars['ReactantFile'],
                        productfile = params.IOvars['ProductFile'],
                        pathfile = params.IOvars['PathFile'],
                        nimages = params.IOvars['Nimages'],
                        OptimizeEnds=params.IOvars['OptimizeEnds'],
                        EnergyCalc = EnergyCalcs,
                        AlignEndPoints = params.IOvars['AlignEndPoints'])

        # If the pathfile was empty, generate an initial interpolated reaction-path. This interpolation
        # in Cartesian Coordinates is not great as an initial guess, but can be improved by IDPP later.
        params.OutputFile.write("\n*** Generating initial linear interpolation...\n")
        rp.LinearInterpolation()

        # If required, performed a NEB calculation using IDPP:
        if params.IOvars['InterpolationType'] == "IDPP":
            rp.NEB(PotentialType='IDPP',
                    spring = params.IOvars['NEBspring'],
                    Iterations = 10 * params.IOvars['NEBiterations'],
                    ForceTolerance = params.IOvars['NEBForceTolerance'],
                    MaxForce = params.IOvars['NEBMaxForce'],
                    ClimbingImage = False,
                    Fileroot=inputfile,
                    File = params.OutputFile)

            rp.PrintToFile(inputfile+"_idpp_path.xyz",format="xyz")


        # Calculate the path energy:
        rp.GetFullPathEnergy(EnergyCalcs)


        # Report the path energy.
        rp.ReportFullPathEnergy(File=params.OutputFile,
                                comment="Reaction path energy after initial interpolation")


        # Run NEB analysis using PES.
        rp.NEB(PotentialType='PES',
               spring=params.IOvars['NEBspring'],
               Iterations= params.IOvars['NEBiterations'],
               ForceTolerance=params.IOvars['NEBForceTolerance'],
               MaxForce=params.IOvars['NEBMaxForce'],
               ClimbingImage=False,
               Fileroot=inputfile,
               File=params.OutputFile,
               EnergyCalc=EnergyCalcs)


    # - Double-ended graph-driven search (DEDGDS).
    elif params.IOvars['CalculationType'] == "DEGDS":

        # Initialize a mechanism object.
        mech = mfind.Mechanism( reactantfile = params.IOvars['ReactantFile'],
                                productfile = params.IOvars['ProductFile'],
                                maxreactions = params.IOvars['MaxReactions'],
                                valenceranges = params.IOvars['ValenceRanges'])

        # Initialize the reaction moves:
        #moveset = moves.ReadMoveSet( movefile = params.IOvars['MoveFile'] )

        # mech.RunMechanismFinder(moves = moveset,
                        #        nsaiterations = params.IOvars['NSAiterations'],
                        #        SAInitialTemperature = params.IOvars['SAInitialTemperature'])

# Main driver.
#
if __name__ == '__main__':
    main(sys.argv[1:])

