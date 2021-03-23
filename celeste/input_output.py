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

"""General file I/O functions for celeste package.
"""

import sys
import getopt
import celeste.params as params

def ParseCommandLine(argv):
    """Parses the command-line arguments to retrieve the input and output files when
    using ``RunCeleste.py``.

    Args:
         argv (str) : Command line arguments for parsing.

    Returns:
        (tuple) containing

                - inputfile (str) : input file from command line
                - outputfile (str) : output file from command line
    """

    inputfile = ''
    outputfile = sys.stdout

    # Check input variables.
    usage1 = "\nUsage: RunCeleste.py -i <inputfile> -o <outputfile>\n"
    usage2 = "\nUsage: RunCeleste.py -i <inputfile> [OUTPUT WILL BE SENT TO SYS.STDOUT]\n"
    if len(sys.argv) != 5 and len(sys.argv) != 3:
        print("\nUsage cases:\n")
        print(usage1)
        print(usage2)
        sys.exit(1)
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])

    except getopt.GetoptError:
        print(usage)
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit(1)
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    return inputfile, outputfile


def ReadInputFile(inputfile):
    """Reads information from inputfile and stores it in the params.py module.

    Args:
        inputfile (str) : Reads the file in ``inputfile`` and parses input information based on keywords in
        ``celeste.params.py``.
    """

    # Check that the input file exists.
    try:
        f = open(inputfile)
    except IOError:
        print("Input file named " + inputfile + " does not exist...")
        sys.exit(1)

    # Read the input file.
    data = f.readlines()
    for line in data:

        # If line begins with #, then ignore it as a comment line.
        if line.startswith('#'):
            continue
        if line.startswith('\n'):
            continue

        # Read the input files - some entries have special formats, so they are checked first.
        # If the input string does not correspond to a special format entry, it is read as a
        # simple member of the dictionary params.IOvars.
        #
        key, value = line.split(":" , maxsplit=1)

        # If the key value is "PEScalc", then we need to store the dictionary which is
        # given in the input file.
        if key.strip() == 'PEScalc':
            params.IOpesinput.append(value.strip())

        # Check that key is in the input dictionary.
        if key not in params.IOvars:
            print("\n** Harsh Realm... " + key + " is not an allowed input parameter in params.IOvars\n")
            print("Allowed input keys are:")
            for key in params.IOvars:
                print(key)
            sys.exit(1)

        #Excluding the free-form parameters (e.g. filenames), sense-check input parameters.
        if key not in params.IOvarsFreeForm:

            if key in params.PositiveIntegerInputs:
                try:
                    int(value) < 0
                except ValueError:
                    print("\n** Harsh Realm..." + key + " must be greater than zero, currently "+value+"\n")
                    sys.exit(1)

            if key in params.IOvarsFixedValues:
                if (value.strip() in params.IOvarsAllowed[key]) == False:
                    print("\n** Harsh Realm..." + key + " must be one of: "+str( params.IOvarsAllowed[key]))
                    sys.exit(1)

        # Assign the key, value pair in the dictionary.
        params.IOvars[key.strip()] = value.strip()

        # Convert some of the input values to integers to avoid having to repeat this later.
        if key in params.PositiveIntegerInputs:
            params.IOvars[key.strip()] = int(value)

        # Convert some of the input values to integers to avoid repeat later.
        if key in params.FloatInputs:
            params.IOvars[key.strip()] = float(value)

    # Close the input file.
    f.close()


def StartupOutput(outputfile, inputfile = None):
    """Opens a generic logfile for calculation output and returns the corresponding filestrea as f.

    Args:
        outputfile (str) : File to be opened for output.
        inputfile (str) : Name of input file.
    """

    # Open output file.
    if outputfile == sys.stdout:
        outf = outputfile
    else:
        outf = open(outputfile, "w")

    # Store output file for use across other modules.
    params.OutputFile = outf

    str_line ="\n*============================================================================\n"

    outf.write(str_line)
    outf.write("* Celeste - Tools for automated reaction discovery\n")
    outf.write("\n* Main author: Scott Habershon\n\n")
    outf.write("* Contributing authors: Christopher Robertson (Warwick)\n")
    outf.write("                        Idil Ismail (Warwick)\n")
    outf.write("                        James Gardner (Warwick)\n")
    outf.write("                        Sebastian Havens (Warwick)\n")
    outf.write("\n* Address: Department of Chemistry, University of Warwick, Coventry, UK")
    outf.write(str_line)

    #TODO Need to add something about licence here.

    # Output the list of IOvars parameters.
    outf.write("*\n*** CALCULATION PARAMETER SUMMARY\n*")
    if inputfile is not None:
        outf.write("\n* Input filename: " + inputfile + "\n")
    for x, y in params.IOvars.items():
        if str(x) != 'PEScalc':
            outf.write("* " + str(x) + " :    " + str(y) + "\n")
    outf.write("*\n*** END OF CALCULATION PARAMETER SUMMARY\n*")

    # Note that we don't close the output file - it's needed for other output.