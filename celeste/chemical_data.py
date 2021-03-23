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

"""Contains dictionaries and other data structures holding chemical data such as
atomic masses and covalent radii.

Notes:
    All values are given in Hartree atomic units, where the mass of the electron
    is 1, the charge on an electron is -1, :math:`\\hbar = 1` and :math:`4 \\pi \\epsilon = 1`.

Attributes:
    nelements (int) : Total number of elements with data entries.
    ElementSymbol (dict) : Dictionary mapping element symbols to atomic numbers.
    AtomicMass (dict) : Dictionary mapping element symbol to atomic masses.
    CovalentRadii (dict) : Dictionary mapping element symbol to covalent radii.

"""

import sys

# Set the atomic label dictionary.
nelements = 100
ElementSymbol = {'H': 1,
                 'C':6,
                 'N':7,
                 'O':8,
                 }

# Set the atomic mass dictionary.
#
# See: IUPAC Technical Report, Pure Appl. Chem. 2016; 88(3): 265–291.
# Atomic weights taken from Tables 2/3. To convert from Daltons to Hartree atomic units,
#
# We use:
# 1 Da = 1 unified atomic mass unit = 1.660 5402 ×10−27 kg
# 1 electron mass = 9.109 383 7015 x 10-31 kg (NIST).
# So, to convert from atomic mass units to Hartree atomic units, multiply by
# (1 Da) / (1 electron mass) = 1822.889731.
#
AtomicMass = { 'H': 1837.47285,
               'C': 21894.72856,
               'N': 25533.21646,
               'O': 29164.41281,
               'Fe': 101799.277,
               }

# Set the atomic covalent radii in atomic units (Bohr).
#
# See: Covalent radii revisited, Cordero et al, Dalton Trans., 2008, 2832-2838
#
CovalentRadii = { 'H' : 0.59,  # 0.31 Angstroms
                  'C' : 1.44,  # 0.76 Angstroms
                  'N' : 1.34,  # 0.71 Angstroms
                  'O' : 1.25,  # 0.66 Angstroms
                  'Fe' : 2.87,  # 1.52 Angstroms, Take the larger high-spin value to be safe.
                }


def PrintChemicalDataReport(file=None):
    """Prints a report on the Chemical Data to file.
    Args:
        file (str) : File in which to print output. If ``file is None``, print to stdout.
    """

    if file is None:
        out = sys.stdout
    else:
        out = file

    # Mass reports.
    out.write("\n* Atomic mass dictionary (in Hartree atomic units):\n")
    for x, y in AtomicMass.items():
        out.write("Element "+x+", mass = "+str(y)+"\n")

    # Atomic covalent radii reports.
    out.write("\n* Atomic covalent radii dictionary (in Hartree atomic units [Bohr]):\n")
    for x, y in CovalentRadii.items():
        out.write("Element "+x+", Radius = "+str(y)+"\n")