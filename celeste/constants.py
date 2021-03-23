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

"""Contains useful numerical constants such as bonding scale factor, threshold parameters and more.

Attributes:
    BondingScaleFactor (float) : numerical factor used to multiply sum of covalent radii in
        determining cut-off distance for calculating adjacency matrix.
    BIG_BONDING (float) : An arbitrary large distance given to pairs of atoms which are not
        accessible by walking along bonds. Used during calculation of molecules from a graph.
    NMOLMAX (int) : Maximum recommended treatable number of molecules.
    NAMOLMAX (int) : Maximum recommended number of atoms in molecules.
"""

import sys

# Conversion from Bohr to Angstroms - multiply a number in Bohr by bohr_to_ang to get Angstroms.
bohr_to_ang = 0.52917721080

# Conversion from Angstroms to Bohr - multiply a number in Angstroms by ang_to_bohr to get Bohr.
ang_to_bohr = 1.88972612499359

# Conversion from Hartrees to kJ/mol - multiply an energy in Hartrees by hartrees_to_kjmol to get kJ/mol
hartrees_to_kjmol = 2625.4996300

# BondingScaleFactor influences the calculation of adjacency matrices. Two atoms are considered
# to be bonded if the distance between them is (CovalentRadius_1 + CovalentRadius_2)*BondingScaleFactor.
#
BondingScaleFactor = 1.10

# This is a BIG number, used to make sure that the Floyd-Warshel shortest-path algorithm in cxs.py
# avoids finding paths between non-bonded molecules.
BIG_BONDING = 10000

# NMOLMAX is the maximum number of independent molecules we're happy to deal with.
NMOLMAX = 50

# NAMOLMAX is the maximum number of atoms we're happy to have in each molecule.
NAMOLMAX = 25


def PrintConstants(file=None):
    """Prints a report on the current values of important parameters defined in ``constants.py``

    Args:
        file (str) : File in which to print output. If ``file is None``, print to stdout.
    """
    if file is None:
        out = sys.stdout()
    else:
        out = file

    out.write("\n* Useful numerial constants in constants.py *\n")
    out.write("- BondingScaleFactor = "+str(BondingScaleFactor)+"\n")
    out.write("- BIG_BONDING = "+str(BIG_BONDING)+"\n")

