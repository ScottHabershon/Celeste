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

"""Container module for parameter access across Celeste.

Sets default parameters, and makes them available module-wide.

Attributes:
    IOvars (dict): Dictionary containing possible input keywords.
    PositiveIntegerInputs (str): List of keywords which should correspond to positive integers.
    IOvarsFreeForm (str): List of parameters whcih are free-form strings.
    IOvarsAllowed (dict): Dictionary of allowed keywords.

"""

# Initialize params.IOvars with some default values here.
IOvars = {}
IOvars['CalculationType'] = None
IOvars['Description'] = None
IOvars['ReactantFile'] = None
IOvars['ProductFile'] = None
IOvars['PathFile'] = None
IOvars['InterpolationType'] = "LinearCartesian"
IOvars['Nimages'] = 0
IOvars['PEScalc'] = None
IOvars['OptimizeEnds'] = False
IOvars['AlignEndPoints'] = False
IOvars['NEBspring'] = 0.050
IOvars['NEBiterations'] = 10
IOvars['NEBForceTolerance'] = 1E-3
IOvars['NEBMaxForce'] = 5E-3
IOvars['ClimbingImage'] = False

# Mechanism-finding.
IOvars['NSAiterations'] = 100
IOvars['SAInitialTemperature'] = 1E3
IOvars['MoveFile'] = None
IOvars['MaxReactions'] = 5
IOvars['ValenceRange'] = [['C', 1, 4],['O',1,2],['H',0,1],]
IOvars['ReactiveAtoms'] = ["all"]
IOvars['SAOutputFrequency'] = 500
IOvars['ErrorType'] = 'FullStructure'


# Identify input parameters which should be positive integers.
PositiveIntegerInputs = ['Nimages', 'NSAiterations', 'MaxReactions', 'Ninterpolate', 'NEBiterations',
                         'SAOutputFrequency']

# Identify input parameters which should be floats.
FloatInputs = ['NavigationSpring', 'SAInitialTemperature', 'NavigationStep', 'NEBspring',
               'NEBForceTolerance', 'NEBMaxForce']

# Identify which parameters are free-form strings.
IOvarsFreeForm = ["Description", "ReactantFile", "ProductFile", "PathFile", "MoveFile"]

# Identify allowed values of some input parameters
IOvarsFixedValues = ["CalculationType", "InterpolationType"]
IOvarsAllowed = {}
IOvarsAllowed['CalculationType'] = ['NEB', "Navigation", "SGDS", "DEGDS"]
IOvarsAllowed['InterpolationType'] = [None, "LinearCartesian", "IDPP"]

# The following are IO storage related to PES calculations.
IOpesinput = []
IOpesdict = []