=======================
Input files for Celeste
=======================

When running ``Celeste`` as a stand-alone executable code, you'll need an input file.

The input file directs the code in deciding which calculation to perform, and defines the simulation parameters.

In addition, the input file can be used to define which energy calculation types should be used; this is explained in more detail in the [section on energy evaluations](energy.md)

Input file structure
--------------------

.. note::

    A complete description of all possible input keywords can be found :doc:`here <completeinput>`.

A typical ``Celeste`` input file is a plain-text file like the following:

.. code-block::

    # Example input file.

    CalculationType: NEB
    AlignEndPoints: True
    Nimages: 7
    ReactantFile: start.xyz
    ProductFile: end.xyz
    InterpolationType: IDPP
    OptimizeEnds: True
    NEBspring: 0.050
    NEBiterations: 100

    PEScalc: {'calctype' : "Energy", 'code' : "orca", 'executable' : "/Users/scott/programs/orca/orca", 'arglist' : ["! PM3 engrad"]}
  
    PEScalc: {'calctype' : "Optimize", 'code' : "orca", 'executable' : "/Users/scott/programs/orca/orca", 'arglist' : ["! PM3 Opt"]}


Note the following features of the input file:

* The first line, beginning `#`, shows a comment line. Comment lines can be added anywhere in the file, but must begin with `#`.

* Blank lines are ignored.

* Input directives are identifed using:

.. code-block::

      Keyword: value

* The list of available keywords is described in more detail [here.](completeinput)

* ``Celeste`` parses the input once initialized; useful errors should be raised if inconsistent data are input for a keyword.

* The ``PEScalc`` keywords can be used to define potential energy evaluation methods for energy calculation or geometry optimization.

.. important::
    The only way to define energy evaluations at the moment is through the ``PEScalc`` keyword.


Chemical structure formats
--------------------------

.. important::
    Currently, ``Celeste`` only supports input of chemical structures in ``xyz`` format.

In a standard ``xyz`` file, the structure is defined as a simple list of *(x,y,z)* atomic coordinates
(in Angstroms), along with the atom label, as follows:

.. code-block:: bash

    2
    This is a comment line
    C   0.0   0.0   0.0
    O   1.3   0.0   0.0

Note the following:

* The first line of the `xyz` file **MUST** contain the number of atoms

* The second line is a free-text comment. This can be used to store things like structural information, energies, or anything else. It can also be a blank line.

* Each subsequent line is of the format:

.. code-block:: bash

        AtomLabel  X  Y  Z

* The coordinates should be given in Angstrom units.


Input of multiple structure images
----------------------------------

When dealing with reaction-paths, we need to use files which tell us about the coordinates of all atoms at a series of
images (or frames) along the reaction path. Fortunately, ``xyz`` files can deal with this easily.

For example, a ``xyz`` file describing a "reaction path" simply contains a single ``xyz`` frame for each image,
concatenated together in a single file. So, a simple example containing three images is as follows:

.. code-block:: bash

    2
    This is the first image
    C   0.0   0.0   0.0
    O   1.3   0.0   0.0
    2
    This is the second image, slightly modified coordinates
    C   0.0   0.0   0.0
    O   1.4   0.0   0.0
    2
    This is the third image, slightly modified coordinates
    C   0.0   0.0   0.0
    O   1.4   0.0   0.0

