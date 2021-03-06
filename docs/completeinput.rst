===================
Input specification
===================

``Celeste`` input files contain directives telling the code which calculations to perform and how.
This page contains a list of available input directives; further information can be found in the :doc:`API documentation <source/modules>`.

.. important::

    Note that the input of directives relating to energy evaluations is discussed :doc:`here <energy>`.

A note on units
---------------

.. important::
    ``Celeste`` uses atomic units, where the mass of the electron
    is 1.0, the charge on an electron is -1.0, :math:`\hbar=1` and :math:`4\pi\epsilon=1`.

    The exception to this is the input of atomic coordinates which, according to the usual convention of the
    ``xyz`` files used as input, are in Angstroms.


List of available input directives
----------------------------------

Core input directives
^^^^^^^^^^^^^^^^^^^^^

* ``CalculationType`` (**default**: None).

Defines the calculation to be performed. Allowed values are:

    1. ``NEB``: nudged elastic band refinement of reaction-path. Also includes climbing-image NEB.
    2. ``DEGDS``: double-ended reaction-mechanism search.
    3. ``SGDS``: single-ended reaction-network generation.

.. important::
    Currently, only ``NEB`` is supported as calculation type.

* ``Description``(**default**: None)

An optional string describing the calculation.

* ``ReactantFile`` (**default**: None)

An ``xyz`` file describing the reactant structure. This is used in initializing ``NEB``, ``DEGDS`` and ``SGDS`` calculations.

* ``ProductFile`` (**default**: None)

An ``xyz`` file describing the reactant structure. This is used in initializing ``NEB``, ``DEGDS`` and ``SGDS`` calculations.

* ``PathFile`` (**default**: None)
An ``xyz`` file containing a full reaction path, including reactants, products and intermediate images.

.. important::
    You should only define ``ReactantFile`` and ``ProductFile`` **or** ``PathFile``. If you provide
    all three at once, the code doesn't know exactly where to extract structure information from and will
    quit!


Minimum-energy path directives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``InterpolationType`` (**default**: *LinearCartesian*)

When performing a MEP calculation, this directive defines the method to be used to generate the internal images
by interpolation between the reaction end-points.

Available options are:

    1. ``LinearCartesian``: Performs linear interpolation between the end-points in Cartesian coordinates.
    This is usually not a good idea, especially for reactions which involve significant molecular distortions.

    2. ``IDPP``: Performs image-dependent pair-potential (IDPP) interpolation, as proposed by `Smidstrup et al <https://doi.org/10.1063/1.4878664>`_.

.. important::
    Note that ``InterpolationType`` is only relevant if ``PathFile`` is NOT defined. If ``PathFile`` is defined,
    then the MEP calculation will start from the provided reaction-path without interpolation.

* ``Nimages`` (**default**: 0)

This is the number of images in the reaction-path.

.. important::
    Numbering of images starts at 0, and runs to ``Nimages`` internally in ``Celeste``. This means that the total
    number of images in the reaction-path, including end-points, is ``Nimages + 1``.

* ``OptimizeEnds`` (**default**: *False*)

This is a logical flag which dictates whether the reactant and product structures are geometry-optimized
before the MEP calculation proceeds.

* ``AlignEndPoints`` (**default**: *False*)

This is a logical flag which dictates whether the reactant and product structures are orientationally
aligned (by the Kabsch algorithm) before path-interpolation and MEP.

* ``NEBspring`` (**default**: 0.050 Eh/Bohr)

This is the traditional NEB spring-constant, in atomic units.

* ``NEBiterations`` (**default**: 10)

This is the **maximum** number of iterations that a NEB calculation will run for; if convergence is signalled
before this maximum number of iterations, the calculation will stop.

* ``NEBForceTolerance`` (**default**: 1E-3 Eh/Bohr)

This is the RMS tolerance on the **perpendicular forces** along the reaction-path demanded to signal convergence.
This convergence criteria is based on the projected perpendicular forces along the reaction-path, without the NEB spring forces.

* ``NEBMaxForce`` (**default**: 1E-3 Eh/Bohr)

This is the target tolerance on the **maximum perpendicular force-component** along the reaction-path, calculated
without the NEB spring forces.

.. note::
    The ``NEBMaxForce`` value should typically be a bit larger than the ``NEBForceTolerance`` value, otherwise
    you'll be waiting a long-time for convergence....

* ``ClimbingImage`` (**default**: *False*)

This is a logical flag dictating whether or not to use the ``ClimbingImage`` modification of NEB.

.. warning::
    Currently, CINEB is not implemented.

Mechanism-finding directives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``MoveFile`` (**default**: *None*)

Defines the file which stores the available set of chemical reaction families (or *moves*) which can
be used during mechanism-finding or reaction-network generation simulations.

* ``MaxReactions`` (**default**: 1)

Maximum number of reaction-steps which will be used in mechanism-finding (``DEGDS``) or reaction-tree (``SGDS``) generation.

* ``ValenceRange`` (**default**: [['C', 0, 4],['O',1,2],['H',0,1],['N',1,3]])

This is a python tuple which gives minimum and maximum valence ranges for elements in the system. These valence ranges are used during mechanism-generation calculations to constrain the search; molecular structurs which contain atomic valences lying outside this range will be rejected.

For example, to set the carbon atomic valence range to [0,4] and the oxygen atomic valence range to [1,2], we would input:

.. code-block::

    ValenceRange: [['C', 0, 4],['O',1,2],['H',0,1],['N',1,3]]


* ``SAInitialTemperature`` (**default**: 1000)

Initial temperature (in K) for the simulated-annealing search for reaction-mechanisms connecting
input reactant and product configurations.

* ``NSAiterations`` (**default**: 1000)

Maximum number of simulated-annealing iterations which will be performed during search for
double-ended reaction mechanism. The search will stop once a mechanism with graph-error function value
less than ``constants.GCONV`` is found.

* ``SAOutputFrequency`` (**default**: 500)

Frequency with which output is produced during simulated-annealing mechanism search.

* ``ReactiveAtoms`` (**default**: ["all"])

A list of strings defining which atom-types are allowed to participate in reactions during mechanism-searching. For example:

.. code-block::

    ReactiveAtoms: ["C", "O"]


.. warning::
    The ``ReactiveAtoms`` list **must** be defined as a JSON list, so you must use **"all"**, rather than **'all'**.


* ``ErrorType`` (**default**: 'FullStructure')

Defines the graph error function-type to be used during optimization. Options are ``FullStructure`` for comparison based on the eigenvalues of the mass-weighted adjacency matrix for the full structure, or ``SingleMolecule`` for comparison of the eigenvalues of the mass-weighted adjacency matrix of a single molecule, namely that which most closely matches the molecule in the target product.

.. important::

    In the case of ``SingleMolecule`` as the graph error function, you should make sure that the product structure only contains a single target molecular product.


PES evaluation directives
^^^^^^^^^^^^^^^^^^^^^^^^^

* ``PEScalc`` (**default**: None)

This directive provides a definition of a potential energy calculation method; usually this is an interface
to an external *ab initio* electronic structure code, which can be used as a driver to perform energy
evaluations or geometry optimization.

The setup of ``PEScalc`` objects is further detailed in :doc:`here <energy>`.

