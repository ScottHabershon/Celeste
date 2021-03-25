==============================
Setting up energy calculations
==============================

The ``Celeste`` package generally uses external codes to perform energy evaluations and geometry optimization.

The external energy calculation can be set-up in one of two ways:

    1. Through the ``PEScalc`` directive in the input file;
    2. Through the use of a template file, as specified in the input file;

.. important::
    Currently, only the ``PEScalc`` method for setting-up energy evaluations is employed, and the only current interface
    is to the `ORCA quantum chemistry code. <https://orcaforum.kofo.mpg.de/app.php/portal>`_.

Set-up using the ``PEScalc`` directive
--------------------------------------

The input file to ``Celeste`` can contain a series of ``PEScalc`` directives which are used to specify an
energy evalution method.

For example, a typical ``PEScalc`` directive in an input file is shown below:

.. code-block::

    PEScalc: {'calctype' : "Energy", 'code' : "orca", 'executable' : "/Users/scott/programs/orca/orca", 'arglist' : ["! PM3 engrad"]}


Note that the the input value given is in the form of a python dictionary, with each entry of the form:

.. code-block::

    key : value


The entries given in the dictionary are used to create a ``PES`` object, which contains information about the
type of calculation to be performed and the information required to run the calculation.

During a simulation, whenever an energy evaluation job is requested, the ``PES`` object is used to perform the corresponding
calculation on a provided ``CXS`` object, and the results are returned in the relevant entries of the ``CXS`` object
(such as ``CXS.energy`` for the total energy, or ``CXS.dvdr`` for the atomic derivatives).

Entries in the ``PEScalc`` directive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each ``PEScalc`` directive should contain the following dictionary entries:

    * ``calctype``: This defines the calculation type. Allowed values are:

        1. **Energy**: To perform a single-point energy and gradient evaluation.
        2. **Optimize**: To perform a geometry optimization to the nearest local minimum.

    * ``code``: Specifies the external code which will be used. Allowed values are:

        1. **orca**: To use the `ORCA quantum chemistry code. <https://orcaforum.kofo.mpg.de/app.php/portal>`_

    * ``executable``: To define the **full path** to the requested executable.

    .. warning::
        The relevant executable for energy evaluations obviously needs to be installed and available on your system!


    * ``arglist``: This is the list of arguments which will be passed to the energy evaluation code, in the form of a python list.
    In general, the format of this list is code-dependent; this information is used to create an input file for the
    target executable.

    In the case of using **ORCA** to perform external energy evaluation, the ``arglist`` should contain the
    usual input lines which would go in the header of an **ORCA** input file. For example, the ``arglist`` of:


    .. code-block::

        'arglist' : ["! PM3 engrad"]


    would result in production of an **ORCA** input file like this:


    .. code-block::

        ! PM3 engrad
        *xyzfile 0 1 temp.xyz


    Note that the ``! PM3 engrad`` has been input in the first-line of the created input file; the
    rest of the input file contains structure information for the ``CXS`` object on which a ``PES``
    calculation has been requested.


    .. important::
        The ``Celeste`` code handles the generation of the **ORCA** input file here; the user only needs to provide the
        calculation-type information for the header line.


    As a second example, consider the following ``PEScalc`` directive:


    .. code-block::

        PEScalc: { 'calctype' : "Energy", 'code' : "orca", 'executable' : "~/orca", 'arglist' : ["! DFT B3LYP def2-SVP", "! engrad"] }


    Note that, in this case, the ``arglist`` list contains two entries; each of these will become a separate
    line in the generated **ORCA** input file, as follows:


    .. code-block::

        ! DFT B3LYP def2-SVP
        ! engrad
        *xyzfile 0 1 temp.xyz


.. important::
    Note that any method which is accessible through **ORCA**, such as *PM3*, *AM1*, *DFT* and so on, can
    be accessed through ``Celeste.


.. warning::
    A ``Celeste`` input file can contain ``PEScalc`` directives with both ``calctype = Energy`` and ``calctype = Optimize``.
    However, if more than one ``PEScalc`` with ``calctype = Energy`` (or ``Optimize``) is defined, then this might cause internal problems - you should avoid this in your input files!


What is ``Celeste`` doing when it calculates energy?
----------------------------------------------------

For each dictionary (say, ``d``) created through defining a ``PEScalc`` directive, ``Celeste`` generates a
``PES`` object:

.. code-block:: python

    d={ 'calctype' : "Energy", 'code' : "orca", 'executable' : "~/orca",
        'arglist' : ["! DFT B3LYP def2-SVP", "! engrad"] }
    ecalc = PES(d)


During a ``Celeste`` calculation, suppose we generate a chemical structure (``CXS``) object ``a``; to perform PES calculation ``ecalc`` on the ``CXS`` object ``a``, the following command is issued:

.. code-block :: python3

    E = ecalc.RunPES( a )

When ``RunPES`` is executed, ``Celeste`` does the following general steps:

    1. Create a temporary input file for the external executable associated with the ``ecalc`` object (such as **ORCA**)

    2. Run the external code on the temporary input file using the ``executable`` defined in the ``PES`` object.

    3. Read the generated input back into the ``CXS`` object, including some indication of whether the
    external code failed.


Which external codes are supported?
-----------------------------------

Current external energy evaluation codes which are supported are:

    * **ORCA**

The following will be added at some point:

    * **DFTB+**
    * **PSI4**


Useful notes about energy calculations
--------------------------------------

* Once an external code fails to execute correctly (for example, non-convergence of SCF calculation), this will typically flag an error in ``Celeste`` (although, depending on context, the code might be able to struggle on).

* The temporary input files generated for external codes typically have the suffix ``temp``; so, you should avoid having similar filenames!

.. warning::
    ``Celeste`` will occasionally delete temporary files, such as ``temp.xyz``; make sure this does not
    result in your input files being deleted.

* The python modules ``subprocess`` is used to perform calculations using external codes, and is also used to remove temporary files from previous calculations.


Troubleshooting
---------------

If energy calculations fail...

* Have you defined the correct ``executable`` and ``code`` in the ``PEScalc``?

* Is there something wrong with the structure? Sometimes the last-generated output files from the external code will be saved, so can be viewed to get an idea of what's going on?

* Have you really got the external executable installed?

* If you're running on a HPC queue system, is the external executable installed on the compute nodes?