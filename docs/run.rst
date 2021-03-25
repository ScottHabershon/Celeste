===============
Running Celeste
===============

``Celeste`` can be used in two different ways:

1. As a stand-alone executable code;
2. As an imported python package.

Running as an executable
------------------------

To run ``Celeste`` as a stand-alone executable code, type the following on the command line:

.. code-block:: bash

    RunCeleste.py -i [input file] -o [output file]

.. important::

    The above assumes that ``RunCeleste.py`` is in your ``PATH`` variable - you might need to change the above
    command to match your system!


The format of the input file is discussed in :doc:`input files section <input>`.

Importing as a python package
-----------------------------

To import ``Celeste``, you can do the usual:

.. code-block:: python

    import celeste

The :doc:`API documentation <source/modules>` gives much more information on available modules,
classes and functions.