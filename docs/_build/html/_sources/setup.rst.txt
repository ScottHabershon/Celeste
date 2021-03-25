==================
Setting up Celeste
==================

``Celeste`` is provided as a python package on `GitHub <https://github.com/ScottHabershon/Celeste>`_.

To install ``Celeste``, simply clone the ``Celeste`` package from `GitHub <https://github.com/ScottHabershon/Celeste>`_. This will create a `Celeste` directory with the following structure.

.. code-block::

    .
    └── Celeste
        ├── LICENSE.TXT
        ├── README.md
        ├── celeste
        ├── docs
        └── test

The directories contain the following:

* ``celeste`` contains the main python source-code.
* ``docs`` contains the documentation (you're reading it!)
* ``test`` contains a simple test example used during development.

Python path
-----------

If you want to ``import celeste``, you should add the path to the ``celeste`` package to your ``PYTHON_PATH``
environment variable.

Dependencies
------------

``Celeste`` uses several external python packages when performing calculations. You need to make sure that these
external packages are installed (and accessible in your ``PYTHON_PATH``).

Required packages are:

* `numpy <https://numpy.org>`_
* `scipy <https://www.scipy.org>`_
