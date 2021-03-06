==================================
Defining allowed reactions (moves)
==================================

An important part of the graph-driven reaction-mechanism schemes encoded in ``Celeste`` is the
definition of *reaction moves*.

.. important::

    The input reaction moves are used to generate trial mechanisms; a mechanism is defined as a
    sequence of reaction moves, plus the indices of the atoms which the reaction is applied to.


When running ``Celeste`` (either as a standalone application from the terminal, or in an interactive python terminal),
reaction moves are input using a **MoveFile**; this is a plain text file which contains a series of blocks, with each block
defining an allowed reaction.


Format of Move File
-------------------

The **MoveFile** contains a set of input blocks, one for each input reaction move (or family).

The general format of each move block is as follows:

.. code-block::

    *move
    id: Association reaction
    namove: 2
    graph:
    0 1
    1 0
    *end


Note the following:

* Each move-definition begins with ``*move`` and ends with ``*end``.

* The ``id`` is a string describing this reaction - it is not used for anything other than identification.

* The integer ``namove`` defines the number of atoms involved in this reaction; in the example above, it is a two-atom association, so we have ``namove: 2``.

* The ``graph:`` directive indicates that the following lines containing the move-graph for the reaction. The move graph is a ``namove x namove`` matrix containing integers identifying the changes in bonding for atoms subject to the given reaction.

.. important::

    The ``graph`` contains integers -1, 0 or +1 which define whether bonds between the atoms remain unchanged (0), are formed (+1) or broken (-1).

    Some examples are as follows:


    - A two-atom move in which a bond forms between atoms *i* and *j*.

    .. code-block::

        graph:
        0 1
        1 0


    - A two-atom move in which a bond between atoms *i* and *j* breaks.

    .. code-block::

        graph:
        0 -1
        -1 0


    - A three-atom move in which a bond forms between atoms *i* and *j*, but the bond between *j* and *k* breaks.

    .. code-block::

        graph:
        0 1 0
        1 0 -1
        0 -1 0


    **Note:** The bonding graphs must be symmetric.

