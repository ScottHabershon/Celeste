=====================
Algorithms in Celeste
=====================

``Celeste`` implements a series of different algorithms, several developed in the Habershon group at Warwick, and several developed elsewhere.
For clarity, this page summarizes the key algorithms.

Minimum-energy path finding
---------------------------

* **Interpolation using IDPP**

We use the standard image-dependent pair-potential (IDPP) algorithm, as described `here <https://doi.org/10.1063/1.4878664>`_,
to generate initial reaction-path approximations when ``InterpolationType`` is set to ``IDPP``.

* **Orientation of reactants and products**

When ``AlignEndPoints`` is ``True`` in the input file for a calculation, the product structure is rotated to give the
best-fit to the product structure. This is achieved using the `Kabsch algorithm <https://en.wikipedia.org/wiki/Kabsch_algorithm>`_, implemented using ``scipy.spatial.transform.Rotation.align_vectors``
to give the optimal rotation matrix.

.. note::
    This algorithm also shifts the structures so that their centroids coincide.

* **Nudged elastic band (NEB)**

To perform NEB simulations, we use the improved NEB approach described by `Henkelman <https://doi.org/10.1063/1.1323224>`_.

We also use the **QuickMin** algorithm to perform NEB optimization, as described `here <https://doi.org/10.1063/1.2841941>`_.




