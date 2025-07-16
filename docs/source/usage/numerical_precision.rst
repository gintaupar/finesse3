.. include:: ../defs.hrst

.. _numerical_precision:

Numerical precision
===================

In elements
-----------

Normal integers and floats are represented in elements using Python data types. That
means integers have arbitrary precision and are only limited by system memory, while
floats are usually limited to 53 bits of precision, equivalent to about 16 or 17 digits
(a good discussion can be found `here
<https://docs.python.org/3/tutorial/floatingpoint.html>`_).

Numerical arrays in |Finesse| elements, e.g. for the computation of ABCD matrices, are
represented internally using Numpy data structures. Most integers are stored using
``int64`` representation. Larger numbers use ``uint64`` then default back to Python
``int`` objects.

The precision limits for integers and floats for your machine can be found in the
following way using Numpy:

.. jupyter-execute::

    import numpy as np

For integers:

.. jupyter-execute::

    print(np.iinfo(np.uint64))

For floats:

.. jupyter-execute::

    print(np.finfo(np.float64))

As an example, the following shows the switch from ``uint64`` to ``object`` in Numpy:

.. jupyter-execute::

    # Using the maximimum value, we get a uint64.
    print(repr(np.array(18446744073709551615)))

    # Using one more than the maximum value, we get an object.
    print(repr(np.array(18446744073709551616)))

**Only** numbers that can be represented by Numpy data types (either as scalars or
arrays) are fully supported in |Finesse|. |Finesse| needs access to certain Numpy
functions (e.g. :data:`numpy.cos`) and these are only implemented for Numpy data types.
*This effectively places a limit on the minimum and maximum numbers in |Finesse| to
those of Numpy.*

In simulations
--------------

.. todo:: describe numerical precision in |Finesse| simulations (KLU etc.).
