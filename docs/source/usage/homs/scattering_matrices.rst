.. include:: /defs.hrst

.. _arbitrary_scatter_matrices:

Computing arbitrary scattering matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The internal functions for computing coupling coefficient matrices
during a simulation are also exposed via specific methods in the
:mod:`finesse.knm` sub-module. The most-general function for computing
these matrices is :func:`.make_scatter_matrix` - this requires the matrix
type as the first argument, mapping to one of the specific scattering
matrix functions in this module.

As these Python-facing functions internally call the same C functions used
to compute coupling coefficients during simulations, you can rest assured
that these routines are highly optimised - allowing for fast computation of
these matrices.

Calculating coupling matrices via the Bayer-Helms formalism
```````````````````````````````````````````````````````````

To compute a scattering matrix based on a system with mode-mismatches and
tilt (angular) misalignments, one can use :func:`.make_bayerhelms_matrix` (or
:func:`.make_scatter_matrix` with "bayerhelms" as the type, as noted above).

This function requires the input and output beam parameters, in both tangential
and sagittal planes, as well as the misalignment angles in both planes. One can
optionally provide a medium refractive index (defaults to :math:`n_r = 1`) and /
or a beam wavelength (defaults to :math:`\lambda = 1064\,\mathrm{nm}`) too.

As an example, below, we compute a scattering matrix with only mismatches present,
up to a `maxtem` of 5 (see :ref:`selecting_modes`):

.. jupyter-execute::

    import finesse
    finesse.configure(plotting=True)

    from finesse.knm.tools import make_scatter_matrix

    # Non-astigmatic beam with 20% mismatch in w0
    # from input to output
    qx1 = qy1 = finesse.BeamParam(w0=1e-3, z=0)
    qx2 = qy2 = finesse.BeamParam(w0=1.2e-3, z=0)

    # No angular misalignments
    xgamma = ygamma = 0

    # Compute the scattering matrix, returning a KnmMatrix object
    kmat = make_scatter_matrix(
        "bayerhelms",
        qx1, qx2, qy1, qy2,
        xgamma, ygamma,
        maxtem=5,
    )

    # Plot the absolute value of the coupling
    # coefficients as a colormesh
    kmat.plot(cmap="bone");

As expected, we get a plot showing non-zero couplings for even mode orders. Note
that the ``kmat`` object here will be a :class:`.KnmMatrix` object, providing
a convenient interface for accessing specific coupling coefficients, printing the
matrix and plotting the matrix. The :meth:`.KnmMatrix.plot` method computes and
plots the absolute value of the coefficients by default; this can be altered as
shown in the documentation for this method.

This function also accepts mode selection keyword arguments, allowing for custom
mode indices to be used (see :ref:`selecting_modes` and :func:`.make_modes`). For
example, in this case we could select just even modes as we know that no angular
misalignments are present:

.. jupyter-execute::

    kmat = make_scatter_matrix(
        "bayerhelms",
        qx1, qx2, qy1, qy2,
        xgamma, ygamma,
        select="even", maxtem=8,
    )
    kmat.plot(cmap="bone");

.. rubric:: Accessing specific couplings

As implied previously, the return type of :func:`.make_scatter_matrix` is a
:class:`.KnmMatrix` object. Along with plotting and stringifying, this object
provides a dict-like interface for accessing coupling coefficients in the matrix.

Here are a few convenient ways of accessing these coefficients:

.. jupyter-execute::

    # Using the "nm->n'm'" format
    print("Coupling coeff. from 00 -> 02 = ", kmat["00->02"])

    # Using the "nmn'm'" format
    print("Coupling coeff. from 02 -> 40 = ", kmat["0240"])

    # Using the (n, m, n', m') format
    print("Coupling coeff. from 26 -> 44 = ", kmat[(2, 6, 4, 4)])
