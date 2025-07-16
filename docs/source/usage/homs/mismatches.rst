.. include:: /defs.hrst

.. _mismatches_usage:

Creating and detecting mode mismatches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mode mismatches in a |Finesse| model arise during the beam tracing when a beam parameter
set at a node is not equal to the beam parameter at a predecessor node transformed by
the associated components ABCD matrix. See :ref:`transforming_beam_param` for the
transformation formula and a table of component ABCD matrices.

These mismatches occur as a result of:

- manually set beam parameters at nodes which are mismatched from each other or any
  cavities in the model when propagating them,
- and / or multiple cavities in the model which are not matched to each other.

The former case here can be generated in a convenient way using
:meth:`.Model.create_mismatch` and both cases can be detected with the methods
:meth:`.Model.mismatches_table` and :meth:`.Model.detect_mismatches`. Before these
methods are covered, however, it is important to define how |Finesse| computes a
figure-of-merit for a mismatch value.

Defining a figure of merit
``````````````````````````

The mismatch between two beam parameters, :math:`q_1` and :math:`q_2`, is computed with
:meth:`.BeamParam.mismatch` and is defined as,

.. math::

    \mathcal{M} = \frac{\left|q_1 - q_2\right|^2}{\left|q_1 - q_2^*\right|^2}.

This quantity gives a value :math:`\mathcal{M} \in [0, 1]` where 0 represents no
mismatch (i.e. perfectly matched beam parameters) and 1 represents total mismatch.

Creating a custom mismatch
``````````````````````````

The :class:`.Model` class in |Finesse| exposes a method :meth:`.Model.create_mismatch`
which allows you to "inject" an arbitrary mode mismatch at any node in a configuration.
As a beam parameter is described completely by the waist-size :math:`w_0` and distance
to waist :math:`z`, this method provides separate mismatch arguments for these two
parameters in terms of the percentage mismatch you want to inject.

An example of this is shown below where we define an initially mode-matched cavity model
and then inject a mode-mismatch at the first surface of the input mirror.

.. jupyter-execute::

    import finesse
    finesse.configure(plotting=True)

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-2.5
    s sCAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=2.5

    # create the cavity object
    cav FP ITM.p2.o

    pd trns ETM.p2.o

    # only care about even modes as we're just looking
    # at mismatches, so select even modes up to order 4
    modes(even, 4)
    """)

    out_matched = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")
    out_matched.plot(logy=True);

The above plot shows the cavity scan with no mismatches present - as there is just a
single Cavity object present with no manually set beam parameters anywhere. Now we
create a mismatch of 10% in :math:`w_0` and -8% in :math:`z` at ITM:

.. jupyter-execute::

    mm_gauss = model.create_mismatch(model.ITM.p1.i, w0_mm=10, z_mm=-8)

    out_mismatched = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")
    out_mismatched.plot(logy=True);

From this plot we can see that we now get scattering into the 02 and 04 modes as we
would expect from introducing a mode mismatch.

Note that the :meth:`.Model.create_mismatch` method returns the
:class:`~finesse.components.gauss.Gauss` object (i.e. user-set beam parameter object)
which was created / modified in the method. We will use the return value (`mm_gauss`) of
this example later on.

Detecting mode mismatches
`````````````````````````

As mentioned previously, there are two Model methods for detecting / displaying
mismatches - :meth:`.Model.detect_mismatches` and :meth:`.Model.mismatches_table`. The
former is the function which performs the actual computations whilst the latter can be
used to display a nicely formatted table of the mismatch values present.

.. rubric:: Detect mismatches method

The :meth:`.Model.detect_mismatches` method is very simple to use and only takes a
single optional argument which is a flag specifying whether to use the previous beam
trace results (stored in :attr:`.Model.last_trace`) or perform a new beam trace and
compute the mismatches from that. This method will return a dictionary of surface to
mismatch couplings where the mismatch couplings item is itself a dictionary of the
optical node couplings mapping to their mismatch values.

This is easier to visualise using the :meth:`.Model.mismatches_table` method but this
method is also exposed as a programmatic interface for accessing specific mismatches if
you require them.

.. rubric:: Print mismatches method

:meth:`.Model.mismatches_table` calls the above method and prints the resulting
dictionary in a tabulated format for easier viewing of the mismatches. Taking the
example above where we created mismatches in a cavity model, we can then print the
mismatches to verify them:

.. jupyter-execute::

    print(model.mismatches_table())

.. rubric:: Detecting coupling coefficients directly

Another way to detect mismatches indirectly is to probe the coupling coefficients of a
scattering matrix during a simulation run. This will provide you with complex data for
the specific mode coupling that you are probing, over any parameter scan you are
performing.

A simple example of this is shown below, where we take the cavity model from above and
scan over the waist size of the "gauss command" added by the
:meth:`.Model.create_mismatch` call performed earlier.

.. todo:: Example using knm detector whilst scanning mm_gauss.w0

Trace priority
``````````````

.. todo::

    Document how the ``priority`` arguments for cavities and gauss objects affect
    tracing, and how to use them.
