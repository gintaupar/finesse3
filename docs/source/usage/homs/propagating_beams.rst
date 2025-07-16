.. include:: /defs.hrst

.. _propagating_beams:

Beam propagation
~~~~~~~~~~~~~~~~

|Finesse| provides a tool-set of beam tracing and mode matching functions which can be
used outside of a simulation context. These functions allow for arbitrary beam
propagation over any optical path of a given configuration - supporting both numeric and
symbolic calculations.

The primary function recommended for most beam propagation analyses is
:func:`~finesse.tracing.tools.propagate_beam`, also callable from a model via
:meth:`.Model.propagate_beam`.

The 'propagate_beam' function
`````````````````````````````

As noted in the API documentation linked above, the :func:`~finesse.tracing.tools.propagate_beam`
function can take any optical path
of a defined model and propagate an arbitrary beam from the start node to the end node
of the path.

Here we will use a proposed, preliminary design :cite:`PhysRevD.103.023004` of the
Einstein Telescope Low-Frequency (ET-LF) signal recycling cavity (SRC), complete with
arm telescopes, to highlight how this function can be used. The model is defined below
via |Finesse| kat-file syntax (see :ref:`syntax`).

.. jupyter-execute::

    import finesse
    finesse.configure(plotting=True)

    model = finesse.Model()
    model.parse("""
    ### L0 -> BS -> YARM of ET-LF
    # input
    l L0 P=1
    s l0 L0.p1 BS.p1 L=10

    # Main beam splitter
    bs BS T=0.5 L=37.5u alpha=60
    s BSsub1 BS.p3 BSAR1.p1 L=0.07478 nr=nsilica
    s BSsub2 BS.p4 BSAR2.p1 L=0.07478 nr=nsilica
    bs BSAR1 R=50u L=0 alpha=-36.6847
    bs BSAR2 R=50u L=0 alpha=36.6847

    # Y arm telescope
    s lBS_ZM1 BS.p2 ZM1.p1 L=70
    bs ZM1 T=250u L=37.5u Rc=-50
    s lZM1_ZM2 ZM1.p2 ZM2.p1 L=50
    bs ZM2 T=0 L=37.5u Rc=-82.5
    s lZM2_ITMlens ZM2.p2 ITM_lens.p1 L=52.5

    lens ITM_lens 75
    s lITM_th2 ITM_lens.p2 ITMAR.p1 L=0

    # Y arm input mirror
    m ITMAR R=0 L=20u
    s ITMsub ITMAR.p2 ITM.p1 L=0.2 nr=nsilicon
    m ITM T=7000u L=37.5u Rc=-5580

    # Y arm length
    s l_arm ITM.p2 ETM.p1 L=10k

    # Y arm end mirror
    m ETM T=6u L=37.5u Rc=5580
    s ETMsub ETM.p2 ETMAR.p1 L=0.2 nr=nsilicon
    m ETMAR R=0 L=500u

    # SRM
    s lBS_SRM BSAR2.p3 SRM.p1 L=10

    m SRM T=0.2 L=0 Rc=-9410
    s SRMsub SRM.p2 SRMAR.p1 L=0.0749 nr=nsilicon
    m SRMAR R=0 L=50n

    # cavities
    cav cavARM ITM.p2
    cav cavSRC SRM.p1 ITM.p1.i

    var nsilica 1.44963098985906
    var nsilicon 3.42009

    lambda(1550n)
    """)

We can call :func:`~finesse.tracing.tools.propagate_beam` using any two optical nodes as
the targets, in this case we propagate the beam from the ITM to the SRM.

.. jupyter-execute::

    ps = model.propagate_beam(model.ITM.p1.o, model.SRM.p1.i)

Note that we can use :class:`.Port` objects as the end points too, in which case the
`from_node` argument is deduced to be the output optical node of the port and `to_node`
will be the input optical node of the other port. In our example, the line below is
equivalent to the above.

.. jupyter-execute::

    ps = model.propagate_beam(model.ITM.p1, model.SRM.p1)

By specifying the two end nodes / ports, the optical path between these points will be
determined in the function. This is generally a fast operation, however it is likely the
slowest of all operations carried out by the
:func:`~finesse.tracing.tools.propagate_beam` function as a whole. Therefore, in the
rare case where this function is being called in a large loop, it makes sense to
pre-compute the optical path and pass this instead. An example of this is shown below,
again equivalent to the above.

.. jupyter-execute::

    ITM_TO_SRM = model.path(model.ITM.p1, model.SRM.p1)
    ps = model.propagate_beam(path=ITM_TO_SRM)

One thing you may note from this type of call to
:func:`~finesse.tracing.tools.propagate_beam` is that no input beam parameter (`q_in`)
argument has been specified. This means that this argument will be automatically deduced
from the model via an internal :meth:`.Model.beam_trace` call - where the beam parameter
at the input node of the path is then accessed and used as the `q_in` argument. A custom
beam parameter can be used by specifying a `q_in` argument explicity (this can be a
complex number or a :class:`.BeamParam` instance).

Another aspect to take into account is that
:func:`~finesse.tracing.tools.propagate_beam` operates on a single plane (i.e. 'x' for
the tangential plane, 'y' for the sagittal plane). By default, the tangential plane is
used. This can be changed to sagittal by specifying ``direction='y'`` as an argument.

The return value of :func:`~finesse.tracing.tools.propagate_beam` is a
:class:`.PropagationSolution` instance. We will look at the various properties and
methods this class provides below, in the context of our example.

.. rubric:: Plotting properties of the beam propagation

One of the most common use-cases is to plot a beam trace, in order to see how, in
particular, the beam size and accumulated Gouy phase evolve over a path. This is as
simple as calling :meth:`.PropagationSolution.plot` on the resulting solution. An
example is given below for the call in the above section.

.. jupyter-execute::

    ps.plot(
        name_xoffsets={"ITM_lens": 10, "SRM": 10},
        name_yoffsets={"ITM_lens": 10},
    );

The `name_xoffsets` and `name_yoffsets` arguments are optional and provide a way of
shifting the positions of the component names on the figure (in terms of data
coordinates) to avoid clashes.

Just the beam-sizes, for example, can also be plotted.

.. jupyter-execute::

    ps.plot_beamsizes(
        name_xoffsets={"ITM_lens": 10, "SRM": 10},
        name_yoffsets={"ITM_lens": 10},
    );

The wavefront curvature is not plotted by default with a call to
:meth:`.PropagationSolution.plot` but this can be plotted separately too with
:meth:`.PropagationSolution.plot_curvatures` or along with the others via giving
``"all"`` as the first argument to this method.

.. rubric:: Printing propagation data

As well as plotting, one can also display useful data associated with a beam propagation
by simply printing the solution object. For our example here, this results in:

.. jupyter-execute::

    print(ps)

A matrix of distances between all the optics in the propagation can also be printed
using :meth:`.PropagationSolution.distances_matrix_table`:

.. jupyter-execute::

    print(ps.distances_matrix_table())

.. rubric:: Accessing beam properties

:class:`.PropagationSolution` provides various methods for accessing physical properties
of the beam at any point along the computed path. All of these relevant methods
derive from :meth:`.PropagationSolution.q` (i.e. the beam parameter), so here we will
explore the options available for accessing ``q`` at different locations.

The simplest case is to provide an :class:`.OpticalNode` instance to get the beam
parameter at this node:

.. jupyter-execute::

    print(ps.q(model.BS.p2.i))

in the example above we grab the beam parameter at the input node of the second port of
the beam splitter (coming from the ZM1 optic). Whilst ``model.BS.p2.o`` (i.e. the output
node of this port) is technically not traced in the given path, we can still access the
beam parameter at this node as the beam tracing assumes that opposite node beam
parameters are the reverse of the forward node (i.e. :math:`-q^*`):

.. jupyter-execute::

    print(ps.q(model.BS.p2.o))

One can also use the string representation of an optical node:

.. jupyter-execute::

    print(ps.q("BS.p2.i"))

.. note::

    Any beam parameter value returned by :meth:`.PropagationSolution.q` will be a
    :class:`.BeamParam` instance, meaning all the various properties of this can be
    accessed as usual. :class:`.PropagationSolution` does provide some shortcuts for the
    key parameters (such as beam-size with :meth:`.PropagationSolution.w`), however, so
    the choice is up to the user as to whether to access like this:

    .. jupyter-execute::

        ps.q(model.BS.p2.i).w

    or like this:

    .. jupyter-execute::

        ps.w(model.BS.p2.i)

.. rubric:: Obtaining cumulative Gouy phases

It can be useful to see accumulated Gouy phases over segments of the traced path too,
and :class:`.PropagationSolution` provides methods to obtain these. These phases are
always given in degrees.

To obtain the total accumulated Gouy phase over the full traced path one can simply do:

.. jupyter-execute::

    ps.total_acc_gouy

Accumulated Gouy phase over a specific sequence of spaces can be retrieved with
:meth:`.PropagationSolution.acc_gouy`, e.g:

.. jupyter-execute::

    ps.acc_gouy("lZM2_ITMlens", "lZM1_ZM2")

gives the Gouy phase accumulated from the ITM lens to ZM1 (or vice-versa). And finally,
the Gouy phase accumulated up to a specific point can be obtained with
:meth:`.PropagationSolution.acc_gouy_up_to`. This gets the cumulative Gouy phase from
the starting node of the propagation up to the specified point. Similarly to accessing
beam properties, this point can be an :class:`.OpticalNode` instance or a component or
name of a component. For example, the below retrieves the accumulated Gouy phase from
the ITM to BS:

.. jupyter-execute::

    ps.acc_gouy_up_to("BS")


Symbolic beam propagation
`````````````````````````

The default behaviour of :func:`~finesse.tracing.tools.propagate_beam` is to use the
current values of each parameter, returning non-symbolic data in the
:class:`.PropagationSolution` instance. Another way to use this function is to switch on
the `symbolic` flag - resulting in symbolic expressions for all beam parameters, Gouy
phases and ABCD matrices at each point in the solution.

Using the example file from above again, we can simply call:

.. jupyter-execute::

    ps_sym = model.propagate_beam(path=ITM_TO_SRM, symbolic=True)

to compute a symbolic propagated beam solution. Accessing beam properties (such as the
beam parameter via :meth:`.PropagationSolution.q`) will now return a symbolic expression
rather than just a numeric value. Note that in the case of
:meth:`.PropagationSolution.q`, a :class:`.BeamParam` instance is still returned but
this will be a symbolic beam parameter (as indicated by the :attr:`.BeamParam.symbolic`
flag).

Symbolic beam propagation can result in very long symbolic expressions. Rudimentary
simplification routines are provided in Finesse but for large models keeping every
parameter symbolic is unnecessary. If you only need specific parameters to be kept
as symbols then you can specify which to keep, for example, `symbolic=('ITMX.Rcx',)`
or `symbolic=('ITMX.Rcx', 'ITMYlens.f')`. If simplification of symbolic beam propagation
is required you must specify `simplify=True`, it does not happen by default. This is
because simplifying when using every single parameter takes too long and is to
complicated to work with.

.. jupyter-execute:

    prop = model.propagate_beam(path=ITM_TO_SRM, symbolic=('lBS_SRM.L', ), simplify=True);
    prop.abcd()

.. rubric:: Plotting symbolic beam propagations

Plotting of the solution is supported, with the added option of substituting parameters
in the :meth:`.PropagationSolution.plot` call. An example is given below, where the RoC
of ZM1 is changed from :math:`R_c = -50` m to :math:`R_c = -70` m, and the distance
between the two telescope mirrors (ZM1 , ZM2) is reduced to 30 metres:

.. jupyter-execute::

    ps_sym.plot(subs={model.ZM1.Rcx: -70, model.lZM1_ZM2.L: 30});

If ``subs`` is not given, then the current value of each parameter will be used.

Animation of a symbolic :class:`.PropagationSolution` is also supported via the method
:meth:`.PropagationSolution.animate`. This method expects at least one array-like
parameter substitution in the ``subs`` dict in order to perform the animation over this
parameter scan. For example, we can see how the beam sizes and accumulated Gouy phases
change as we scan over the telescope length via:

.. jupyter-execute::
    :hide-output:

    import numpy as np

    # Vary telescope length from 20 m to 70 m
    tel_zs = np.linspace(20, 70, 60)
    fig, axes, anim = ps_sym.animate(
        {model.lZM1_ZM2.L: tel_zs}, interval=50,
    )

.. jupyter-execute::
    :hide-code:

    from IPython.display import HTML
    HTML(anim.to_jshtml())

.. rubric:: Evaluating symbolic beam properties

The true power of symbolic :class:`.PropagationSolution` objects comes from the ability
to obtain symbolic expressions for any geometric property of the beam at any point in
the traced path. Accessing these properties is performed in exactly the same way as
detailed in the "Accessing beam properties" section above.

As an example, we can obtain a symbolic expression for the beam size on ZM1 with:

.. jupyter-execute::

    # Symbolic expr of w at ZM1
    w_zm1_sym = ps_sym.w(model.ZM1.p2.i)

    # Evaluate this with current parameter values
    print(f"w_ZM1 = {w_zm1_sym.eval() / 1e-3} mm")

Similarly to the plotting routines, we can pass a ``subs`` dict to the ``eval`` method
of the symbolic expression to find the beam size at ZM1 using different dependent optic
parameters:

.. jupyter-execute::

    # Evaluate w_ZM1 with ZM2 RoC changed to -70 m
    print(f"w_ZM1 = {w_zm1_sym.eval(subs={model.ZM2.Rcx: -70}) / 1e-3} mm")

    # Evaluate w_ZM1 with ITM lens focal length changed to 100 m
    # and telescope length changed to 60 m
    print(f"w_ZM1 = {w_zm1_sym.eval(subs={model.ITM_lens.f: 100, model.lZM1_ZM2.L: 60}) / 1e-3} mm")

These symbolic expressions also support :class:`numpy.ndarray` type arguments, allowing
for evaluation over N-dimensional arrays. For example, we can evaluate the beam size at
ZM1 over a range of ITM lens focal lengths and plot the results:

.. jupyter-execute::

    import matplotlib.pyplot as plt

    fs = np.linspace(50, 90, 100)
    w_zm1s = w_zm1_sym.eval(subs={model.ITM_lens.f: fs}) / 1e-3 # scale to mm

    plt.plot(fs, w_zm1s)
    plt.xlabel("ITM lens focal length [m]")
    plt.ylabel("Beam size at ZM1 [mm]");

Or we could, for example, compute the Gouy phase accumulated up to ZM1 for the same
parameter space:

.. jupyter-execute::

    acc_gouy_zm1_sym = ps_sym.acc_gouy_up_to("ZM1")

    acc_gouys_zm1 = acc_gouy_zm1_sym.eval(subs={model.ITM_lens.f: fs})

    plt.plot(fs, acc_gouys_zm1)
    plt.xlabel("ITM lens focal length [m]")
    plt.ylabel("Accumulated Gouy phase from\nITM to ZM1 [deg]");

This is a very powerful feature, as it leverages the speed of NumPy array calculations
allowing for fast, high-dimensional grid calculations of arbitrary beam properties over
any dependent model parameter.


Reverse beam propagation
````````````````````````

Some components in a model may not have forwards and backwards couplings between
ports, such as a directional beamsplitter. These elements act as a form of optical
isolation, such as a faraday isolator would. This means you can only trace beams
through a physically viable path, as the beam tracer only follows the connections
at each component.

To get around this, it is possible to "reverse propagate" and let FINESSE know you
want to get around this limitation and follow unphysical paths. Here is a simple
example:

.. jupyter-execute::

    model = finesse.script.parse(
        """
    l l1
    dbs isolator
    lens L1 f=1
    m m1
    link(l1.p1, 1, isolator.p1, isolator.p3, 2, L1, 3, m1)
    gauss g1 l1.p1.o w0=1e-3 z=1e-2
    """
    )

    fwd = model.propagate_beam("l1.p1.o", "m1.p1.i")
    rev = model.propagate_beam("m1.p1.o", "l1.p1.i", reverse_propagate=True)

    fwd.plot()
    rev.plot()
