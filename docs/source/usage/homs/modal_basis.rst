.. include:: /defs.hrst
.. _modal_basis:

Defining the modal basis
~~~~~~~~~~~~~~~~~~~~~~~~

Using the correct modal basis (via beam parameters) at each point in the model
configuration is important for computing couplings through :ref:`scatter_matrices` and
other quantities such as the spatial distribution :math:`u_{nm}` used for, e.g.,
calculating beam profiles. For these reasons, it is necessary to understand how
|Finesse| sets beam parameters at nodes and how you can use this information to model
the systems that you are investigating.

Complete details on the beam tracing algorithm of |Finesse| are included in
:ref:`tracing_manual`, this page will focus on how to use cavities and manual setting of
beam parameters.

Using cavity objects
````````````````````

It is important in any analysis involving cavities that HOM couplings are computed using
a (propagated) cavity eigenmode as the basis. If the beam parameters at nodes inside a
cavity get set to values which are not equal to the corresponding cavity eigenmode then
we can get mismatches where there shouldn't be any and / or matches where there should
be mismatches.

A key role of the :class:`.Cavity` object is to prevent this from happening by
influencing the beam tracing algorithm when instances of this class are present in a
model. The eigenmodes of all cavities in a model are computed upon adding it to the
model and these are then used to set beam parameters in the rest of the model via
propagating them through the ABCD matrix formalism, see :ref:`transforming_beam_param`
for details on this.

.. rubric:: Adding cavities to a model

:class:`.Cavity` objects can be added to a model in the standard way via a call to
:meth:`.Model.add` using the Python API. The constructor of the Cavity class requires a
`name` for the cavity and the `source` node (i.e. starting point) of the cavity; with an
optional `via` node if the cavity path should traverse via a specific node.

For simple cavities where the surfaces are not bounded (on the inside of the cavity
path) by AR / HR coating surfaces, you only ever need to specify the `source` node or
port. The cavity path will then be determined automatically from just this node. An
example of this is shown below, using a simple Fabry-Perot type cavity model:

.. jupyter-execute::

    import finesse
    import finesse.components as components

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-2.5
    s sCAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=2.5

    # Add a cavity object with source node as output node
    # of the second port of mirror ITM
    cav FP ITM.p2
    """)

    # showing the cavity path
    print(model.FP.draw())

If you have a model with AR / HR coatings at surfaces then you will need to specify the
`via` node argument in order for the cavity path to include the surface itself. See the
example below:

.. jupyter-execute::

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s l0 L0.p1 PRM.p1 L=1

    m PRM T=0.046 L=37.5u

    s lpr1 PRM.p2 BS.p1 L=9

    bs BS T=0.5 L=37.5u alpha=60

    s ly BS.p2 ITMAR.p1 L=100

    # Y arm input mirror with AR surface
    m ITMAR R=0 L=20u
    s ITMsub ITMAR.p2 ITM.p1 L=0.2 nr=1.44963098985906
    m ITM T=7000u L=37.5u Rc=-5580

    lambda(1550n)
    """)

    # INCORRECT -- doesn't include ITM
    model.add(components.Cavity("PRC_wrong", model.PRM.p2.o))

    print(model.PRC_wrong.draw())

You can see from this cavity path that the cavity "end" point is ITMAR, **not** ITM, as
by default the cavity will find the shortest path back to the source node. To correct
this, you must use the `via` node:

.. jupyter-execute::

    # CORRECT -- now we specify that the path must traverse via ITM
    model.add(components.Cavity("PRC", model.PRM.p2.o, model.ITM.p1.i))

    print(model.PRC.draw())

.. rubric:: Showing the importance of including Cavity objects

Below is an example highlighting the importance of including Cavity objects in a model
for HOM modelling. In this example, the detuning of the input mirror of a Fabry-Perot
cavity is scanned whilst the power on transmission is detected.

.. jupyter-execute::

    import finesse

    model = finesse.Model()
    # define a configuration with no Cavity object present
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-2.5
    s sCAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=2.5

    # set a beam parameter with w0 = 1 mm, z-z0 = -0.3 m at the laser
    gauss gL0 L0.p1.o w0=1m z=-0.3

    # detect the transmitted power
    pd trns ETM.p2.o

    modes(maxtem=4)
    """)

    out_no_cav = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")

Now we add a Cavity object to the model and run it again, we also display the mismatches
that are now present in the system due to the manually set beam parameter at the laser
and the beam parameters automatically set in the cavity based on its eigenmode. See
:ref:`mismatches_usage` for details on the :meth:`.Model.mismatches_table` method and
more.

.. jupyter-execute::

    # source node is output of second port of ITM
    # -> cavity path is then determined automatically from this
    model.add(finesse.components.Cavity("FP", model.ITM.p2.o))

    print(model.mismatches_table())

    out_cav = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")

And then plot the results of each to show the difference:

.. jupyter-execute::

    import matplotlib.pyplot as plt
    finesse.configure(plotting=True)

    fig = plt.figure()
    plt.semilogy(out_no_cav.x1, out_no_cav["trns"], label="No cavity tracing")
    plt.semilogy(out_cav.x1, out_cav["trns"], label="With cavity tracing")
    plt.xlabel("ETM phi [deg]")
    plt.ylabel("Transmitted power [W]")
    plt.legend();

From this plot you can see that with the cavity included we get the expected resonance
peaks corresponding to the fundamental (TEM00) mode and the 02/20 and 04/40 modes; the
latter arising due to the mismatch present within the system. Without the cavity object
present, the transmitted power trace obtained is nonsensical.


Defining manual beam parameters
```````````````````````````````

Before we demonstrated how we can use the cavity command to automatically set
the beam parameters in a model. In certain cases you may also need to specify
what beam parameter to use manually. For example, if you wanted the beam shape
coming out of a laser to be of a particular shape and see what mismatch there is
to a cavity.

Here we start with a cavity and defines the mode shape using the cavity command.

.. jupyter-execute::

    import finesse
    import finesse.components as components
    import matplotlib.pyplot as plt

    finesse.init_plotting()

    model = finesse.Model()
    model.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-2.5
    s sCAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=2.5
    cav FP ITM.p2
    # detect the transmitted power
    pd trns ETM.p2.o
    modes(maxtem=4)
    """)

    # run the simulation to get the ideal case
    out0 = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")

Next we use :kat:command:`gauss` command to set what shape the beam is at the
output of the laser, this will add a new object to the model called `g1`.

.. jupyter-execute::

    model.parse("gauss g1 L0.p1.o w0=1m z=0")
    out1 = model.run("xaxis(ETM.phi, lin, -180, 180, 500)")

We can see what beam parameters are being set by a "gauss" command using:

.. jupyter-execute::

    model.gausses

Lastly we can plot the result and see that we have scattering into higher order
modes due to the mismatch.

.. jupyter-execute::

    fig = plt.figure()
    plt.semilogy(out0.x1, out0["trns"], label="Ideal mode matching")
    plt.semilogy(out1.x1, out1["trns"], label="With manual laser shape")
    plt.xlabel("ETM phi [deg]")
    plt.ylabel("Transmitted power [W]")
    plt.legend();

We can also look at the table of mismatches that are ocurring in the ABCD beam
tracing calculations. With this we should see siginifcant mismatches around the
ITM mirror.

.. jupyter-execute::

    model.mismatches_table()

Remove a manual beam parameter or cavity
````````````````````````````````````````

To remove a manually set beam parameter or cavity command you can use the remove
command with the name of the target:

.. jupyter-execute::

    model.remove('g1')


Defining manual beam parameters with the Python API
```````````````````````````````````````````````````

Above we defined the beam parameter using the KatScript command. To achieve the
same result from using the Python API we need to look at the optical nodes `q`
property. We can print these easily:

.. jupyter-execute::

    model = finesse.script.parse("""
    l L0 P=1
    s s0 L0.p1 ITM.p1

    m ITM R=0.99 T=0.01 Rc=-2.5
    s sCAV ITM.p2 ETM.p1 L=1
    m ETM R=0.99 T=0.01 Rc=2.5
    cav FP ITM.p2
    # detect the transmitted power
    pd trns ETM.p2.o
    modes(maxtem=4)
    """)
    model.beam_trace()
    print(model.L0.p1.o.qx)
    print(model.L0.p1.o.qy)

We can also easily set them as well:

.. jupyter-execute::

    # qx and qy can be set separately or you can set .q to set both at the same time
    model.L0.p1.o.q = finesse.BeamParam(w0=1e-3, z=0)
    model.beam_trace()

Setting the beam parameter this way will not automatically perform a new beam
trace through the model. If you want to query the new beam parameters due to
your change you must run a :code:`model.beam_trace()`.

Doing the above will have the same effect as using a gaussian beam command as
seen before. If you check you will see a :class:`Gauss` is actually added to
the model with an auto-generated name and can be removed if needed.

.. jupyter-execute::

    model.gausses
