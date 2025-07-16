.. include:: /defs.hrst
.. _shiftedbeam_convergence:

Simulating off-axis beams: convergence and scaling
--------------------------------------------------

In order to properly describe off-axis beams in the paraxial approximation, a
number of higher order modes will need to be taken into account. How many modes
need to be included (by setting ``maxtem`` using the :kat:command:`modes`
command) will depend on the optical layout and geometry of the beam such as the
waist position :math:`z_0` and the relative distance of the different optical
components to each other. Typically, if we choose ``maxtem`` too small, we
expect to miss part of the beam in simulations, resulting in missing total beam
energy. In this section we will therefore study how this total beam energy - as
measured using a :kat:element:`power_detector_dc` detector - behaves when
changing the various parameters. In general, when increasing ``maxtem`` we
expect convergence towards the total energy put into the beam by the
:kat:element:`laser`.

We will use the same optical layout as used in the previous section
:ref:`shiftedbeam`, it is sketched again in :numref:`fig_shiftbeam_basic_conv`.

.. _fig_shiftbeam_basic_conv:
.. figure:: images/telescope.*
   :align: center
   :width: 80%

   Basic setup

Both beamsplitters are tilted over the same small angle :math:`\beta`, resulting
in a parallel off-axis beam. In order to correctly describe the physics at
increasingly larger :math:`\beta`, more and more higher order modes will need to
be taken into account. We will start by investigating this behaviour for
different values of ``maxtem``. We will then study the effects of changing the
two gaussian beam parameters waist position :math:`z_0` and waist size
:math:`w_0`. For :math:`z_0` we will find an optimal position within the optical
layout, providing the largest valid range of :math:`\beta`. Once we have found
the optimal :math:`z_0`, we will use it to study the valid :math:`\beta` range
as a function of the distance between the two beamsplitters, looking in
particular at the scaling under changes of the waist size :math:`w_0` and
linking it to the diffraction angle and Rayleigh range.

Base script for the different simulations
"""""""""""""""""""""""""""""""""""""""""

Each of the simulations will be using the following base script, consisting of a
:kat:element:`laser`, two :kat:element:`beamsplitter` elements and a
:kat:element:`power_detector_dc`. The different (default) parameters are defined
as python variables directly above the |Finesse| KatScript. The layout is the
same as that used in :ref:`example10`.

.. jupyter-execute::

    import matplotlib.pyplot as plt
    import finesse
    from finesse.analysis.actions import Xaxis
    finesse.configure(plotting=True)

    power = 1.5  # default laser power in Watt
    w0 = 10e-3   # default gaussian waist size in meter
    alpha = 30   # angle of incidence for both BS in degrees
    xbeta = 1e-5 # default tilt of both beam splitters in radians
    s1 = 1000    # distance till 1st beamsplitter
    s2 = 400     # distance between beamsplitters
    s3 = 600     # default distance between 2nd beamsplitter and detector
    z0 = -1200   # default waist position, 1200 meter right of laser
    maxtem = 7   # default maxtem

    basescript = f"""
    # Define laser and gaussian beam
    laser l1 P={power}
    gauss g1 l1.p1.o w0={w0} z={z0}
    modes(maxtem={maxtem})

    # define beamsplitters and their positions
    space s1 l1.p1 bs1.p1 L={s1}
    beamsplitter bs1 R=1 T=0 alpha={alpha} xbeta={xbeta}
    space s2 bs1.p2 bs2.p1 L={s2}
    beamsplitter bs2 R=1 T=0 alpha=bs1.alpha xbeta=bs1.xbeta
    space s3 bs2.p2 n1.p1 L={s3}

    # use photodiode to measure integrated power
    nothing n1
    power_detector_dc pd1 node=n1.p1.i pdtype=none
    """

    basekat = finesse.Model()
    basekat.parse(basescript)

Tilt angle dependency for various maxtem
""""""""""""""""""""""""""""""""""""""""

We run an appropriate action :kat:analysis:`xaxis` / :class:`.Xaxis`
for the model to measure the total energy as a function :math:`\beta` for 4
different values of ``maxtem``: 3, 5, 7 and 9 (set using the
:kat:command:`modes` / :meth:`.Model.modes` method)

.. jupyter-execute::

    # Set of maxtem values (and their plot styles)
    maxtem_vals   = [3, 5, 7, 9]
    maxtem_styles = ['r', 'b', 'g', 'm']

    # Define the model and set sweep parameter to xbeta
    kat1 = finesse.Model()
    kat1.parse(basescript)

    out1a = []
    for maxtem_val in maxtem_vals:
        kat1.modes(maxtem=maxtem_val)
        out1a.append(kat1.run(Xaxis(kat1.bs1.xbeta, "lin", 0, 4e-5, 40)))

In order to verify that the convergence is independent of the position of the
photodiode, i.e. to verify that the (diverging) gaussian beam emerging from the
second beamsplitter remains stable over larger distances, we repeat the
simulation at fixed ``maxtem`` = 7 at two different values of :math:`s_3`: 600
and 1000 meter:

.. jupyter-execute::

    # Set of s3 values (and their plot styles)
    s3_vals   = [600, 1000]
    s3_styles = ['g', 'b.']

    # Reset maxtem to its default
    kat1.modes(maxtem=maxtem)

    out1b = []
    for s3_val in s3_vals:
        kat1.spaces.s3.L = s3_val
        out1b.append(kat1.run(Xaxis(kat1.bs1.xbeta, "lin", 0, 4e-5, 40)))

We plot both results below:

.. jupyter-execute::

    f,ax = plt.subplots(ncols=2, figsize=(12, 5))

    ax[0].set_title(f"Effect maxtem on total power (distance {s3})")
    for sim in range(len(out1a)):
	ax[0].plot(out1a[sim].x[0], out1a[sim]['pd1'],
		   maxtem_styles[sim], label=f"maxtem {maxtem_vals[sim]}")

    ax[1].set_title(f"Effect distance on total power (maxtem {maxtem})")
    for sim in range(len(out1b)):
	ax[1].plot(out1b[sim].x[0], out1b[sim]['pd1'],
		   s3_styles[sim], label=f"{s3_vals[sim]} meter")

    for i in 0,1:
	ax[i].set_xlabel("xbeta (radian)")
	ax[i].set_ylabel("power (Watt)")
	ax[i].legend(loc="lower left")

From the left figure, it is clear that we can simulate larger and larger
:math:`\beta` if we include more and more higher order modes in the simulation,
i.e. by increasing the ``maxtem`` value. The valid range of :math:`\beta` scales
roughly linearly with ``maxtem``, meaning the number of modes to be includes
roughly scales quadratically (note that for ``maxtem`` = :math:`N` there are
:math:`\frac{1}{2}(N+1)(N+2)` modes).

.. todo:: do we expect convergence to break down at very high beta?

From the right figure we see that the convergence is not dependent on the
distance behind the second beamsplitter at which we measure the total power in
the bundle indicating a stable beam.

Waist position dependency
"""""""""""""""""""""""""

Above we showed that the total power in the beam can give a clear indication for
the range of validity of the tilt angle :math:`\beta` (for given ``maxtem``).
We will now use this dependency to find the optimal position of the waist of the
Gaussian beam :math:`z_0`.

We run simulations for two different ranges in :math:`z_0`: one for a larger
range that starts at the laser and ends at the detector, and a second for a
smaller range and closer to the two beamsplitters, where :math:`z_0` varies from
600 till 1600 meter to the right of the laser (i.e from 400 meter left of the
first beamsplitter till 200 meter right of the second beamsplitter).

Each simulation is run using an :kat:analysis:`x2axis` action varying both
:math:`\beta` and :math:`z_0`.
Note that as sweep parameter in the :kat:analysis:`x2axis` action we have to
specify either the x- or the y-coordinate of :math:`z_0`, so we need to make
sure that they change simultaneously. We can easily enforce this by making one a
reference to the other (see :ref:`expressions`).

.. jupyter-execute::

    # Define the model
    kat2 = finesse.Model()
    kat2.parse(basescript)

    # Make sure zy equals zx
    kat2.g1.zy = kat2.g1.zx.ref

    # vary beta and z0
    out2 = kat2.run(f"""
        series(
            x2axis(bs1.xbeta, lin, 1e-5, 4e-5, 40, g1.zx, lin, 0, -2000, 60, name="full"),
            x2axis(bs1.xbeta, lin, 1e-5, 4e-5, 40, g1.zx, lin, -600, -1600, 60, name="zoom")
        )
    """)

    # Plot the full and zoom results
    f,ax = plt.subplots(ncols=2, figsize=(11.8, 5))
    for (i, name) in ([0, 'full'], [1, 'zoom']):
	pxy_extent = (out2[name].x[0].min(), out2[name].x[0].max(),
		      out2[name].x[1].max(), out2[name].x[1].min())
	p = ax[i].imshow(out2[name]['pd1'].T, aspect='auto', extent = pxy_extent)
	plt.colorbar(p, ax = ax[i])
	ax[i].contour(out2[name]['pd1'].T, colors='white', extent = pxy_extent)
	ax[i].set_xlabel("xbeta (radian)")
	ax[i].set_ylabel("z0 (w.r.t. laser)")
	ax[i].set_title("total power (W)")
	ax[i].ticklabel_format(style='sci', scilimits=(0,0))

The left plot shows the full range from the laser till the detector, the right
plot the close-up range around the beamsplitters. We see that the optimal
location for the waist is around 1 km right of the laser, i.e. at the location
of the first beamsplitter.

Effect of waist size
""""""""""""""""""""

Now that we know the optimal location for the Gaussian waist, we will fix it
that position and then study the effect of changing the distance between
the two beamsplitters :math:`s_2` and the scaling behaviour with respect to the
waist size :math:`w_0` of the beam.

We do two simulations, one at the default waist size
:math:`w_0` = 10 mm and one at a slightly larger waist size
:math:`w_0` = 12 mm. We will vary :math:`s_2` between 0 and 1 km for
the first simulation and between 0 and about 1.4 km for the second (we will
discuss below why we choose this range). We do both simulations in a single
:kat:analysis:`series`, using the :kat:analysis:`change` action to adapt the
waist size. We again use a reference to ensure that the x- and y-components of
both :math:`w_0` and :math:`z_0` are changing simultaneusly (see
:ref:`expressions`).

.. jupyter-execute::

    # w0, beta-range and L for 10mm run
    w0a = 10
    beta_a = (1e-5, 4e-5)
    La = 1000

    # w0, beta-range and L for 12mm run
    w0b = 12
    beta_b = (1e-5/1.2, 4e-5/1.2)
    Lb = 1000*1.2**2

    # Define model
    kat3 = finesse.Model()
    kat3.parse(basescript)

    # could set w0x and w0y separately, easier to use ref for y
    kat3.g1.w0y = kat3.g1.w0x.ref

    # Position gauss at optimal z0: i.e. position bs1
    kat3.g1.zx = -1000
    kat3.g1.zy = kat3.g1.zx.ref

    # run both subsimulations using series()
    out3 = kat3.run(f"""
        series(
            change(g1.w0x={w0a}e-03),
            x2axis(bs1.xbeta, lin, {beta_a[0]}, {beta_a[1]}, 40, s2.L, lin, 0, {La}, 60, name="a"),
            change(g1.w0x={w0b}e-03),
            x2axis(bs1.xbeta, lin, {beta_b[0]}, {beta_b[1]}, 40, s2.L, lin, 0, {Lb}, 60, name="b")
        )
    """)

And we plot the results

.. jupyter-execute::

    # Plot results for both waist sizes
    f,ax = plt.subplots(ncols=2, figsize=(11.8, 5))
    for (i, name, w0val) in ([0, 'a', w0a], [1, 'b', w0b]):
	pxy_extent = (out3[name].x[0].min(), out3[name].x[0].max(),
		      out3[name].x[1].min(), out3[name].x[1].max())
	p = ax[i].imshow(out3[name]['pd1'].T, aspect='auto', extent = pxy_extent)
	plt.colorbar(p, ax = ax[i])
	ax[i].contour(out3[name]['pd1'].T, colors='white', extent = pxy_extent)
	ax[i].set_xlabel("xbeta (radian)")
	ax[i].set_ylabel("length s2 (meter)")
	ax[i].set_title(f"w0={w0val}mm, total power (W)")
	ax[i].ticklabel_format(axis='y', style='plain', scilimits=(0,0))
    # Also add rescaled contours from 'a' to 'b'
    pxy_extent_scaled = (out3['a'].x[0].min()/1.2,   out3['a'].x[0].max()/1.2,
			 out3['a'].x[1].min()*1.2**2,out3['a'].x[1].max()*1.2**2)
    ax[1].contour(out3['a']['pd1'].T, colors='green', linestyles='dotted',
		  extent = pxy_extent_scaled);

Looking first at the left picture, we see that the usable range of :math:`\beta`
roughly remains constant till a distance of about 200-300 meter and then quickly
decreases for larger distances. This distance of 200-300 meter corresponds
roughly to the Rayleigh range :math:`z_R`

.. math::
   z_R = \frac{\pi w_0^2}{\lambda}
   :label: eq_shiftbeam_rayleigh

which for the parameters in the left plot (:math:`w_0 = 10\textrm{mm}, \lambda =
1.064 \mu\textrm{m}`) is 295 meter. It would seem that our default value of 400
meter is not ideal, but we will discuss this further below when looking at the
maximum attainable shift :math:`\Delta`.

Although in general the largest usable :math:`\beta` increases as a function of
``maxtem`` we expect that it will be proportional to the diffraction angle
:math:`\Theta` given by (see Eq. (9.19) in :cite:`LivingReview`)

.. math::
   \Theta = \arctan \left(\frac{w_0}{z_R}\right) \approx \frac{w_0}{z_R}
   = \frac{\lambda}{\pi w_0}
   :label: eq_shiftbeam_Theta

which - again for the parameters used in the left plot - equals
:math:`3.4 \cdot 10^{-5}` radians (which happens to be very close to the maximum
:math:`\beta` at the smaller distances for the ``maxtem`` = 7 used in this
simulation).

We can now discuss the scaling behaviour as a function of the waist size
:math:`w_0`. From the diffraction angle :math:`\Theta` Eq.
:eq:`eq_shiftbeam_Theta` we expect the useable range of :math:`\beta` to scale
as :math:`1/w_0`, while from the Rayleigh range :math:`z_R` Eq.
:eq:`eq_shiftbeam_rayleigh` we expect the vertical range to scale as
:math:`w_0^2`.
To verify these assumptions, we use - for the right plot with the results of
:math:`w_0` = 12 mm - a 1.2 times smaller horizontal range and a 1.44 times
larger vertical range. We also plot in the same plot the contours (dotted
green) from the left plot, scaled with the same factors and see two identical
pictures, confirming that :math:`z_R` and :math:`\Theta` set the relevant
scales.

From Eq. :eq:`eq_shiftbeam_shift` in :ref:`shiftedbeam` it now follows that for
small angles :math:`\beta` the largest attainable shift :math:`\Delta` scales
(roughly) linearly with :math:`w_0`, i.e. the dimensionless shift in units of
:math:`w_0` is actually independent of that waist size :math:`w_0`. It is
therefore interesting to also look at the beam energy as a function of
:math:`\Delta/w_0` instead of as a function of :math:`\beta`. This will allow us
to find the optimal distance between the two beamsplitters, leading to the
largest range of validity for :math:`\Delta/w_0`. We could refactor the results
of the above simulation, but it is easier to do a simulation using
:math:`\Delta/w_0` directly in the :kat:analysis:`x2axis`.

.. jupyter-execute::

    import numpy as np

    # define model
    kat4 = finesse.Model()
    kat4.parse(basescript)

    # Position gaussian waist at the optimal location
    kat4.g1.zx = -1000
    kat4.g1.zy = kat3.g1.zx.ref

    # Introduce variable delta
    kat4.parse("var Delta 1.0")
    # xbeta follows from Delta/w0 = s2*sin(2*beta)
    kat4.bs1.xbeta = 0.5*np.arcsin(w0*kat4.Delta.ref/kat4.spaces.s2.L.ref)

    # Run actual simulation
    out4 = kat4.run(f"x2axis(Delta, lin, 0.0, 4.0, 40, s2.L, lin, 0, 1500, 60)")

and plot the result

.. jupyter-execute::

    f,ax = plt.subplots()
    pxy_extent = (out4.x[0].min(), out4.x[0].max(),
		  out4.x[1].min(), out4.x[1].max())
    p = ax.imshow(out4['pd1'].T, aspect='auto', extent = pxy_extent)
    plt.colorbar(p, ax = ax)
    ax.contour(out4['pd1'].T, colors='white', extent = pxy_extent)
    ax.set_xlabel("Delta/w0")
    ax.set_ylabel("length s2 (meter)")
    ax.set_title(f"w0={w0a}mm, total power (W)")
    ax.ticklabel_format(axis='y', style='plain', scilimits=(0,0))

We see that for larger distances between the two beamsplitters, the range of
validity of :math:`\Delta/w_0` becomes constant, it is then only dependent on
the number of higher order modes included.

To summarize, we have seen that the best location for the first beamsplitter is
at the waist of the Gaussian beam. Furthermore, the optimal position for the
second beamsplitter is just beyond the Rayleigh range. Going to larger distances
decreases the maximum tilt, going to smaller distances decreases the maximum
(dimensionless) beamshift. The maximum value of this beamshift is not dependent
on the waist size and for larger inter-beamsplitter distances depends only on
the number of higher order modes included in the simulation.
