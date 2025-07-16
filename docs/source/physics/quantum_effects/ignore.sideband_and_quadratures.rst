.. include:: /defs.hrst
.. _sidebands_and_quadratures:

Sidebands and quadratures
~~~~~~~~~~~~~~~~~~~~~~~~~

The transfer functions that |Finesse| calculates (see :ref:`frequency_response_actions`)
are in terms of the upper and (conjugate) lower sidebands of optical fields. The
conjugate is returned as this is what is commonly used to solve linearised
systems in the sideband picture when radiation pressure (non-linear) effects are involved.

When analyzing quantum noise it is often more convenient to work with the
quadratures of the fields, i.e. the "two-photon" formalism. In the usual notation,
:math:`a_+(\omega)\equiv a(\omega_0 + \omega)` is the upper sideband around the carrier
frequency :math:`\omega_0` and :math:`a_-^\dag(\omega) \equiv a(\omega_0 - \omega)` is the
conjugate of the lower sideband. The cosine quadrature :math:`a_1(\omega)` and the sine
quadrature :math:`a_2(\omega)` are given by

.. math::
    \begin{bmatrix}
      a_1 \\
      a_2
    \end{bmatrix}
    =
    A
    \begin{bmatrix}
      a_+ \\
      a_-^\dag
    \end{bmatrix},
    \qquad
    A = \frac{1}{\sqrt{2}} \begin{bmatrix}
      1 & 1 \\
      -\mathrm{i} & \mathrm{i}
    \end{bmatrix},
    \qquad
    A^{-1} = \frac{1}{\sqrt{2}} \begin{bmatrix}
      1 & \mathrm{i} \\
      1 & -\mathrm{i}
    \end{bmatrix}.

.. note::

    There is nothing particularly quantum about this: we are just using the sine and
    cosine quadratures instead of the upper and lower sidebands. The only quantum part of
    "the quantum calculation" is when we use the same transfer functions to propagate a
    quantum state throughout the optomechanical system. We do need to be careful about
    normalizations, however. All of the calculations that |Finesse| does are in terms of
    electric fields (in units of :math:`\sqrt{\mathrm{W}}`).


Converting sidebands to quadratures: free mass mirror
-----------------------------------------------------

Whilst the above seems fairly straight forward we provide an example with a perfectly
reflecting free mass mirror. Consider a constant carrier field impinging on a mirror,
:math:`E_c`. We use the definition of power in an optical field as :math:`P = |E|^2`.
The result AC power fluctuation in an optical field is then :math:`P = 2 (E^{\ast}_c
E_{+} + E_c E^{\ast}_{-})`.

We then solve the problem of an incident upper and lower sideband present on this carrier
and how they beat with the carrier field to produce a power fluctuation. This power
fluctuation results in a radiation pressure force which moves the mirror. This mirror
motion then phase modulates the carrier field producing phase sidebands on reflection.
We need to solve for the output optical fields, as well as the motion and force on the
optic (Note that you can roll the force and motion into one equation but we separate
here for clarity).

.. jupyter-execute::

  import sympy as sy

  k, P, m, omega, lambda0, c = sy.var("k P M Omega lambda0 c", real=True, positive=True)

  E1u, E1l, E2u, E2l, z, F = sy.var(
      "E_{1+} E^{\\dagger}_{1-} E_{2+} E^{\\dagger}_{2-} z F"
  )

  Ec = sy.sqrt(P)  # Carrier just sqrt power

  # Solve the output optical fields, motion, and force
  eqs = sy.FiniteSet(
      E2u - (E1u + sy.I * 2 * k * z * Ec),
      E2l - (E1l - sy.I * 2 * k * z * Ec),
      z - (F * 1 / (-m * omega**2)),
      # 2 for momentum change of the light on reflection
      # Another 2 because
      F - (-2 * 2 * (sy.conjugate(Ec) * E2u + Ec * E2l) / c),
  )

  analytic = sy.solve(eqs, (E2u, E2l, z, F))


.. note::

  This is more complicated if the mirror is not perfectly reflective and more
  care should be taken to compute the momentum of all four fields incoming and
  outgoing from the mirror in practice, as |Finesse| does.

We find the reflected sidebands and mirror motion to be

.. jupyter-execute::

  display(sy.Eq(E2l, analytic[E2l].collect(-4 * sy.I * P * k)));
  display(sy.Eq(E2u, analytic[E2u].collect(-4 * sy.I * P * k)));
  display(sy.Eq(z, analytic[z]));

Converting to quadratures in then a case of applying the above transformation matrix

.. jupyter-execute::

  A = 1 / sy.sqrt(2) * sy.Matrix([[1, 1], [-sy.I, sy.I]])
  Q = A @ sy.Matrix([[analytic[E2u]], [analytic[E2l]]])
  Q.simplify();

  Q

If we substitute in 1W of total sidebands (1/2W in each of the upper and lower)
we find the standard result (See Eq 2.29 :cite:`Miao_2010`) for the reflected
quadrature from a free mass mirror:

.. jupyter-execute::

  Q1 = Q.subs({E1l: 1 / sy.sqrt(2), E1u: 1 / sy.sqrt(2)})

  display(sy.Eq(sy.var('a_{1}'), Q1[0]));
  display(sy.Eq(sy.var('a_{2}'), Q1[1]));



Comparison
----------

We can do a comparison with these analytics and a |Finesse| model using the
:class:`.FrequencyResponse2` and :class:`.FrequencyResponse3` actions. These
can be used to compute how optical fields couple into signal
(electrical/mechanical) states and vice-versa. The model below is a simple
free-mass mirror with some carrier field from a laser incident on it.
We plot the transfer functions for the sidebands and motions and compare them
to the analytics above.

.. jupyter-execute::

  import numpy as np
  import finesse

  import matplotlib.pyplot as plt
  from finesse.analysis.actions import FrequencyResponse3, FrequencyResponse2, Series
  import finesse.components as fc
  import scipy.constants as sc

  finesse.init_plotting()

  model = finesse.Model()
  model.add(fc.Laser("l1", P=1))
  model.add(fc.Mirror("m1", R=1, T=0))
  model.connect(model.l1.p1, model.m1.p1, L=0)
  model.add(fc.FreeMass("m1_sus", model.m1.mech, mass=1))
  model.fsig.f = 1

  f = np.geomspace(0.01, 1, 21)
  fsig = model.fsig.f.ref
  sol = model.run(
      Series(
          FrequencyResponse3(
              f,
              [("m1.p1.i", fsig), ("m1.p1.i", -fsig)],
              [("m1.p1.o", fsig), ("m1.p1.o", -fsig)],
              name="3",
          ),
          FrequencyResponse2(
              f,
              [("m1.p1.i", fsig), ("m1.p1.i", -fsig)],
              ["m1.mech.z", "m1.mech.F_z"],
              name="2",
          ),
      )
  )

  # propagate unit amount of upper and lower sidebands
  # for the simulated results
  ROOT2 = np.sqrt(2)  # need to inject 0.5W in each sideband
  Eout = sol["3"].out.squeeze() @ [1 / ROOT2, 1 / ROOT2]
  zout = sol["2"].out.squeeze() @ [1 / ROOT2, 1 / ROOT2]

  # Substitution values for this simulation state
  values = {
      m: model.m1_sus.mass.value,
      P: model.l1.P.value,
      k: model.k0,
      E1u: 1 / ROOT2,
      E1l: 1 / ROOT2,
      c: sc.c,
  }

  analytic_upper = sy.lambdify(
      omega,
      analytic[E2u].subs(values),
  )(f * 2 * np.pi)

  analytic_z = sy.lambdify(
      omega,
      analytic[z].subs(values),
  )(f * 2 * np.pi)

  analytic_F = sy.lambdify(
      omega,
      analytic[F].subs(values),
  )(f * 2 * np.pi)

  # assert np.allclose(Eout[:, 0], analytic_upper, atol=1e-14)
  # assert np.allclose(zout[:, 0], analytic_z, atol=1e-14)
  # assert np.allclose(zout[:, 1], analytic_F, atol=1e-14)

  plt.figure()
  plt.loglog(f, abs(Eout[:, 0]), label="upper")
  plt.loglog(f, abs(analytic_upper), label="Analytic upper", ls="--", c="r")
  plt.legend()
  plt.xlabel("Frequency [Hz]")
  plt.ylabel("Amplitude [$\sqrt{\mathrm{W}}$]")
  plt.title("Optical field outputs")

  plt.figure()
  plt.loglog(f, abs(zout[:, 0]), label="$z$")
  plt.loglog(f, abs(analytic_z), label="Analytic $z$", ls="--", c="r")
  plt.legend()
  plt.xlabel("Frequency [Hz]")
  plt.ylabel("Amplitude [m]")
  plt.title("Optical fields to motion")

  plt.figure()
  plt.semilogx(f, abs(zout[:, 1]) / 1e-9, label="$F_z$")
  plt.semilogx(
      f, abs(analytic_F) * np.ones_like(f) / 1e-9, label="Analytic $F_z$",
      ls="--", c="r",
  )
  plt.legend()
  plt.xlabel("Frequency [Hz]")
  plt.ylabel("Amplitude [nano-Newton]")
  plt.title("Optical fields to force")
  plt.margins(y=0.01);


Finally we can compute the quadratures from these sidebands and compare to
above.

.. jupyter-execute::

  nA = np.array(A).astype(complex)
  nQ = nA @ Eout.T

  analytic_a1 = sy.lambdify(
      omega,
      Q1[0].subs(values),
  )(f * 2 * np.pi)

  analytic_a2 = sy.lambdify(
      omega,
      Q1[1].subs(values),
  )(f * 2 * np.pi)

  plt.loglog(f, nQ[0], label='Finesse $a_1$')
  plt.loglog(f, nQ[1], label='Finesse $a_2$')
  plt.loglog(f, analytic_a1 * np.ones_like(f), label='Analytic $a_1$', ls=':', lw=3, c='g')
  plt.loglog(f, analytic_a2 * np.ones_like(f), label='Analytic $a_2$', ls=':', lw=3, c='k')
  plt.legend()
  plt.xlabel("Frequency [Hz]")
  plt.ylabel("Quadrature [$\sqrt{\mathrm{W}}$]")
  plt.title("Sideband to quadrature vs analytic");
