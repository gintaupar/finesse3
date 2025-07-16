.. include:: /defs.hrst

.. _ignore_two_photon_cavity:

Two-photon transfer functions and a Fabry-Perot cavity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we analyze a simple Fabry-Perot cavity where the ETM is a
free mass of mass :math:`m` and the ITM is fixed.

.. jupyter-execute::

    import numpy as np
    import scipy.constants as scc
    import finesse
    import finesse.components as fc
    import finesse.analysis.actions as fa
    import finesse.detectors as fd
    from finesse.plotting import bode

    finesse.init_plotting()

    L_arm = 4000
    model = finesse.Model()
    LASER = fc.Laser("Laser", P=5e3)
    ETM = fc.Mirror("ETM", T=0, L=0, Rc=2245)
    ITM = fc.Mirror("ITM", T=0.01, L=0, Rc=1940)
    ITMAR = fc.Mirror("ITMAR", T=1, L=0, Rc=np.inf)

    # Link and add all the components to the model
    model.link(LASER, ITMAR.p2, ITMAR.p1, ITM.p2, ITM.p1, L_arm, ETM.p1)

    cavARM = model.add(fc.Cavity("cavARM", model.ETM.p1))
    ETM_sus = model.add(fc.FreeMass("ETM_sus", model.ETM.mech, mass=40))
    model.add(fd.PowerDetector("Parm", model.ETM.p1.o))


Since |Finesse| does not (yet) have real mirrors, you have to manually add an HR surface
and a perfectly transmisive AR surface. Most realistic models do this. It is important
in this case because the convention in |Finesse| is that transmission through each
surface is imaginary :math:`\mathrm{i}t`. The convention for real thick mirrors is thus
that transmission is real and negative :math:`-t`. Having the transmission be imaginary
would lead to non-intuitive and non-standard two-photon transfer matrices.

.. warning::

    The AR and HR surfaces should be appropriately linked. We are not moving or rotating
    the ITM or having radiation pressure act on it, so we can be lazy and not do it here,
    but to be safe the AR parameters `phi`, `xbeta`, and `ybeta` should be symbolically
    linked to the corresponding HR parameters and all of the mechanical degrees of
    freedom of the AR surface should be connected to the corresponding degrees of freedom
    of the HR surface. Thick components will eventually be added which will make this step
    unnecessary.

.. todo::
    Add a section explaining how to do this being careful about the gains depending on
    whether the surfaces are pointing towards each other or not.

Calculation of the transfer functions
*************************************

As described in :ref:`radiation_pressure_mirror`, the light exerts a radiation pressure
force on the free mirror. In the quadrature picture, this induces a coupling from the
amplitude to the phase quadrature with strength

.. math::
    \mathcal{K} = \frac{4k(2R + L)\chi_0P}{c} \rightarrow - \frac{8kP}{m\omega^2c},

where :math:`P` is the power in the cavity, :math:`R` and :math:`L` are the reflectivity
and loss of the mirror, respectively, :math:`\chi_0` is the free susceptibility of the
mirror, and :math:`k=2\pi/\lambda_0` is the wavenumber of the laser. The mirror in this
example is a free perfectly reflecting mirror with no loss.

If :math:`a_1` and :math:`a_2` are the quadratures of the fields incident on the ITM,
:math:`b_1` and :math:`b_2` are the quadratures of the fields reflected from the ITM,
:math:`x` is the position of the ETM, and :math:`F` is an external force exerted on the
ETM, the transfer functions are

.. math::
    \begin{aligned}
    \begin{bmatrix}
      b_1 \\
      b_2
    \end{bmatrix}
    &=
    \underbrace{\begin{bmatrix}
      r_a & 0 \\
      -t_a^2 \mathcal{K} & r_a
    \end{bmatrix}}_{\texttt{FrequencyResponse3}}
    \begin{bmatrix}
      a_1 \\
      a_2
    \end{bmatrix}
    +
    \underbrace{2kt_a\sqrt{P}
    \begin{bmatrix}
      0 \\
      1
    \end{bmatrix}}_{\texttt{FrequencyResponse4}}
    x\\
    x &= \underbrace{-\frac{4t_a\sqrt{P}}{m\omega^2 c}
    \begin{bmatrix}
      1 & 0
    \end{bmatrix}}_{\texttt{FrequencyResponse2}}
    \begin{bmatrix}
      a_1 \\
      a_2
    \end{bmatrix}
    \underbrace{- \frac{1}{m\omega^2}}_{\texttt{FrequencyResponse}} F
    \end{aligned}

where :math:`r_a` and :math:`t_a` are the reflection from the arm and the transmission
through the arm, respectively. When the cavity is strongly overcoupled, as it is here,
and for frequencies far below the FSR

.. math::
    r_a = \frac{1 - \mathrm{i}\omega/\gamma}{1 + \mathrm{i}\omega/\gamma}, \qquad
    t_a = -\sqrt{\frac{2\mathcal{F}}{\pi}} \frac{1}{1 + \mathrm{i}\omega/\gamma}

where :math:`\gamma` is the cavity pole and :math:`\mathcal{F}` is the cavity finesse.

To compute these four transfer functions in |Finesse| we need all four of the
``FrequencyResponse`` actions. We will calculate the transfer functions for both a
cavity on resonance and for a detuned cavity. We therefore define a function
``analysis`` that returns all necessary actions to simplify the calculation.

.. jupyter-execute::

    F_Hz = np.geomspace(1, 1e3, 200)
    model.fsig.f = 1
    fsig = model.fsig.f.ref

    def analysis(key):
        """Actions needed for computing the two-photon transfer functions"""
        return [
            # for DC arm power
            fa.Noxaxis(name=f"DC_{key}"),
            # for field reflection off of and transmission through the cavity
            fa.FrequencyResponse3(
                F_Hz,
                [
                    (ITMAR.p2.i, +fsig),
                    (ITMAR.p2.i, -fsig),
                    (ETM.p1.o, +fsig),
                    (ETM.p1.o, -fsig),
                ],
                [
                    (ITMAR.p2.o, +fsig),
                    (ITMAR.p2.o, -fsig),
                ],
                name=f"fresp3_{key}",
            ),
            # for ETM motion to fields
            fa.FrequencyResponse4(
                F_Hz,
                [ETM.mech.z],
                [
                    (ITMAR.p2.o, +fsig),
                    (ITMAR.p2.o, -fsig),
                ],
                name=f"fresp4_{key}",
            ),
            # for mechanical modification
            fa.FrequencyResponse(
                F_Hz,
                [ETM.mech.z, ETM.mech.F_z],
                [ETM.mech.z],
                name=f"fresp_{key}",
            ),
            # for radiation pressure on mirror
            fa.FrequencyResponse2(
                F_Hz,
                [
                    (ITMAR.p2.i, +fsig),
                    (ITMAR.p2.i, -fsig),
                ],
                [ETM.mech.z],
                name=f"fresp2_{key}",
            ),
        ]

    sol = model.run(
        fa.Series(
            fa.PrintModelAttr("ETM.phi"),
            fa.Change({ETM.phi: 0}),
            # transfer functions for the cavity on resonance
            *analysis("tuned"),
            # detune the cavity by 1 degree
            fa.Change({ETM.phi: -1}, relative=True),
            # transfer functions for the detuned cavity
            *analysis("detuned"),
            fa.PrintModelAttr("ETM.phi"),
        )
    )

Note that the state of the model has changed at the end of the calculation so that the
cavity is now detuned. To prevent that, we could have made the detuned analysis
temporary using either :func:`.temporary`, :class:`.Temporary`, or
:class:`.TemporaryParameters`, but did not in this case to simplify accessing the
results.

The results of :class:`.FrequencyResponse3` are a matrix of sideband transfer functions
:math:`M_\mathrm{sb}`, those of :class:`.FrequencyResponse4` a vector
:math:`v_\mathrm{sb}` of sideband transfer functions, those of
:class:`.FrequencyResponse2` an adjoint vector of sideband transfer functions, and that
of :class:`.FrequencyResponse` a scalar not involving optical fields. (Actually two
scalars in this case because we asked for the response of mirror motion to both mirror
motion and an external force.) The conversions to two-photon transfer functions are

.. math::
   M_{2\mathrm{p}} = A M_\mathrm{sb} A^{-1}, \qquad
   v_{2\mathrm{p}} = A v_\mathrm{sb}, \qquad
   v_{2\mathrm{p}}^\dag = v_\mathrm{sb}^\dag A^{-1}

.. jupyter-execute::

    A2 = np.array([
        [1,  1],
        [-1j, 1j],
    ]) / np.sqrt(2)
    A2i = np.array([
        [1, 1j],
        [1, -1j],
    ]) / np.sqrt(2)

    def extract_2p_tfs(key):
        from types import SimpleNamespace
        tfs = dict(
            refl = A2 @ sol[f"fresp3_{key}"].out[..., :2, 0, 0] @ A2i,
            trans = A2 @ sol[f"fresp3_{key}"].out[..., 2:, 0, 0] @ A2i,
            motion = A2 @ sol[f"fresp4_{key}"].out[..., 0],
            mech = sol[f"fresp_{key}"].out,
            rp = sol[f"fresp2_{key}"].out[..., 0] @ A2i,
        )
        return SimpleNamespace(**tfs)

    tuned = extract_2p_tfs("tuned")
    detuned = extract_2p_tfs("detuned")

Cavity on resonance
*******************

First consider the behavior of a cavity on resonance described by the above equations.
We did not need to add a :class:`.Cavity` to this model since we are not modeling HOMs
and could easily calculate the cavity pole and finesse by hand, but added the
:class:`.Cavity` anyway to illustrate how you can use it to quickly get these parameters
(and others) from a more complicated model.

.. jupyter-execute::

    k = model.k0
    M = ETM_sus.mass.value
    Parm = sol["DC_tuned"]["Parm"]
    Krp = -8 * k * Parm / (M * (2 * np.pi * F_Hz)**2 * scc.c)
    Fp_Hz = cavARM.pole
    Fa = cavARM.finesse
    print(f"Arm power: {Parm * 1e-6:0.1f} MW")
    print(f"Arm cavity pole: {Fp_Hz:0.1f} Hz")
    print(f"Arm finesse: {Fa:0.0f}")
    r_arm = (1 - 1j * F_Hz / Fp_Hz) / (1 + 1j * F_Hz / Fp_Hz)
    t_arm = -np.sqrt(2 * Fa / np.pi) * 1 / (1 + 1j * F_Hz / Fp_Hz)

First look at the reflection of fields from the cavity as calculated by
:class:`.FrequencyResponse3`. The reflection of the phase to phase and amplitude to
amplitude quadratures are the same. Amplitude fluctuations are converted to phase
fluctuations through radiation pressure, but there is no conversion of phase
fluctuations to amplitude fluctuations.

.. jupyter-execute::

    axs = bode(F_Hz, tuned.refl[..., 0, 0], db=False, label="simulation")
    bode(F_Hz, r_arm, axs=axs, db=False, ls="--", label="theory")
    axs[0].set_ylim(0.1, 10)
    axs[0].set_ylabel("Magnitude [$\sqrt{\mathrm{W}}$ / $\sqrt{\mathrm{W}}$]")
    axs[0].set_title("Cavity amplitude to amplitude reflection")

    axs = bode(F_Hz, tuned.refl[..., 1, 0], db=False, label="simulation")
    bode(F_Hz, -t_arm**2 * Krp, axs=axs, db=False, ls="--", label="theory")
    axs[0].set_ylabel("Magnitude [$\sqrt{\mathrm{W}}$ / $\sqrt{\mathrm{W}}$]")
    axs[0].set_title("Cavity amplitude to phase reflection")
    print(
        "Cavity amplitude reflection equal to phase reflection?",
        np.allclose(tuned.refl[..., 1, 1], tuned.refl[..., 0, 0]),
    )
    print(
        "Cavity phase to amplitude reflection zero?",
        np.allclose(tuned.refl[..., 0, 1], 0),
    )

Next look at the response of the optical fields reflected from the cavity to mirror
motion as calculated by :class:`.FrequencyResponse4`.

.. jupyter-execute::

    axs = bode(F_Hz, tuned.motion[..., 1, 0], db=False, label="simulation")
    bode(F_Hz, -2 * k * t_arm * np.sqrt(Parm), axs=axs, db=False, ls="--", label="theory")
    axs[0].set_ylabel("Magnitude [$\sqrt{\mathrm{W}}$ / m]")
    axs[0].set_title("Mirror motion to reflected phase");

There is no modification to the mirror dynamics and so its susceptibility is still that
of a free mass as shown by the results of :class:`.FrequencyResponse`:

.. jupyter-execute::

    print(
        "Unmodified?",
        np.allclose(tuned.mech[..., 0, 0], np.ones_like(tuned.mech[..., 0, 0])),
    )
    print(
        "Free mass?",
        np.allclose(tuned.mech[..., 0, 1], -1 / (M * (2 * np.pi * F_Hz)**2)),
    )

Finally, the radiation pressure from the amplitude quadrature of the incident field is
calculated by :class:`.FrequencyResponse4`. The phase quadrature does not source any
radiation pressure in this case.

.. jupyter-execute::

    axs = bode(F_Hz, tuned.rp[..., 0, 0], db=False, label="simulation")
    bode(
        F_Hz, -4 * t_arm * np.sqrt(Parm) / (M * (2 * np.pi * F_Hz)**2 * scc.c),
        axs=axs, db=False, ls="--", label="theory",
    )
    axs[0].set_ylabel("Magnitude [m / $\sqrt{\mathrm{W}}$]")
    axs[0].set_title("Amplitude to mirror motion")
    print("Phase to mirror motion zero?", np.allclose(tuned.rp[..., 0, 1], 0))

Detuned cavity
**************

Detuning the cavity creates an optical spring. Notably, the mechanics of the mirror are
modified and the susceptibility is no longer that of a free mass

.. jupyter-execute::

    axs = bode(F_Hz, detuned.mech[..., 0, 0], db=False)
    axs[0].set_title("Mechanical modification")
    axs = bode(F_Hz, detuned.mech[..., 0, 1], db=False)
    axs[0].set_ylabel("Magnitude [m / N]")
    axs[0].set_title("Radiation pressured modified susceptibility");

Both quadratures now cause mirror motion

.. jupyter-execute::

    axs = bode(F_Hz, detuned.rp[..., 0, 0], db=False, label="amplitude")
    bode(F_Hz, detuned.rp[..., 0, 1], axs=axs, db=False, label="phase")
    axs[0].set_ylabel("Magnitude [m / $\sqrt{\mathrm{W}}$]")
    axs[0].set_title("Field to mirror motion");

and radiation pressure causes both quadratures to mix on reflection of the cavity

.. jupyter-execute::

    axs = bode(
        F_Hz, detuned.refl[..., 1, 0], db=False, label="amplitude to phase",
    )
    bode(
        F_Hz, detuned.refl[..., 0, 1], axs=axs, db=False, ls="--",
        label="phase to amplitude",
    )
    axs[0].set_ylabel("Magnitude [$\sqrt{\mathrm{W}}$ / $\sqrt{\mathrm{W}}$]")
    axs[0].set_title("Field reflection off of the cavity");
