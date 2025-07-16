.. include:: /defs.hrst
.. _beamsplitter_phase:


Mirror and beamsplitter phase relationships
*******************************************

As follows from Fresnel equations, passive optical components, such as mirrors,
beamsplitters and lenses, can be described as flat thin layers linearly coupling
with the incident light. When light impinges on that surface, both reflection
and refraction of the light may occur. A coupling coefficient is the ratio of
the reflected, or transmitted, light to the incident light. These ratios are
generally complex, describing not only the relative amplitudes but also the
phase shifts at the interface.

We define a mirror as the contact plane between two media. We derive the
coupling coefficients under normal incidence. To do that, we calculate the phase
jump for the transmitted and reflected beam. Finally, we obtain the phase for
the transmitted beam as a function pf the position of the mirror. We repeat
this for a beamsplitter, which is similar to a mirror except for an angle of
incidence.

The phase relationship used is important to conserve energy. A variety exist in
literature, the two most common are a symmetric and an anti-symmetric
relationship between the two sides of the optic. Many theory papers use the
anti-symmetric one which has no phase accumulated on transmission while one of the
sides accquires a 180 degree phase shift on reflection. The other - symmetric case - is where
neither reflection accumulates any phase, but the transmission from either side
has a 90 degree shift, or a factor of complex :math:`i`.

Both cases should be physically consistent, the differences will often be in the
detunings you need to use in a model to place an interferometer at the correct
operating point. In more complicated scenarios however, this phase relationship
must be expanded. Below also highlights the relationship between an arbitrary
reflecting surface with a different refractive index and angle of incidence
on each side. Since such a problem is no longer symmetric it must be treated more
carefully. By default |Finesse| will use the below relationships when it sees
such cases. This should be taken into account if you are doing analytic
comparisons.

.. note::

    You can revert back to |Finesse| v2 behaviour by changing the model settings.
    E.g. to use |Finesse| v2 phase convention on transmission, use
    `self._settings.phase_config.v2_transmission_phase = True`.

.. _mirror_coupling:

Mirror
######

The coupling of light field amplitudes with a mirror under normal incidence can
be described as follows: there are two coherent input fields, :math:`a_{1}`
impinging on the mirror on the front and :math:`a_{2}` on the back surface. Two
output fields leave the mirror, :math:`b_{1}` and :math:`b_{2}`.

.. _fig_flat_surface:
.. figure:: images/flat_surface.*
   :align: center

   Schematic for two coherent beams falling on the front and back surface of a
   mirror, :math:`a_{1}` and :math:`a_{2}`,  along with the outcoming beams
   reflected off either surface, :math:`b_{1}` and :math:`b_{2}`.

The following linear equations can be used to describe the coupling:

.. math::
    \begin{pmatrix}
        b_{1}\\
        b_{2}
    \end{pmatrix}
    =
    \begin{pmatrix}
        M_{11} & M_{12}\\
        M_{21} & M_{22}
    \end{pmatrix}
    \begin{pmatrix}
        a_{1} \\
        a_{2}
    \end{pmatrix}

Tuning
******

We define the *tuning* :math:`\phi_{0}` of a surface as the shift in the mirror
position (expressed in radians) with respect to the reference plane. A tuning
of :math:`\phi_{0}=2\pi` translates the mirror by one vacuum wavelength (default
wavelength set in |Finesse|): :math:`x=\lambda_{0}`. A positive tuning is
defined to be in the direction of the normal vector on the front surface.

When the mirror shifts its position :math:`x` meters, the corresponding
tuning becomes

.. math::
    \phi_{0} = k_{0} x = \frac{2\pi}{\lambda_{0}} x = \frac{\omega_{0}}{c} x

A certain displacement results in different changes in the optical path for
light fields with different frequencies. To take that into account, :math:`\phi`
can be generalised as follows:

.. math::
    \phi = \phi_{0} \frac{\omega}{\omega_{0}}

.. _fig_flat_surface_tuning:
.. figure:: images/flat_surface_tuning.*
   :align: center

   Tuning of a mirror under normal incidence. The solid line is the reference
   untuned mirror, the dashed one is tuned mirror.

**********

We define the light-mirror coupling coefficients as:

.. math::
    M_{11} = r e^{i\varphi_{11}(\phi)}\\
    M_{22} = r e^{i\varphi_{22}(\phi)}\\
    M_{12} = M_{21} = t e^{i\varphi_{12}(\phi)} = t e^{i\varphi_{21}(\phi)}

where :math:`r` is the amplitude reflectance of the mirror, :math:`t` the mirror
transmittance and :math:`\varphi_{ij}(\phi)` the phase jump at the surface for
the impinging light, which depends on the surface tuning. We assume
:math:`\varphi_{12}(\phi)` to be equal to :math:`\varphi_{21}(\phi)` not to
introduce a preferred direction of propagation.
For a loss-less surface we can compute conditions for :math:`\varphi_{ij}(\phi)`
from energy conservation:

.. math::
    |b_1|^{2} + |b_2|^{2} = a_1^{2} + a_2^{2}

while :math:`|b_1|^{2}` and :math:`|b_2|^{2}` can be expressed using the
definitions of the light-mirror coupling coefficients above as:

.. math::
    |b_{1}|^{2} = r^{2} a_{1}^{2} + t^{2} a_{2}^{2} + 2rt\, a_{1} a_{2}\,
    \cos(\varphi_{12}-\varphi_{11})\\

    |b_{2}|^{2} = t^{2} a_{1}^{2} + r^{2} a_{2}^{2} + 2rt\, a_{1} a_{2}\,
    \cos(\varphi_{12}-\varphi_{22})

Hence energy conservation requires:

.. math::
    \cos(\varphi_{12}-\varphi_{11}) = -\cos(\varphi_{12}-\varphi_{22})

which in turn requires:

.. math::
    \varphi_{12}-\varphi_{11} = (2N+1)\pi - (\varphi_{12}-\varphi_{22})

where :math:`N` is an integer. After some simple algebraic steps, we obtain:

.. math::
    \varphi_{12} = (2N+1)\frac{\pi}{2} + \frac{\varphi_{11}+\varphi_{22}}{2}

We arbitrarily set :math:`N=0` and we will follow this convention throughout
this modeling. In general, the conditions for :math:`\varphi_{ij}(\phi)` are
given by the following equation:

.. math::
    \varphi_{12}(\phi) = \frac{\pi}{2} + \frac{\varphi_{11}(\phi)+\varphi_{22}(\phi)}{2}


Phase jumps in the untuned case
*******************************

At the reference position of the mirror, we arbitrarily set
:math:`\varphi_{11}(0)=\varphi_{22}(0)=0` and we will follow this convention for
the rest of the modeling. Hence the phase gain for a beam transmitted through the
mirror at the reference position is:

.. math::
    \varphi_{12}(0) = \varphi_{21}(0) = \frac{\pi}{2}


Phase jumps in the tuned case
*****************************

As the mirror is tuned, the phase of a beam reflected off the surface
:math:`\varphi_{ii}(\phi)` is given as:

.. math::
    \varphi_{11}(\phi) = 2n_{1}\phi\\
    \varphi_{22}(\phi) = -2n_{2}\phi

where :math:`n_{1}` and :math:`n_{2}` are the indices of refraction of the media
on either side of the surface. Substituting :math:`\varphi_{11}(\phi)` and
:math:`\varphi_{22}(\phi)` in the  equation for :math:`\varphi_{ij}(\phi)`, we
obtain :math:`\varphi_{12}(\phi)`:

.. math::
    \varphi_{12}(\phi) = \varphi_{21}(\phi) = \frac{\pi}{2} + (n_{1}-n_{2})\phi


Beamsplitter
############

A beamsplitter is similar to a mirror except for the extra parameter
:math:`\alpha` which indicates the tilt angle relative to the incoming beams.

Reflection
**********

.. _fig_flat_surface_tilted_Refl:
.. figure:: images/flat_surface_tilted_Refl.*
   :align: center

   Schematic for the beam reflected off a beamsplitter in the reference (solid lines) and tuned case (dashed lines).

Referring to the figure above, we define the following geometrical paths as:

.. math::
    a = \frac{x}{\cos(\alpha)}\\[10pt]
    b = a\cos(2\alpha)\\[10pt]
    c = \frac{x}{\cos(\beta)}\\[10pt]
    d = c\cos(2\beta)\\[10pt]

where :math:`\beta` is the refraction angle given by Snell's law:

.. math::
    n_{1}\sin(\alpha) = n_{2}\sin(\beta)

The phase change for a beam reflected on one side of a beamsplitter is:

.. math::
    \left(\frac{\omega}{\omega_{0}}k_{0}\right)\, n_{1}|a + b| = 2 n_{1}\phi \cos{\alpha}\\[10pt]

As for a beam reflected on the other side, the phase change is:

.. math::
    \left(\frac{\omega}{\omega_{0}}k_{0}\right)\, n_{2}|c + d| = 2 n_{2}\phi \cos{\beta}\\[10pt]

**********

As was done for the mirror, we model the beamsplitter-light coupling via linear coefficients :math:`M_{ij}`:

.. math::
    M_{11} = r e^{i\varphi_{11}(\phi)}\\
    M_{22} = r e^{i\varphi_{22}(\phi)}\\
    M_{12} = M_{21} = t e^{i\varphi_{12}(\phi)} = t e^{i\varphi_{21}(\phi)}

The conditions for :math:`\varphi_{ij}(\phi)` are given by:

.. math::
    \varphi_{12}(\phi) = \frac{\pi}{2} + \frac{\varphi_{11}(\phi)+\varphi_{22}(\phi)}{2}

The phase change for a beam reflected off either side of the beamsplitter is given as:

.. math::
    \varphi_{11}(\phi) = 2n_{1}\phi\cos\alpha\\
    \varphi_{22}(\phi) = -2n_{2}\phi\cos\beta

Substituting :math:`\varphi_{11}(\phi)` and :math:`\varphi_{22}(\phi)` in the
equation for :math:`\varphi_{ij}(\phi)`, we obtain :math:`\varphi_{12}(\phi)`:

.. math::
    \varphi_{12}(\phi) = \varphi_{21}(\phi) = \frac{\pi}{2} + (n_{1}\cos\alpha-n_{2}\cos\beta)\phi
