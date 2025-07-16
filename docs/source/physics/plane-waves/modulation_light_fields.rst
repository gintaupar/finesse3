.. include:: /defs.hrst
.. _modulation_of_light_fields:

Modulation of light fields
**************************

In principle, all parameters of a light field can be modulated. This section describes
the modulation of the amplitude, phase and frequency of the light. Any sinusoidal
modulation of amplitude or phase generates new field components that are shifted in
frequency with respect to the initial field. Basically, light power is shifted from one
frequency component, the carrier, to several others, the sidebands. The relative
amplitudes and phases of these sidebands differ for different types of modulation and
different modulation strenghts.

.. contents:: :local:

.. _phase_mod:

Phase Modulation
~~~~~~~~~~~~~~~~

Phase modulation can create a large number of sidebands. The number of sidebands with
noticeable power depends on the modulation strength (or depth) given by the *modulation
index* :math:`m`. Assuming an input field

.. math::
    E_{\mathrm{in}}=E_0~\exp{\left(\mathrm{i}\,\omega_0\,t\right)},


a sinusoidal phase modulation of the field can be described as

.. math::
    E=E_0~\exp\Bigl(\mathrm{i}\,(\omega_0\,t + m \cos{\left(\omega_m\,t\right)})\Bigr).


This equation can be expanded using the *Bessel functions of the first kind*
:math:`J_k(m)`. We can write

.. math::
    E = E_0 ~ \exp{\left(\mathrm{i}\,\omega_0\,t\right)} ~ \sum_{k=-\infty}^{\infty}
    \mathrm{i}^{\,k} ~ J_k(m) ~ \exp{\left( \mathrm{i} \, k \, \omega_m \, t \right)}.
    :label: bessel0


The field for :math:`k=0`, oscillating with the frequency of the input field
:math:`\omega_0`, represents the carrier. The sidebands can be divided into *upper*
(:math:`k>0`) and *lower* (:math:`k<0`) sidebands. These sidebands are light fields that
have been shifted in frequency by :math:`k\, \omega_m`. The upper and lower sidebands
with the same absolute value of :math:`k` are called a pair of sidebands of order
:math:`k`.

Equation :eq:`bessel0` shows that the carrier is surrounded by an infinite number of
sidebands.  However, for small modulation indices (:math:`m<1`) the Bessel functions
rapidly decrease with increasing :math:`k`, so we can use the approximation:

.. math::
    J_k(m)~=\frac{1}{k!}\left(\frac{m}{2}\right)^k+O\left(m^{k+2}\right).


In this case, only a few sidebands have to be taken into account. For
:math:`m\ll1` we can write

.. math::
    \begin{array}{lcl}
    E&=&E_0~\exp{\left(\mathrm{i}\,\omega_0\,t\right)}\\
    & & \times\Bigl(J_0(m)-\mathrm{i}\, J_{-1}(m)~\exp{\left(-\mathrm{i}\,
    \omega_m\,t\right)}+\mathrm{i}\, J_{1}(m)~\exp{\left(\mathrm{i}\,
    \omega_m\,t\right)}\Bigr),
    \end{array}


and with

.. math::
    J_{-k}(m)=(-1)^kJ_k(m),


we obtain

.. math::

    E=E_0~\exp{\left(\mathrm{i}\,\omega_0\,t\right)}~\left(1+\mathrm{i}\,
    \frac{m}{2}\Bigl(\exp{\left(-\mathrm{i}\,
    \omega_m\,t\right)}+\exp{\left(\mathrm{i}\, \omega_m\,t\right)}\Bigr)\right),


as the first-order approximation in :math:`m`. In the above equation the carrier field
remains unchanged by the modulation, therefore this approximation is not the most
intuitive. It is clearer if the approximation up to the second order in :math:`m` is
given:

.. math::

    E=E_0~\exp{\left(\mathrm{i}\,\omega_0\,t\right)}~\left(1-\frac{m^2}{4}+\mathrm{i}\,
    \frac{m}{2}\Bigl(\exp{\left(-\mathrm{i}\,
    \omega_m\,t\right)}+\exp{\left(\mathrm{i}\, \omega_m\,t\right)}\Bigr)\right),


which shows that power is transferred from the carrier to the sideband fields.

Higher-order expansions in :math:`m` can be performed simply by specifying the highest
order of Bessel function, which is to be used in the sum in Equation :eq:`bessel0`, i.e.

.. math::

    E=E_0~\exp{\left(\mathrm{i}\,\omega_0\,t\right)}
    ~\sum_{k=-order}^{order}i^{\,k}~J_k(m)
    ~\exp{\left(\mathrm{i}\, k \omega_m\,t\right)}.

.. _freq_mod:

Frequency modulation
~~~~~~~~~~~~~~~~~~~~

For small modulation, indices, phase modulation and frequency modulation can be
understood as different descriptions of the same effect :cite:`Heinzel99`. With the
frequency defined as :math:`f = d\varphi/dt` a sinusoidal frequency modulation can be
written as:

.. math::

    E=E_0~\mEx{\mathrm{i}\,\left(\omega_0\,t + \frac{\Delta\omega}{\omega_m}
    \cos{\left(\omega_m\,t\right)}\right)},

with :math:`\Delta\omega` as the frequency swing (how *far* the frequency is shifted by
the modulation) and :math:`\omega_m` the modulation frequency (how *fast* the frequency
is shifted). The modulation index is defined as:

.. math::
   m = \frac{\Delta\omega}{\omega_m}

.. _amp_mod:

Amplitude modulation
~~~~~~~~~~~~~~~~~~~~

In contrast to phase modulation, (sinusoidal) amplitude modulation always generates
exactly two sidebands.  Furthermore, a natural maximum modulation index exists: the
modulation index is defined to be one (:math:`m=1`) when the amplitude is modulated
between zero and the amplitude of the unmodulated field.

If the amplitude modulation is performed by an active element, for example by modulating
the current of a laser diode, the following equation can be used to describe the output
field:

.. math::
    \begin{array}{lcl}
    E&=&E_0~\mEx{\mathrm{i}\,\omega_0\,t}
    ~\Bigl(1+m\cos{\left(\omega_m\,t\right)}\Bigr)\\
    &=&E_0~\mEx{\mathrm{i}\,\omega_0\,t}~\Bigl(1+\frac{m}{2}~\mEx{\mathrm{i}\, \omega_m
   \,t}+\frac{m}{2}~\mEx{-\mathrm{i}\, \omega_m\,t}\Bigr).
    \end{array}


However, passive amplitude modulators (like acousto-optic modulators or electro-optic
modulators with polarisers) can only reduce the amplitude. In these cases, the following
equation is more useful:

.. math::
    \begin{array}{lcl}
    E&=&E_0~\mEx{\mathrm{i}\,\omega_0\,t}~\left(1-\frac{m}{2}\Bigl(1-\cos{\left(\omega_m
   \,t\right)}\Bigr)\right)\\
    &=&E_0~\mEx{\mathrm{i}\,\omega_0\,t}~\Bigl(1-\frac{m}{2}
    +\frac{m}{4}~\mEx{\mathrm{i}\, \omega_m
   \,t}+\frac{m}{4}~\mEx{-\mathrm{i}\, \omega_m\,t}\Bigr).
    \end{array}
