.. include:: /defs.hrst

.. _legacy_syntax:

=================
Finesse 2 support
=================

Supported legacy syntax
=======================

|Finesse| 3 supports most common |Finesse| 2 syntax via its legacy parsing mode. Legacy
parsing is available either via the :ref:`cli` (using :option:`--legacy <kat3
--legacy>`) or in Python using :meth:`.Model.parse_legacy` or
:meth:`.Model.parse_legacy_file`.

The |Finesse| 2 syntax supported by the legacy parser is listed below, along with links
to the relevent |Finesse| 3 KatScript syntax and Python API. Directives not listed below
are not supported. The supported parameters for each directive can be found in the
|Finesse| 2 `manual <http://www.gwoptics.org/finesse/download/manual.pdf>`__.

.. warning::

    Support for |Finesse| 2 script support is only provided as a convenience to assist
    in transitioning to |Finesse| 3. New scripts should always be written in |Finesse| 3
    script.

Components
----------

.. list-table::
    :header-rows: 1

    * - Legacy name(s)
      - New name(s)
      - Python API

    * - | ``bs``
        | ``bs1``
        | ``bs2``
        | ``beamsplitter``
      - | :kat:element:`bs`
      - | :class:`~.Beamsplitter`

    * - | ``cav``
        | ``cavity``
      - | :kat:element:`cav`
        | :kat:element:`cavity`
      - | :class:`~.Cavity`

    * - | ``dbs``
      - | :kat:element:`dbs`
      - | :class:`~.DirectionalBeamsplitter`

    * - | ``gauss``
        | ``gauss*``
        | ``gauss**``
      - | :kat:element:`gauss`
      - | :class:`~.Gauss`

    * - | ``gouy``
      - | :kat:element:`gouy`
      - | :class:`~.Gouy`

    * - | ``isol``
        | ``diode``
      - | :kat:element:`isol`
        | :kat:element:`isolator`
      - | :class:`~.Isolator`

    * - | ``l``
        | ``laser``
        | ``light``
      - | :kat:element:`l`
        | :kat:element:`laser`
      - | :class:`~.Laser`

    * - | ``lens``
        | ``lens*``
        | ``lens**``
        | ``lens***``
      - | :kat:element:`lens`
      - | :class:`~.Lens`

    * - | ``m``
        | ``m1``
        | ``m2``
        | ``mirror``
      - | :kat:element:`m`
        | :kat:element:`mirror`
      - | :class:`~.Mirror`

    * - | ``mod``
      - | :kat:element:`mod`
        | :kat:element:`modulator`
      - | :class:`~.Modulator`

    * - | ``sq``
        | ``squeezer``
      - | :kat:element:`sq`
        | :kat:element:`squeezer`
      - | :class:`~.Squeezer`

    * - | ``s``
        | ``space``
      - | :kat:element:`s`
        | :kat:element:`space`
      - | :class:`~.Space`

Detectors
---------

.. list-table::
    :header-rows: 1

    * - Legacy name(s)
      - New name(s)
      - Python API

    * - | ``ad``
      - | :kat:element:`ad`
      - | :class:`~.AmplitudeDetector`

    * - | ``beam``
      - | :kat:element:`ccd`
        | :kat:element:`ccdline`
        | :kat:element:`ccdpx`
      - | :class:`~.CCD`
        | :class:`~.CCDScanLine`
        | :class:`~.CCDPixel`

    * - | ``bp``
      - | :kat:element:`bp`
      - | :class:`~.BeamPropertyDetector`

    * - | ``xd``
      - | :kat:element:`xd`
      - | :class:`~.MotionDetector`

    * - | ``pd[n]``
      - | :kat:element:`pd`
        | :kat:element:`pd1`
        | :kat:element:`pd2`
      - | :class:`~.PowerDetector`
        | :class:`~.PowerDetectorDemod1`
        | :class:`~.PowerDetectorDemod2`

    * - | ``qnoised[S/N]``
      - | :kat:element:`qnoised`
        | :kat:element:`qnoised1`
        | :kat:element:`qnoised2`
      - | :class:`~.QuantumNoiseDetector`
        | :class:`~.QuantumNoiseDetectorDemod1`
        | :class:`~.QuantumNoiseDetectorDemod2`

    * - | ``qshot[S/N]``
      - | :kat:element:`qshot`
        | :kat:element:`qshot1`
        | :kat:element:`qshot2`
      - | :class:`~.QuantumShotNoiseDetector`
        | :class:`~.QuantumShotNoiseDetectorDemod1`
        | :class:`~.QuantumShotNoiseDetectorDemod2`

Commands
--------

.. list-table::
    :header-rows: 1

    * - Legacy name(s)
      - New name(s)
      - Python API

    * - | ``attr``
      - | Replaced with arguments to the relevent components.
      - | See relevent component constructors.

    * - | ``attr mass``
        | ``attr tf``
      - | Replaced with dedicated mechanics components e.g. :kat:element:`pendulum`, :kat:element:`free_mass`.
      - | See components under :mod:`components.mechanical<finesse.components.mechanical>`.

    * - | ``fadd``
        | (Implemented in |Finesse| 3 as
        | ``freq name frequency``)
      - | Not implemented.
      - | :meth:`Model.add_frequency<finesse.model.Model.add_frequency>`

    * - | ``fsig name f``
        | (Setting the frequency)
      - | :kat:command:`fsig`
      - | ``Model.fsig``

    * - | ``fsig name comp ...``
        | (Applying a signal)
      - | :kat:element:`sgen`
        | :kat:element:`signal_generator`
      - | :class:`~.SignalGenerator`

    * - | ``func``
        | ``set``
        | ``put``
        | ``put*``
      - | Replaced with KatScript :ref:`expressions`.
      - | :class:`~.Symbol`
        | :attr:`Parameter.ref<finesse.parameter.Parameter.ref>`

    * - | ``lambda``
        | ``lambda0``
      - | :kat:command:`lambda`
      - | :attr:`Model.lambda0<finesse.model.Model.lambda0>`

    * - | ``lock``
        | ``lock*``
      - | :kat:element:`lock`
      - | :class:`~.Lock`

    * - | ``mask``
      - | Replaced by ``mask`` property of the relevant detector.
      - | :attr:`MaskedDetector.mask<finesse.detectors.general.MaskedDetector.mask>`

    * - | ``maxtem``
      - | :kat:command:`modes`
      - | :meth:`Model.modes<finesse.model.Model.modes>`

    * - | ``pdtype``
      - | Use `pdtype` keyword argument for detectors
      - | Use `pdtype` keyword argument for detectors

    * - | ``phase``
      - | Not implemented.
      - | :meth:`Model.phase_level<finesse.model.Model.phase_level>`
        | ``Model.phase_config``

    * - | ``retrace``
      - | Not implemented.
      - | :attr:`Model.sim_trace_config<finesse.model.Model.sim_trace_config>`

    * - | ``startnode``
      - | Replaced by `priority` argument to :kat:element:`gauss` / :kat:element:`cav`.
      - | :attr:`TraceDependency.priority<finesse.components.trace_dependency.TraceDependency.priority>`

    * - | ``tem``
      - | :kat:command:`tem`
      - | :meth:`Laser.tem<finesse.components.laser.Laser.tem>`

    * - | ``var``
        | ``variable``
        | ``const``
      - | :kat:element:`var`
        | :kat:element:`variable`
      - | :class:`Variable<finesse.components.general.Variable>`

Axes
----

.. list-table::
    :header-rows: 1

    * - Legacy name(s)
      - New name(s)
      - Python API

    * - | ``noxaxis``
      - | :kat:analysis:`noxaxis`
        | `N.B.` This is now the default when no analyses are specified.
      - | :class:`~.Noxaxis`

    * - | ``xaxis``
        | ``x2axis``
        | ``x3axis``
      - | :kat:analysis:`xaxis`
        | :kat:analysis:`x2axis`
        | :kat:analysis:`x3axis`
      - | :class:`~.Xaxis`
        | :class:`~.X2axis`
        | :class:`~.X3axis`

    * - | ``yaxis``
      - | Not implemented.
      - | ``Model.yaxis``

    * - | ``scale``
      - | Not implemented.
      - | Not implemented.

.. _internal_differences:

Internal differences from Finesse 2
===================================

Modulator sideband amplitudes
-----------------------------

|Finesse| 2 uses a Bessel function routine optimised for speed of calculation. This
causes modulator field amplitudes to be incorrect by a small amount only for extreme
values of modulation depth when more than one modulation is performed. |Finesse| 3 uses
a correct implementation, and so will differ from |Finesse| 2 slightly.

The following table lists the modulation indices for which |Finesse| 2's error is
distinguishable from numerical noise (greater than :math:`10^{-14}`):

================  ===========================
Order, :math:`n`  Modulation index, :math:`m`
================  ===========================
1                 Always Correct
2                 :math:`0.482 \le m \le 2.0`
3                 :math:`0.958 \le m \le 3.0`
4                 :math:`2.128 \le m \le 4.0`
5                 :math:`2.920 \le m \le 5.0`
================  ===========================
