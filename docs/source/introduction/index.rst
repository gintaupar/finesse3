.. include:: /defs.hrst

.. _intro:

Introduction to Finesse
=======================

|Finesse| is a simulation program for interferometers.  The user can build any kind of
virtual laser interferometer using the following components:

- lasers, with user-defined power, wavelength and shape of the output beam;
- free spaces with arbitrary index of refraction;
- mirrors and beam splitters, with flat or spherical surfaces;
- modulators to change amplitude and phase of the laser light;
- amplitude or power detectors with the possibility of demodulating the detected signal
  with one or more given demodulation frequencies;
- lenses and isolators.

For a given optical setup, the program computes the light field amplitudes at every
point in the interferometer assuming a steady state.  To do so, the interferometer
description is translated into a set of linear equations that are solved numerically.
For convenience, a number of standard analyses can be performed automatically by the
program, namely computing modulation-demodulation error signals and transfer functions.
|Finesse| can perform the analysis using plane waves or Hermite-Gauss modes. The latter
allows computation of the effects of mode matching and misalignments. In addition, error
signals for automatic alignment systems can be simulated.

.. figure:: /images/pound_drever_hall01.*
    :align: center

    A schematic diagram of a laser interferometer which can be modelled using |Finesse|
    (in this case a Fabry-Perot cavity with a Pound-Drever-Hall control scheme).


Literally every parameter of the interferometer description can be tuned during the
simulation. The typical output is a plot of a photodetector signal as a function of one
or two parameters of the interferometer (e.g. arm length, mirror reflectivity,
modulation frequency, mirror alignment).  Optional text output provides information
about the optical setup including, but not limited to, mode mismatch coefficients,
eigenmodes of cavities and beam sizes.

|Finesse| provides a fast and versatile tool that has proven to be very useful during
design and commissioning of interferometric gravitational wave detectors.  However, the
program has been designed to allow the analysis of arbitrary, user-defined optical
setups. In addition, it is easy to install and easy to use.  Therefore |Finesse| is very
well suited to study basic optical properties, like, for example, the power enhancement
in a resonating cavity or modulation-demodulation methods.

Motivation and History
----------------------

The search for gravitational waves with interferometric detectors has led to a new type
of laser interferometer: new topologies are formed combining known interferometer types.
In addition, the search for gravitational waves requires optical systems with a very
long baseline, large circulating power and an enormous stability. The properties of this
new class of laser interferometers have been the subject of extensive research for
several decades.

|Finesse| has been used to support the research on laser interferometers for
gravitational wave detection since 1999 :cite:`Finesse04`, and since 2013 |Finesse| is
continuously
developed as an open source project :cite:`phdBrown`. More about the background
and the early years of |Finesse| (and Pykat) are available in the :ref:`History <history>` 
section.

|Finesse| has become an important tool for the commissioning of Advanced LIGO :cite:`AdvancedLIGO`, 
Advanced Virgo :cite:`AdvancedVirgo` and KAGRA :cite:`KAGRA` and is used for the design of future
detectors such as the Einstein Telescope :cite:`ET2010`. The :ref:`Impact <impact>`
section lists more than 100 documents citing |Finesse|.

|Finesse| version 3 is a complete
re-development, started in 2017, of both the original software, and its eventual wrapper
and utility code Pykat :cite:`Pykat,Brown2020`. The main aim of the redevelopment was to
transform our well tested and established tool with a large active user base into a
modern software package and to make |Finesse| ready for the next 20 years of active
research in laser interferometry.

.. cssclass:: nobg
.. figure:: /images/geo600_birds_eye.jpg
    :align: center

    Bird's eye view of the |GEO600| gravitational wave detector
    near Hannover, Germany. Image courtesy of Harald Lück, Albert
    Einstein Institute Hannover.

Several prototype interferometers had been developed to investigate laser-interferometer
technologies for detecting gravitational waves. This was followed by the work on the
large-scale laser interferometric gravitational wave detectors that led to the first
direct detection of a gravitational wave by the LIGO interferometers in 2015
:cite:`gw150914`. Gravitational-wave astronomy is now an established field in science in
which instrument science remains a major challenge.

The optical systems involved, Fabry-Perot cavities, a Michelson interferometer and
combinations thereof are in principle simple and have been used in many fields of
science for many decades. The sensitivity required for the detection of the expected
small signal amplitudes of gravitational waves, however, has put new constraints on the
design of laser interferometers. The work of the gravitational wave research groups has
led to a new exploration of the theoretical analysis of laser interferometers.
Especially, the clever combination of known interferometers has produced new types of
interferometric detectors that offer an optimised sensitivity for detecting
gravitational waves. We have shown that the models describing the optical system become
very complex even though they are based on simple principles. Consequently, computer
programs have been developed to automate the computational part of the analysis. To
date, several custom-made programs for analysing optical systems are available to the
gravitational wave community, and |Finesse| is one of the most widely used 
tools in this field.





