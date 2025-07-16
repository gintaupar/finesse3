.. include:: /defs.hrst

.. _physics:

Physics of Finesse
------------------

The following sections provide information about how the various aspects of an
interferometer simulation are coded within the |Finesse| source code. The analysis of
optical systems described here is based on the principle of superposition of light
fields: a laser beam can be described as the sum of different light fields. The possible
degrees of freedom are:

   - frequency,
   - geometrical shape and position,
   - polarisation.

In the analysis of interferometric gravitational wave detectors, the amplitudes and
frequencies of light fields are of principal interest. The polarisation is neglected in
the analysis given here, but the formalism can in principle be easily extended to
include polarisation also.

This chapter describes the mathematical formalism based on plane waves only. In Chapter
4 the formalism with respect to Hermite-Gauss modes will be given; it is a
straightforward extension of the plane wave analysis and makes use of the methods
described

.. toctree::
   :maxdepth: 1

   plane-waves/index
   higher_order_modes/index
   radiation_pressure/index
   thermal_effects/index
   .. quantum_effects/index
