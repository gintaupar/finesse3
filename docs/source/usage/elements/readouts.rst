.. include:: /defs.hrst

.. _readouts:

========
Readouts
========

Readouts are a new component in |Finesse| 3 which combine together both a detector and
an optical component. In essence, they represent a photodetector or some other
measurement device that produces an electronic signal. Therefore readout components will
always have one or more optical ports which generate some electrical signal which is
outputted at one or more electronic ports. These components are typically used for
computing transfer functions or AC models of closed loop systems, where the electrical
output can be filtered and then fed back into another component.

DC readout
**********
.. kat:element:: readout_dc

RF photodetector readout
************************
.. kat:element:: readout_rf
