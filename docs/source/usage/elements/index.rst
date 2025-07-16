.. include:: /defs.hrst

.. _elements:

Model elements
==============

A model in |Finesse| is built up of multiple elements. Each element represents some
real-world components such as mirrors, lenses, photodetectors, suspensions, control
filters, etc. This section contains descriptions of each element the physics behind
|Finesse|, with equations provided that match the expressions used within the code
itself.

Every model element has a `name` given to it. These can all be referenced in KatScript
and through the Model Python object.

Components and detectors
************************

Components and detectors are model elements that can be connected together to build your
model.

.. toctree::
    :maxdepth: 2

    optics
    sources
    connectors
    detectors
    readouts
    control_and_filtering
    mechanical

Other elements
**************

There are also other elements which do not represent physical components that can be
used in a model.

.. toctree::
    :maxdepth: 2

    variables
    gaussian_beams
