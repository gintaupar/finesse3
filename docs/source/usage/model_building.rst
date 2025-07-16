.. include:: ../defs.hrst

Building and interacting with models
====================================

Before you can perform simulations of your optical system you must first build a model
of it. |Finesse| includes a variety of tools and features to make it straight forward to
build models ranging from simple single Fabry-Perot cavities to full gravitational wave
interferometers.

We offer two ways in which you can build a model. Firstly using KatScript, which is a
declarative domain specific language we have developed for building models of
interferometers. This removes the complexity of dealing with Python and keeps you
focused on the physics of your model. Secondly, we offer a full Python application
programming interface (API) that can also be used for building a model for those wishing
to use a more programmatic object-orientated interface.

Neither is better or worse, however most of the examples in the documentation will be
shown using KatScript. Any feature in KatScript has an equivalent in the Python API,
however not all features in |Finesse| are guaranteed to have an equivalent KatScript
option. This will likely be true for newer or more complex features that users would
only really use from a Python interface.

The Python API section also details how you can interact with your model once it has
been built. It allows you to easily amend and update the model and parameters.


.. toctree::
    :maxdepth: 3

    kat_script
    python_api/index.rst
    symbolics
