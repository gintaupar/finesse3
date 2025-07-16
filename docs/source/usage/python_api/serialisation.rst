.. include:: /defs.hrst

.. _serialisation_python:

Saving and loading models
=========================

.. note::

    The Python API is still in development and the serialisation format may
    change in the future. It is recommended to only use Pickle for temporary
    storage and not for long term storage.

KatScript offers one way to easily store a particular optical layout which
describes the components, detectors, and analyses to be performed. Whilst this
is convenient and can be stored as a simple ASCII file for long term storage,
it is limited in how complex a model can be described. The Python API offers
full access to make models which use advanced features in |Finesse| as well as
custom user code. In such cases the :class:`finesse.model.Model` objects much be serialised and
deserialised and saved somewhere.

In Python, :mod:`Pickle <pickle>` is the standard library for serialising objects. Pickling
tries to save the object as a binary representation of the object, which can
then be loaded back into memory. The benefit of this is that arbitrarily complex
objects can be saved and loaded. However, there are several downsides to Pickle:

- Pickle is not human readable, so it is not suitable for long term storage
- Pickle is not secure, loading a Pickle file can execute arbitrary code so care
  should be taken on the origin of the Pickle file
- Pickle is not guaranteed to be compatible between different versions of Python
- Pickle is not guaranteed to be compatible between different versions of
  |Finesse|
- Pickle is not guaranteed to be compatible between different versions of
  the Python API
- Pickle is not guaranteed to be compatible between different operating systems,
  such as Windows and Linux, or different architectures, such as x86 and ARM.

With this being said, it is still very useful to be able to save and load
models for temporary short term storage, such as when running an optimisation
and saving the state of the model at each iteration. Or serialising models to
distribute to other processors in a cluster or multi core system.

The interface to load and save models is very simple. The
:meth:`.finesse.model.Model.save` and :meth:`.finesse.model.load` methods can be
used like so:

.. jupyter-execute::

    import finesse
    from finesse import Model
    from finesse.components import Laser, Mirror, Space, Beamsplitter
    from finesse.detectors import PowerDetector

    model = Model()
    l1 = model.add(Laser("l1", P=10))
    pd1 = model.add(PowerDetector("pd1", l1.p1.o))

    model.save("model.pkl")
    loaded_model = finesse.model.load('model.pkl')

    print(model, loaded_model)

.. warning::

    Pickling is not always a fool proof way to save and load models. Sometimes
    the saving works and the loading does not. This is usually due to internal
    issues in  |Finesse| that can be fixed. Please report any issues you have
    to the Gitlab project issue tracker.
