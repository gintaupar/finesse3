.. include:: /defs.hrst

.. _selecting_modes:

Selecting the modes to model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Models in |Finesse| are, by default, created as plane-wave representations of the
underlying configuration. This can be switched to a modal basis by specifying the mode
indices that should be modelled by the system. There are a few key ways in which this
can be done, as outlined in the sections below.

Maximum TEM order
`````````````````

Set a value for the maximum TEM (termed `maxtem` for short) order. Users of previous
version of |Finesse| will be familiar with this option. By setting a value for maxtem,
all of the modes up to and including this order will be included in the model. See an
example below:

.. jupyter-execute::

    import finesse

    model = finesse.Model()
    # set maximum mode order (n+m) to 3
    model.modes(maxtem=3)

    print(model.modes())

Note that setting ``maxtem`` to zero will also switch the model to a modal
basis, where the only mode being modelled then is the TEM00 mode.

You can switch a model back to a plane wave representation by using
:meth:`.Model.switch_off_homs` or passing "off" to
:meth:`finesse.model.Model.modes`:

.. jupyter-execute::

    import finesse

    model = finesse.Model()
    model.modes(maxtem=3)
    model.modes("off")

    print(model.modes())



Selection method
````````````````

One of the main new features of modelling higher-order modes in |Finesse| 3 is the ability
to select the modes you want to model rather than including all modes up to
a specific order with `maxtem`. This can be achieved via varied arguments to the
:meth:`finesse.model.Model.modes` method. The documentation for this
function explains how to use it, so here are a few examples of it in practice.

Selecting even modes up to order 4:

.. jupyter-execute::

    model.modes("even", 4)
    print(model.modes())

Odd modes up to order 3:

.. jupyter-execute::

    model.modes("odd", 3)
    print(model.modes())

Tangential (x) modes up to order 5:

.. jupyter-execute::

    model.modes("x", 5)
    print(model.modes())

Sagittal (y) modes up to maxtem 6:

.. jupyter-execute::

    model.modes("y", 6)
    print(model.modes())

Assigning to modes property
```````````````````````````

The :meth:`finesse.model.Model.modes` function can be used to set the modes
more to include more explicitly. A couple of examples are shown below.

Add specific modes to the model via a list of strings:

.. jupyter-execute::

    model.modes(["00", "11", "22"])
    print(model.modes())

Or a list of tuples (any iterable of length two iterables will work):

.. jupyter-execute::

    model.modes([(0, 0), (2, 1), (1, 3)])
    print(model.modes())

Including additional modes
``````````````````````````

The model class also provides a :meth:`finesse.model.Model.include_modes` method for
inserting additional modes (at the correct, sorted position of the mode indices array)
into the existing modes.

For example, one could select even modes up to order 2 but also include the 11 mode
with:

.. jupyter-execute::

    model.modes("even", 2)
    model.include_modes("11")
    print(model.modes())

A note on specifying mode indices explicitly
````````````````````````````````````````````

If the `modes` argument to :meth:`.Model.modes` (or, indeed,
:meth:`.Model.include_modes` and :meth:`.Model.remove_modes`) is specified as an
iterable of mode indices, then care must be taken to define these modes properly. For
example, if one tries to do::

    model.modes(["010", "123"])

then this would result in an exception being raised due to the ambiguity in the mode
indices when they are larger than 9. For these large HOM indices, one should instead use
tuples (or lists), e.g::

    model.modes([(0, 10), (12, 3)])

to select the HG(0,10) and HG(12,3) modes.

This principle applies equally to selecting detector masks via
:meth:`.MaskedDetector.select_mask`.
