Examining components
====================

After parsing a component, you can examine the values of the component parameters by
using the :func:`.finesse.element.ModelElement.info` method:

.. jupyter-execute::

    import finesse
    from finesse.components import Laser

    m = finesse.Model()
    m.add(Laser("l1", P=1))
    print(m.l1.info())

If you have defined a parameter value using a reference, this will be reflected in the
info table:

.. jupyter-execute::

    import finesse
    from finesse.components import Laser

    m = finesse.Model()
    m.add(Laser("l1", P=1))
    m.add(Laser("l2", P=m.l1.P.ref))
    print(m.l2.info())

If you want to see the current numerical value for all the parameters, pass
``eval_refs=True``:

.. jupyter-execute::

    import finesse
    from finesse.components import Laser

    m = finesse.Model()
    m.add(Laser("l1", P=1))
    m.add(Laser("l2", P=m.l1.P.ref))
    print(m.l2.info(eval_refs=True))
