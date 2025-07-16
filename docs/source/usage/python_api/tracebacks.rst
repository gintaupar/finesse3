Tracebacks in Jupyter notebooks
===============================


Andreas Freise, 16/04/2020

The following demonstrates how Finesse handles parse errors. When run from IPython or Jupyter the full traceback is suppressed and only the parser error itself is shown. When run from Python it always returns the full traceback.

You can print the latest full traceback with `finesse.tb()`.

If you want to always see the full traceback, for example, during debugging, you can set this with `finesse.show_tracebacks()`

We use a very simple example model for the demonstration:

.. jupyter-execute::

    from finesse.env import is_interactive

    print(is_interactive())

.. jupyter-execute::

    import finesse

    code = """
    l laser P=1
    pd tmp laser.p1.o
    """

    kat = finesse.Model()
    kat.parse(code)

A syntax error in the katscript gives a parser error:

.. jupyter-execute::
    :raises:

    code = """
    l1 laser P=1
    pd tmp laser.p1.o
    """

    kat = finesse.Model()
    kat.parse(code)

.. We can print the full traceback:

.. .. jupyter-execute::

..     finesse.tb(colors=False)

.. To switch tracebacks on globally:

.. .. jupyter-execute::

..     finesse.show_tracebacks(True)

.. jupyter-execute::
    :raises:

    code = """
    l1 laser P=1
    pd tmp laser.p1.o
    """

    kat = finesse.Model()
    kat.parse(code)

`finesse.show_tracebacks()` can also switch tracebacks off again:

.. jupyter-execute::

    finesse.show_tracebacks(False)
