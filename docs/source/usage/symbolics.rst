.. include:: /defs.hrst

.. _symbolics:

=========
Symbolics
=========

|Finesse| 3 introduces a new symbolic math feature that can be used to easily
define mathematical expressions, transformations, and build mode complex models.
For those coming from |Finesse| 2, this engine has replaced
commands such as `func`, `var`, `put`, and `set` with a more powerful framework.

The symbolic engine provided in |Finesse| is custom written and, whilst similar
to packages such as Mathematica and Sympy, is not as feature complete
or powerful. The symbolics have both a KatScript and Python interface to provide
maximum flexibility. Here we will detail some of the main use cases for symbolic
expressions in |Finesse|.

Referencing parameters
~~~~~~~~~~~~~~~~~~~~~~

The most common use case for symbols is to quickly allow you to link different
elements' parameters together. Take this simple example of two lasers, using
KatScript, we can quickly make `l2` always have the same power as `l1` by
simply writing:

.. jupyter-execute::

    import finesse
    model = finesse.script.parse("""
    l l1 P=1
    l l2 P=l1.P
    """)

    model.l1.P, model.l2.P

We can see that `l1.P` has a value of 1 whereas `l2.P` has a value of `l1.P`.
To explore this more, we can print the value of the parameter:

.. jupyter-execute::

    model.l2.P.value, type(model.l2.P.value)

We can see this is a `Symbolic` object which is a `ParameterRef` type, which
means it is a symbol that references the current value of `l1.P`. If we want to
get the current value of a symbol we can simply evaluate it:

.. jupyter-execute::

    model.l2.P.value.eval()

and we see we get a value of `1`. If `l1.P` changes and we evaluate it again, we
see we get the new value.

.. jupyter-execute::

    model.l1.P = 10
    model.l2.P.value.eval()

We can also achieve this using the Python API:

.. jupyter-execute::

    model = finesse.script.parse("""
    l l1 P=1
    l l2 P=2
    """)

    print(model.l1.P.ref)
    model.l2.P = model.l1.P.ref
    print(model.l2.P)

Note that here we can take any `ModelParameter` and get a symbolic reference to
it by using the `.ref` attribute. These references can be used to make more
complicated expressions if you need to:

.. jupyter-execute::

    import numpy as np
    expression = 10 * model.l2.P.ref + np.cos(model.l1.P.ref)

    print(expression)
    print(expression.eval())

Note above we can easily do standard mathematical operations on symbols in
Python and also in KatScript:

.. jupyter-execute::

    model = finesse.script.parse("""
    l l1 P=1
    l l2 P=10*cos(l1.P)
    """)

KatScript will always assume you are making a symbolic operation, there is no need to
specify `.ref` in KatScript. The benefit of this behaviour is that when you run simulations
the second laser power will automatically follow the expressions when the power of `l1`
is changed.

Another common use case is if you wanted to vary the reflectivity or transmission
of a mirror. For example, you can sweep the transmission of a mirror and always
keep the reflectivity the correct value:

.. jupyter-execute::

    model = finesse.script.parse("""
    l l1 P=1
    m m1 R=1-m1.T T=0 L=0
    link(l1, m1)

    xaxis(m1.T, lin, 0, 1, 100)
    """)

    # or with Python
    model.m1.R = 1 - model.m1.T.ref


Variables vs ParameterRef
~~~~~~~~~~~~~~~~~~~~~~~~~
There are two types of variables in the |Finesse| symbolics, `ParameterRef` and
`finesse.symbol.Variable`. The former is described above and is a symbol whose value is attached
to some model elements parameter. It is also possible to use a simple variable
which is not.

Broadly they both work in the same way. The difference is mostly when evaluating
some expression. A `ParameterRef` will evaluate to the current value of the parameter,
whereas a `Variable` will not. The only way to give a variable a value to numerically
evaluate an expression is to substitute one in.


Substitution
~~~~~~~~~~~~
Substitution involves mapping a symbol with another expression or numeric value.
There are two ways to do this in |Finesse|. The first is a pure substitution:

.. jupyter-execute::

    a = finesse.symbols.Variable('a')
    b = finesse.symbols.Variable('b')
    y = a+b
    y.substitute({a: b+1})

Alternatively you can also evaluate and substitute at the same time

.. jupyter-execute::

    a = finesse.symbols.Variable('a')
    b = finesse.symbols.Variable('b')
    y = a+b
    y.eval(subs={a: b+1})

Note that this is the same result when using `Variables`, as they are not evaluated
to a specific value. Note that when using parameter references this behaviour will
not be the same:

.. jupyter-execute::

    model = finesse.script.parse("""
    var a 1
    var b 2
    """)
    y = model.a.ref + model.b.ref
    y.eval(subs={a: b+1})

`a` is replaced with `b+1` then evaluated as `b = 2`. Using `substitute` you can
avoid the evaluation of the parameter references:

.. jupyter-execute::

    y.substitute({a: b+1})

You can also override `b` in the substitution if needed

.. jupyter-execute::

    y.eval(subs={a: b+1, b: 3})

even with arrays

.. jupyter-execute::

    y.eval(subs={a: b+1, b: np.array([1,2,3,4])})

Some substitutions/mappings are not unique and if applied recursively in different
orders will result in different results. |Finesse| again does not apply the mappings
recursively to avoid this. Therefore a substitution like below has a definite answer
but still has a `b` term in it, despite the mapping `b=>a` being applied.

.. jupyter-execute::

    (a+b).substitute({a: a+b, b: a})


Making a function from an expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is possible to convert a symbolic expression into some evaluable function.
This process is called "lambdifcation", or generating some anonymous
`lambda` function that can be called. This is particular useful if you have some
expression, say from the ABCD tracing, and you want to run it in some optimiser.

.. jupyter-execute::

    model = finesse.script.parse("""
    var a 1
    laser L1 P=10
    """)
    y = 2 * model.a.ref * model.L1.P.ref

    print("y =", y, "=", y.eval())

We can convert this expression into a function we can just call.

.. jupyter-execute::

    f = y.lambdify(model.a.ref, model.L1.P.ref)
    print(f(1, 10))
    print(f(2, 10))
    print(f(2, 20))

If you need something that evaluates faster, you can also look into using `numba`
JIT to compile the generated function into something that evaluates faster.

Or if you want a function that is just using one parameter and the rest are the



.. jupyter-execute::

    f = y.lambdify(model.a.ref, model.L1.P.ref)
    print(f(1, 10))
    print(f(2, 10))
    print(f(2, 20))

Preserving intent
~~~~~~~~~~~~~~~~~

One feature of the |Finesse| symbolics which differs from many other Symbolic
packages is that there is no automated simplification. The key idea here is that
we want to record exactly what the user originally intended. This for example,
is useful for keeping track of factors of two or pi, etc, which when cancelled or
simplified maybe lost.

.. jupyter-execute::

    model = finesse.script.parse("""
    l l1 P=2*2/2-pi+pi
    """)

    model.l1.P


two-arg vs N-arg operators
~~~~~~~~~~~~~~~~~~~~~~~~~~

In our symbolic computation everything is broken down into operations and their
arguments. These in turn then form a tree structure. Two varieties exist, binary-trees
or operators accepting only two arguments, and a general tree structure, where
operators can have N-arguments. The |Finesse| symbolic feature actually support
both. The binary-tree structure is used by default, and generates an expression
tree that fully preserves the original intent of the user. The binary type of tree
is not ideal for performing symbolic simplification procedures though. Therefore
there is also a generic tree structure option as well. General |Finesse| users
should not really be able to tell which is being used at any particular time,
nor should it make large differences in the behaviour of the code. It may be
important for users performing more complex analysis though. More details on
symbolic computation can be read up on in :cite:`cohen2003computer`, which the
|Finesse| implementation is based upon.

In general the binary tree structure is apparent when you print longer expressions
(You can also check the symbols `_is_narg_expression_tree` attribute):

.. jupyter-execute::

    a = finesse.symbols.Variable('a')
    y = 2+a+2+a+2*a*a/2
    y

You should note that the brackets collect together always two terms in each operation.
you can convert between the two tree structures using `.to_nary_add_mul()` and
`.to_binary_add_mul()`. These are named as such as the general tree structure
only really applies to the add and multiplication operator.

.. jupyter-execute::

    a = finesse.symbols.Variable('a')
    y = 2+a+2+a + 2*a*a/2
    y.to_nary_add_mul()

During the above conversion processes some obvious simplification happens. Negation
and division operators also do not exist in the N-argument trees. Those can always
be represented by multiply and power operators which makes simplification and
comparing two expressions significantly easier.

.. jupyter-execute::

    a = finesse.symbols.Variable('a')
    y = 1/a
    print(y)
    print(y.to_nary_add_mul())

This however highlights why the N-argument option is not the default, as we have lost
the intent of a division. This has subtle numeric implications as a floating point
division is not numerically identical to a power operation in all cases.

Basic simplification
~~~~~~~~~~~~~~~~~~~~

In cases where you do want some basic simplification of |Finesse| symbolics there
are only a few options, but they should cover most basic algebraic requirements.
As needs evolve more may be added:

* `expand` - expand operators
* `collect` - collect together like terms

.. jupyter-execute::

    model = finesse.script.parse("""
    var b 1
    l l1 P=b+2-b+(b+b)
    """)

    print(model.l1.P.value)
    print(model.l1.P.value.expand())
    print(model.l1.P.value.collect())

The simplification should also work on functions:

.. jupyter-execute::

    (np.cos(a) - np.cos(a)).collect()
