.. include:: /defs.hrst

Basic usage of ports and nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All components in |Finesse| 3 that can be connected together have various :class:`ports
<finesse.components.node.Port>` associated with them, which in turn contain
:class:`.Node` objects representing some optical, or small signal (electrical, or mechanical)
state. A port can only contain one type of node.

Each node has a directionality associated with it for the purposes of builing a
model and desrcibing the flow of information: input, output, or bidirectional. Input
nodes can only be connected to another output node, whereas bidirectional do not care.
Every optical port has an input and output node, usually called `i` and `o`. Electrical
ports usually just have a single input or output node. Mechanical ports are bidirectional nodes.

Ports and nodes have a standard naming and access convention: `component.port.node`.

*Optical Ports* are nearly always called something like `p1`, `p2`, ..., `pN`. Each optical
port always has one input node `i` and one output node `o`. These describe the optical field
going into and out of the component. Optical ports between two different elements can
only be connected to other optical ports via a `Space` element.

*Electrical and Mechanical Ports* contain *Signal Nodes*. These nodes represent electronic,
mechanical, or other small linearised states in the model.
Signal nodes can all be connected together via `Wires` which allow you to form feedback
loops. For example, you can use a `Readout` to convert some optical signal into an electronic
one, which is then filtered, amplified, and fed back to some mechanical state. These signal
states are all small linearised oscillatory states, they do not represent a static effect
in the model. These signal nodes can then be used to inject in signals and then measure
the corresponding response at another signal node.


Accessing ports
```````````````

Ports of a component can be accessed via the component through the names that they were
assigned during construction. For example, a :class:`.Mirror` has two ports (`p1` and
`p2`) which can be grabbed directly:

.. jupyter-kernel:: python3
    :id: ex1

.. jupyter-execute::

    from finesse.components import Mirror

    M1 = Mirror("M1")
    print(M1.p1)
    print(M1.p2)

You can also get a *read-only tuple* of all the ports available at an element using:

.. jupyter-execute::

    print(M1.ports)

or all of the optical nodes:

.. jupyter-execute::

    print(M1.optical_nodes)

or all of the signal nodes:

.. jupyter-execute::

    print(M1.signal_nodes)


What does a port contain?
`````````````````````````

In the example code below we show the main properties of a :class:`.Port` - namely its
type, the nodes that it holds and the component that it is attached to (which is always
a :class:`.Space` between different optical connections):

.. jupyter-kernel:: python3
    :id: ex2

.. jupyter-execute::

    from finesse import Model
    from finesse.components import Mirror

    M1 = Mirror("M1")
    M2 = Mirror("M2")

    model = Model()
    # connect M1 <-> M2 in a model via a Space of length 1m
    model.chain(M1, {"name": "M1_M2", "L":10}, M2)

    print(f"Port M1.p2 name = {M1.p2.name}")
    print(f"Port M1.p2 type = {M1.p2.type}")
    print(f"Port M1.p2 owning component = {M1.p2.component}")
    print(f"Port M1.p2 attached component = {M1.p2.attached_to}")
    print(f"Port M1.p2 nodes = {M1.p2.nodes}")

Accessing nodes --- via a component
```````````````````````````````````

Nodes can be accessed directly through components or via the :class:`.Model` instance
that their owning component is associated with (see next section for details). To access
all the nodes of a component:

.. jupyter-execute::

    print(f"All nodes of M1 = {M1.nodes}")

Or all the optical nodes:

.. jupyter-execute::

    print(f"Optical nodes of M1 = {M1.optical_nodes}")

Get a single optical node of a component by its *direction*:

.. jupyter-execute::

    print(f"Input node of port M1.n1 = {M1.p1.i}")
    print(f"Output node of port M1.n1 = {M1.p1.o}")

Accessing nodes --- via the model
`````````````````````````````````

Nodes play an important role in the :class:`.Model` class as the :attr:`.Node.full_name`
property forms the `node_type` of the underlying directed graph object (stored in
:attr:`.Model.network`). Thus, we use Node instances to report on graph data as well as
perform operations on this graph - see `the networkx DiGraph documentation
<https://networkx.github.io/documentation/stable/reference/classes/digraph.html>`_ for
details on these methods and attributes.

You can also access all the nodes of a given :class:`.NodeType` in a Model instance with
(e.g. for optical nodes):

.. jupyter-execute::

    print(f"All optical nodes in model: {model.optical_nodes}")

Node names are used as keys in the network of a Model to get data on the node itself and
edges connected to the node:

.. jupyter-execute::

    print(model.network.nodes[model.M2.p1.i.full_name])
    print(model.network[model.M2.p1.i.full_name][model.M2.p1.o.full_name])
