.. include:: /defs.hrst
.. _ports:

The port and node system
************************

One of the key structural improvements made to |Finesse| is the introduction of a
flexible port and node system. This system underpins the code structure of |Finesse| as
all couplings and interactions are performed via ports and nodes of components.

A :class:`.Port` in |Finesse| 3 is fundamentally similar to the definition of a node in
|Finesse| 2; i.e. it represents a point of connections between components in a model.
Each port then has single or multiple `nodes`. There are three types of ports/nodes in
|Finesse| 3, Optical, Electrical, and Mechanical (see :class:`.NodeType`) - |Finesse| 2
only had optical which limited the type of simulations that were possible. Nodes
(:class:`.Node`) at each port always have the same type as the port itself. Each node
represents a state of the optomechanical system we are modelling. Mathematically this
means that each node corresponds to a particular linear equation in the linear set of
equations that describe the behaviour of the system we are modelling.


.. toctree::
    :maxdepth: 1

    basic_usage
