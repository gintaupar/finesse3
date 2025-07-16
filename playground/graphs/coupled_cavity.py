"""Single cavity test from Alexei's thesis, section 4."""
# %%
import finesse
import finesse.graph
import numpy as np
import matplotlib.pyplot as plt

from finesse.graph.operator_graph import OperatorGraph

# %%

model = finesse.script.parse(
    """
lens L1
m m1
lens L2
m m2
lens L3
m m3
link(L1, m1, L2, m2, L3, m3)
modes(maxtem=1)
"""
)
node_2_index = {n.full_name: i for i, n in enumerate(model.optical_nodes)}
index_2_node = {i: n.full_name for i, n in enumerate(model.optical_nodes)}


def make():
    G = OperatorGraph(len(model.optical_nodes))

    N = len(model.homs)

    for component in list(model.components) + list(model.spaces):
        for name, (i, o) in component._registered_connections.items():
            if i in node_2_index and o in node_2_index:
                G.add_edge(
                    node_2_index[i],
                    node_2_index[o],
                    np.eye(N, dtype=complex),
                    f"{component.name}.{name}",
                )
    return G


carrier_graph = make()


def reduce(self, sources, sinks, keep=None, reductions=None):
    if len(sources) == 0:
        raise Exception("Must specify at least 1 source node")
    if len(sinks) == 0:
        raise Exception("Must specify at least 1 sink node")

    for node in sinks:
        self._validate_node(node)
    for node in sources:
        self._validate_node(node)

    for start in sources + list(range(self.number_of_nodes)):
        N = 1
        while N != 0:
            # Keep going whilst there are edges left to explore
            # by doing a product rule we may open up new paths
            # for reduction
            N = 0
            for a, b in self.output_edges(start):
                if a != b:  # no self loops
                    for c, d in self.output_edges(b):
                        if c != d:  # no self loops
                            if c not in keep:
                                try:
                                    carrier_graph.product_rule(start, c, d)
                                    if reductions is not None:
                                        reductions.append([start, c, d])
                                    N += 1
                                except Exception:
                                    pass
                                    # print("error")


plt.figure(figsize=(5, 5))
pos = carrier_graph.plot(
    alpha=0.2,
)
reductions = []
reduce(carrier_graph, list([2]), [0], [7, 15], reductions)

carrier_graph.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    ignore_nodes=carrier_graph.sink_nodes(),
)

# %%
carrier_graph = make()
G = make()

for i, (a, b, c) in enumerate([(None, None, None)] + reductions):
    plt.figure(figsize=(6, 6))
    G.plot(
        pos=pos,
        alpha=0.2,
    )
    if a is not None:
        carrier_graph.product_rule(a, b, c)

    carrier_graph.plot(
        pos=pos,
        connectionstyle="arc3,rad=0.1",
        ignore_nodes=carrier_graph.sink_nodes(),
    )
    plt.title(f"{a},{b},{c}")
    plt.savefig(f"figures/{(i):03d}.jpeg")
    plt.close()

# %%
