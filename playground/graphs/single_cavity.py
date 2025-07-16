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
m m1
m m2
link(m1.p1, m1.p2, m2.p1)
modes(maxtem=1)
"""
)
node_2_index = {n.full_name: i for i, n in enumerate(model.optical_nodes)}
index_2_node = {i: n.full_name for i, n in enumerate(model.optical_nodes)}

# %%
carrier_graph = OperatorGraph(len(model.optical_nodes))

N = len(model.homs)

for component in list(model.components) + list(model.spaces):
    for name, (i, o) in component._registered_connections.items():
        if i in node_2_index and o in node_2_index:
            carrier_graph.add_edge(
                node_2_index[i],
                node_2_index[o],
                np.eye(N, dtype=complex),
                f"{component.name}.{name}",
            )

plt.figure(figsize=(6, 6))
pos = carrier_graph.plot()

carrier_graph.product_rule(3, 4, 5)
carrier_graph.product_rule(5, 2, 3)
carrier_graph.product_rule(3, 5, 3)
carrier_graph.product_rule(3, 4, 7)
carrier_graph.product_rule(6, 5, 2)
carrier_graph.product_rule(6, 2, 1)

plt.figure(figsize=(6, 6))
carrier_graph.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    alpha=0.2,
)
carrier_graph.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    ignore_nodes=carrier_graph.sink_nodes(),
)

# %%

model = finesse.script.parse(
    """
lens L1
lens L2
m m1
m m2
lens L3
lens L4
link(L1, L2, m1.p1, m1.p2, m2.p1, m2.p2, L3, L4)
modes(maxtem=1)
"""
)
node_2_index = {n.full_name: i for i, n in enumerate(model.optical_nodes)}
index_2_node = {i: n.full_name for i, n in enumerate(model.optical_nodes)}

carrier_graph = OperatorGraph(len(model.optical_nodes))

N = len(model.homs)

for component in list(model.components) + list(model.spaces):
    for name, (i, o) in component._registered_connections.items():
        if i in node_2_index and o in node_2_index:
            carrier_graph.add_edge(
                node_2_index[i],
                node_2_index[o],
                np.eye(N, dtype=complex),
                f"{component.name}.{name}",
            )


def reduce(self, keep=[], reductions=None):
    for node in keep:
        self._validate_node(node)

    total_reductions = 0

    for start in self.source_nodes():
        N = 1  # reset to start loop as no do-while in python...
        while N != 0:
            # Keep going whilst there are edges left to explore
            # by doing a product rule we may open up new paths
            # for reduction
            N = 0
            for a, b in self.output_edges(start):
                if a != b:  # don't follow self loops
                    for c, d in self.output_edges(b):
                        if c != d:  # don't follow self loops
                            if c not in keep and not carrier_graph.has_self_loop(c):
                                # if middle node has a self loop we can reduce this
                                carrier_graph.product_rule(start, c, d)
                                total_reductions += 1
                                if reductions is not None:
                                    reductions.append([start, c, d])
                                N += 1
    return total_reductions


plt.figure(figsize=(6, 6))
pos = carrier_graph.plot(
    alpha=0.2,
)

print(reduce(carrier_graph, keep=[]))

carrier_graph.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    ignore_nodes=carrier_graph.sink_nodes(),
)
