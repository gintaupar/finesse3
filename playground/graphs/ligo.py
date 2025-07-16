# %%
import finesse
import finesse.graph
import numpy as np
import matplotlib.pyplot as plt
import finesse.ligo
from finesse.graph.operator_graph import OperatorGraph
from finesse.components.general import InteractionType

# %%
model = finesse.ligo.make_aligo()
# %%
node_2_index = {n.full_name: i for i, n in enumerate(model.optical_nodes)}
index_2_node = {i: n.full_name for i, n in enumerate(model.optical_nodes)}

G = OperatorGraph(len(model.optical_nodes))

N = len(model.homs)

for component in list(model.components) + list(model.spaces):
    for name, (i, o) in component._registered_connections.items():
        if i in node_2_index and o in node_2_index:
            add = True
            if hasattr(component, "R") and component.R.value == 0:
                if component.interaction_type(i, o) == InteractionType.REFLECTION:
                    add = False
            if hasattr(component, "T") and component.T.value == 0:
                if component.interaction_type(i, o) == InteractionType.TRANSMISSION:
                    add = False
            if add:
                G.add_edge(
                    node_2_index[i],
                    node_2_index[o],
                    np.eye(N, dtype=complex),
                    f"{component.name}.{name}",
                )
# %%
plt.figure(figsize=(10, 10))
pos = G.plot(
    alpha=0.2,
    ignore_nodes=list(G.sink_nodes()),
)

G.reduce(
    keep=[
        node_2_index["PRM.p2.o"],
        node_2_index["ITMX.p2.o"],
        node_2_index["ITMY.p2.o"],
        node_2_index["SRM.p1.o"],
    ]
)

G.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    ignore_nodes=list(G.sink_nodes()),
)

# %%
include_nodes = G.evaluation_nodes()
reduced_index = {n: i for i, n in enumerate(include_nodes)}
reduced_index_2_index = {i: n for i, n in enumerate(include_nodes)}

g = OperatorGraph(len(reduced_index))

for n in include_nodes:
    for i, o in G.output_edges(n):
        if o in include_nodes:
            g.add_edge(reduced_index[i], reduced_index[o], np.eye(N, dtype=complex), "")
g.plot()

# eval_nodes = np.arange(g.number_of_nodes)
# M = np.zeros((len(eval_nodes), len(eval_nodes)))
# for i in eval_nodes:
#     for i, o in g.output_edges(i):
#         if o in eval_nodes:
#             M[o, i] = 1

# plt.imshow(M)
# # %%

# %%
