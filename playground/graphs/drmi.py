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
    ###########################################################################
    ###   Variables
    ###########################################################################
    var Larm 3995
    var Mtm  40
    var itmT 0.014
    var lmichx 4.5
    var lmichy 4.45

    ###########################################################################
    ###   Input optics
    ###########################################################################
    l L0 125
    s l_in L0.p1 prm.p1
    # Power recycling mirror
    m prm T=0.03 L=37.5u phi=90
    s prc prm.p2 bs.p1 L=53


    # Central beamsplitter
    bs bs T=0.5 L=0 alpha=45

    ###########################################################################
    ###   X arm
    ###########################################################################
    s lx bs.p3 itmx.p1 L=lmichx
    m itmx T=itmT L=37.5u phi=90
    s LX itmx.p2 etmx.p1 L=Larm
    m etmx T=5u L=37.5u phi=89.999875

    pendulum itmx_sus itmx.mech mass=Mtm fz=1 Qz=1M
    pendulum etmx_sus etmx.mech mass=Mtm fz=1 Qz=1M

    ###########################################################################
    ###   Y arm
    ###########################################################################
    s ly bs.p2 itmy.p1 L=lmichy
    m itmy T=itmT L=37.5u phi=0
    s LY itmy.p2 etmy.p1 L=Larm
    m etmy T=5u L=37.5u phi=0.000125

    pendulum itmy_sus itmy.mech mass=Mtm fz=1 Qz=1M
    pendulum etmy_sus etmy.mech mass=Mtm fz=1 Qz=1M

    ###########################################################################
    ###   Output and squeezing
    ###########################################################################
    s src bs.p4 srm.p1 L=50.525
    m srm T=0.2 L=37.5u phi=-90

    # A squeezed source could be injected into the dark port
    sq sq1 db=0 angle=90
    s lsqz sq1.p1 srm.p2
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


plt.figure(figsize=(7, 7))
pos = carrier_graph.plot(
    alpha=0.2,
)
carrier_graph.reduce(keep=[5, 31, 17, 25])

carrier_graph.plot(
    pos=pos,
    connectionstyle="arc3,rad=0.1",
    ignore_nodes=list(carrier_graph.sink_nodes()) + [20, 28],
)
# %%
