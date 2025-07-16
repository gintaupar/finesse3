# %%
import finesse
import finesse.graph
import numpy as np
import matplotlib.pyplot as plt


from finesse.simulations import BaseSimulation
from finesse.simulations.sparse.KLU import KLUSolver
from finesse.simulations.sparse.solver import SparseSolver
from finesse.simulations.sparse.simulation import SparseMatrixSimulation


import cython
%load_ext Cython

# %%
%%cython -I/Users/ddb/mambaforge/envs/py310/lib/python3.10/site-packages/numpy/core/include
from finesse.simulations.sparse.KLU cimport KLUSolver
from finesse.graph.operator_graph import OperatorGraph
from finesse.simulations.sparse.simulation cimport SparseMatrixSimulation
from finesse.simulations.simulation cimport CNodeInfo
from finesse.simulations.solver cimport MatrixSystemWorkspaces


import numpy as np
cimport numpy as np

cdef class GraphCarrierSolver(KLUSolver):
    """
    Attributes
    ----------
    needs_coupled_matrix : bool
        True if more complicated cycles exist in the graph, essentially means there
        are things like coupled cavities in the system. If True, a larger coupled
        matrix must be made and inverted. Otherwise simple matrix inversions or
        propagations are only needed.
    """
    cdef:
        object graph
        object reduced_graph

        readonly bint plot_graph
        readonly tuple keep_nodes_indices
        readonly tuple input_nodes
        readonly tuple output_nodes
        readonly dict index_2_reduced_index
        readonly dict reduced_index_2_node
        readonly tuple reduced_eval_indices
        readonly tuple reduced_self_loop_nodes

        readonly bint needs_coupled_matrix

        readonly dict self_loop_inversions

    def initialise(self, SparseMatrixSimulation sim):
        self.graph = OperatorGraph(len(self.nodes))
        self.reduced_graph = None
        self.plot_graph = sim.simulation_options.get('plot_graph', False)
        self.keep_nodes_indices = tuple(
            self.node_2_index[name] for name in
            sim.simulation_options.get('keep_nodes', [])
        )

        self.output_nodes = tuple(
            n for _, n in self.nodes.items()
            if len(n.used_in_detector_output) > 0
        )
        self.input_nodes = tuple(
            n for comp in self.input_components() for n in comp.nodes if n in self.nodes
        )

    def make_full_graph(self, connector_workspaces):
        """Using the complete model this makes and populates the self.graph attribute to
        contain all nodes and operators between them. From this a reduced graph will be
        made.

        Parameters
        ----------
        connector_workspaces
            Workspaces for elements connecting other elements that are included in
            this simulation
        """
        cdef Py_ssize_t i_idx, o_idx

        for ws in connector_workspaces:
            component = ws.owner
            if hasattr(component, '_registered_connections'):
                for name, (i, o) in component._registered_connections.items():
                    if i in self.node_2_index and o in self.node_2_index:
                        i_idx = self.node_2_index[i]
                        o_idx = self.node_2_index[o]

                        Neq_i = self._c_node_info[i_idx].nhoms
                        Neq_o = self._c_node_info[o_idx].nhoms

                        self.graph.add_edge(
                            i_idx,
                            o_idx,
                            np.zeros((Neq_o, Neq_i), dtype=complex),
                            f"{component.name}.{name}",
                        )

        if self.plot_graph:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(6, 6))
            self.graph.plot()

    def make_reduced_graph(self):
        """After this step, self.graph will be in a reduced form and self.reduced_graph
        will be a new graph containing the nodes and edges that must be solved for by
        some inversion process."""

        if self.reduced_graph is not None:
            raise Exception("Graph already reduced")

        if self.plot_graph:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(6, 6))
            plt.title("Original vs reduced graph")
            pos = self.graph.plot(
                alpha=0.2
            )

        self.graph.reduce(keep=self.keep_nodes_indices)
        include_nodes = self.graph.evaluation_nodes()
        self.index_2_reduced_index = {n: i for i, n in enumerate(include_nodes)}
        self.reduced_index_2_node = {i: self.index_2_node[n] for i, n in enumerate(include_nodes)}
        self.reduced_graph = OperatorGraph(len(self.index_2_reduced_index))

        if self.plot_graph:
            self.graph.plot(
                pos=pos,
                connectionstyle="arc3,rad=0.1",
                ignore_nodes=list(self.graph.sink_nodes()),
            )

        for n in include_nodes:
            for i, o in self.graph.output_edges(n):
                if o in include_nodes:
                    i_idx = self.index_2_reduced_index[i]
                    o_idx = self.index_2_reduced_index[o]
                    # size of the final operator
                    Neq_i = self._c_node_info[i_idx].nhoms
                    Neq_o = self._c_node_info[o_idx].nhoms

                    self.reduced_graph.add_edge(
                        self.index_2_reduced_index[i],
                        self.index_2_reduced_index[o],
                        np.zeros((Neq_o, Neq_i), dtype=complex),
                        ""
                    )

        if self.plot_graph:
            plt.figure(figsize=(6, 6))
            self.reduced_graph.plot()
            plt.title("New reduced graph")

    def analyse_outputs(self):
        output_plan = {}
        for node in self.output_nodes:
            name = node.full_name
            output_plan[name] = {
                'propagate': [],
                'propagate_inverse': [],
                'propagate_self_loop': [],
            }
            edges = tuple(self.graph.input_edges(self.node_2_index[name]))
            # are the inputs of each of these edges live?
            for i, o in edges:
                if i in self.reduced_eval_indices:
                    # then a matrix inversion result is needed
                    output_plan[name]['propagate_inverse'].append((i, o))
                elif i in self.index_2_reduced_index and self.index_2_reduced_index[i] in self.reduced_self_loop_nodes:
                    # then a matrix inversion result is needed
                    output_plan[name]['propagate_self_loop'].append((i, o))
                else:
                    # a simple propagation is needed
                    if self.index_2_node[i] in self.input_nodes:
                        output_plan[name]['propagate'].append((i,o))

    def analyse(self):
        """Based on the full and reduced graph this will determine what parts of the
        system needs to inverted/sovled or just propagated."""
        if self.reduced_graph is None:
            raise Exception("Reduced graph has not been created")

        # By ignoring self loops, we can see how many nodes are still connected
        # together, as evaluation nodes are those with inputs and outputs.
        self.reduced_eval_indices = self.reduced_graph.evaluation_nodes(ignore_self_loops=True)
        # Nodes that are only self loops and not coupled to others
        self.reduced_self_loop_nodes = tuple(
            set(self.reduced_graph.nodes_with_self_loops) - set(self.reduced_eval_indices)
        )
        if len(self.reduced_eval_indices) == 0:
            # There is either nothing to invert here, a model with no cavities
            self.needs_coupled_matrix = False
        elif len(self.reduced_eval_indices) == 1:
            raise Exception("Unexpected situation, singular node left with no self loop")
        else:
            self.needs_coupled_matrix = True
            # full graph index edges that need including in the matrix
            self.reduced_edges = tuple(
                (
                    self.index_2_reduced_index[i],
                    self.index_2_reduced_index[o]
                )
                for i, o in self.reduced_graph.edges()
            )

        self.self_loop_inversions = {}
        for n in self.reduced_self_loop_nodes:
            self.self_loop_inversions[n] = None

        self.analyse_outputs()

    def define_inverting_matrix(self):
        """"""
        assert(self.needs_coupled_matrix == True)
        node_2_index = {
            self.reduced_index_2_node[n]: self.node_2_index[self.reduced_index_2_node[n]]
            for n in self.reduced_eval_indices
        }
        self._add_matrix_equations(node_2_index)

        for i, o in self.reduced_edges:
            pass

    def get_frequency_coupling_edges(self, connector_workspaces):
        result = {}

        for ws in connector_workspaces:
            component = ws.owner
            if hasattr(component, '_registered_connections'):
                for _, (i, o) in component._registered_connections.items():
                    if i in self.node_2_index and o in self.node_2_index:
                        couples_f, does_f_couple = self._element_couples_frequencies(
                            component, (component.nodes[i], component.nodes[o])
                        )
                        if couples_f:
                            result[
                                (self.node_2_index[i], self.node_2_index[o])
                            ] = does_f_couple

        return result


    cpdef assign_operators(self, connector_workspaces):
        if self._submatrices is not None:
            raise Exception("Submatrices already assigned")

        self._submatrices  = {}
        self._diagonals = {}
        self.connections  = {}
        print(self.get_frequency_coupling_edges(connector_workspaces))
        self.make_full_graph(connector_workspaces)
        self.make_reduced_graph()
        self.analyse()
        # Now we have analysed the system and know what things need inverting,
        # if any, we build the matrix
        if self.needs_coupled_matrix:
            self.define_inverting_matrix()




# %%
model = finesse.script.parse(
    """
laser l1
mod mod1 1M 0
lens L1
lens L2
m m1
m m2
lens L3
lens L4
readout_dc A
link(l1, mod1, L1, L2, m1.p1, m1.p2, m2.p1, m2.p2, L3, L4, A)

modes(maxtem=1)
gauss g1 L1.p1.o w0=1m z=0

pd P m1.p2.o
pd Pt A.p1.i
"""
)

# %%
sim_opts = {
    "simulation_type": SparseMatrixSimulation,
    "carrier_solver": GraphCarrierSolver,
    "signal_solver": KLUSolver,
    "plot_graph": False,
    "keep_nodes": ['m1.p2.o']
}

model.run(simulation_options=sim_opts)
# %%
