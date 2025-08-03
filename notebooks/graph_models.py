# Import necessary modules and functions
from sage.graphs.graph import Graph  # Base graph class in Sage
from sage.graphs.graph_generators import graphs  # For generating predefined graphs
from sage.matrix.constructor import Matrix  # For constructing Sage matrices
from sage.combinat.permutation import (
    Permutation,
)  # For permutations and permutation matrices

# For numerical operations
import numpy as np

# For executing external processes
import subprocess


class FYP_Graph:
    """
    A class to represent a graph for the Final Year Project (FYP) research.
    Includes graph operations, Weisfeiler-Leman (WL) refinement, and coherent rank computation.
    """

    def __init__(self, obj):
        """
        Initialize the graph with a SageMath Graph object.

        Args:
            obj (sage.graphs.graph.Graph): A SageMath Graph object.
        """
        self.obj = obj
        self.dimension = obj.order()
        self.adjacency_matrix = obj.adjacency_matrix()
        self.type_matrix = None
        self.coherent_rank = None
        self.weisfeiler_results = None
        self.P = None
        self.final_config_matrix = None
        self.intervals = None
        self.blocks = None
        self.total = None
        self.points_of_interest = None
        self.config_list = None

    def __repr__(self):
        """
        Returns a string representation of the graph object, showing all attributes.
        """
        return f"Graph object with attributes: {self.__dict__}"

    def get_weisfeiler_results(self):
        """
        Compute and return the Weisfeiler-Leman refinement results.
        Uses an external cpp tool ('weisfeiler') to calculate graph invariants.

        Returns:
            list: Weisfeiler-Leman results split into lines.
                - dimension of graph
                - fibres rank
                - coherent rank
                - array of configurations of the adjacency matrix
        """
        if self.weisfeiler_results is None:
            self.weisfeiler_results = cr(self.obj).split("\n")
        return self.weisfeiler_results

    def get_coherent_rank(self):
        """
        Compute and return the coherent rank of the graph.
        The coherent rank is extracted from the WL refinement results.

        Returns:
            int: Coherent rank of the graph.
        """
        if self.coherent_rank is not None:
            return self.coherent_rank

        if self.weisfeiler_results is None:
            self.get_weisfeiler_results()

        self.coherent_rank = self.weisfeiler_results[2]
        return self.coherent_rank

    def show_graph(self):
        """
        Visualize the graph using SageMath's built-in plot function.
        """
        return self.obj.show()

    def get_adjacency_matrix(self):
        """
        Return the adjacency matrix of the graph.

        Returns:
            sage.matrix: Adjacency matrix of the graph.
        """
        return self.adjacency_matrix

    def get_type_matrix(self):
        """
        Compute and return the type matrix of the graph.
        The type matrix encodes block structure derived from WL refinement.

        Returns:
            sage.matrix.Matrix: The computed type matrix.
        """
        if self.type_matrix is not None:
            return self.type_matrix

        if self.weisfeiler_results is None:
            self.get_weisfeiler_results()

        tm, P1, P2 = type_matrix(self.weisfeiler_results)
        diag_tm = [tm[i][i] for i in range(len(tm))]

        while diag_tm != sorted(diag_tm):
            self.weisfeiler_results = cr(self.obj).split("\n")
            tm, P1, P2 = type_matrix(self.weisfeiler_results)
            diag_tm = [tm[i][i] for i in range(len(tm))]

        self.type_matrix = Matrix(tm)
        self.P = Matrix(P1 * P2)
        matrix_data = list(
            map(int, filter(None, self.weisfeiler_results[3].split(" ")))
        )
        config_matrix = np.array(matrix_data).reshape((self.dimension, self.dimension))
        config_matrix = Matrix(config_matrix)
        self.final_config_matrix = self.P.T * config_matrix * self.P
        return self.type_matrix

    def get_interval(self):
        """
        Compute intervals of consecutive identical elements in the sorted diagonal labels.
        Used to identify blocks.

        Returns:
            list: List of intervals as (start, end) tuples.
        """
        if self.intervals is not None:
            return self.intervals

        if self.weisfeiler_results is None:
            self.get_weisfeiler_results()

        if self.final_config_matrix is None:
            self.get_type_matrix()

        diags = [self.final_config_matrix[i][i] for i in range(self.dimension)]
        self.intervals = get_intervals(diags)
        return self.intervals

    def get_blocks(self, show=False):
        """
        Extract submatrices (blocks) of the adjacency matrix based on intervals.

        Args:
            show (bool): If True, display blocks and their dimensions.

        Returns:
            list: List of submatrices (blocks).
        """
        if self.blocks is None:
            if self.intervals is None:
                self.get_interval()

            permutated_adjacency_matrix = self.P.T * self.adjacency_matrix * self.P

            points_of_interest = set()
            total = []
            for i in range(len(self.intervals)):
                for j in range(len(self.intervals)):
                    total.append((i, j))
                    if self.type_matrix[i][j] != 1:
                        points_of_interest.add((i, j))

            blocks = []
            for i, j in total:
                i_start, i_end = self.intervals[i]
                j_start, j_end = self.intervals[j]
                sub_matrix = permutated_adjacency_matrix[
                    i_start : i_end + 1, j_start : j_end + 1
                ]
                blocks.append(sub_matrix)

            self.blocks = blocks
            self.total = total
            self.points_of_interest = points_of_interest

        if show:
            for matrix, (i, j) in zip(self.blocks, self.total):
                if (i, j) in self.points_of_interest:
                    print(f"Position: {(i, j)}")
                    dims = matrix.dimensions()
                    print(f"Size: {dims}")
                    print(matrix)
                    if dims[0] == dims[1]:
                        visualise_matrix(matrix)
                    print()
        return self.blocks

    def get_config_list(self):
        """
        Return the list of configurations of the adjacency matrix.

        Returns:
            list: List of configurations of the adjacency matrix.
        """
        if self.config_list is not None:
            return self.config_list
        if self.weisfeiler_results is None:
            self.get_weisfeiler_results()

        if self.final_config_matrix is None:
            self.get_type_matrix()

        ret = []
        for i in range(int(self.weisfeiler_results[2].split(" ")[1])):
            temp_matrix = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
            for j in range(self.dimension):
                for k in range(self.dimension):
                    if i == self.final_config_matrix[j][k]:
                        temp_matrix[j][k] = 1
            ret.append(temp_matrix)

        self.config_list = ret
        return ret


    def switch_graph(self, i):
        """
        Perform edge switching for the first vertex and return a new graph.

        Args:
            i: index of vertex to switch edges

        Returns:
            FYP_Graph: A new graph object with switched edges.

        """

        new_g = self.obj.copy()
        vertex = self.obj.vertices()[i]

        for j in self.obj.vertex_iterator():
            if vertex == j:
                continue

            if self.obj.has_edge(vertex, j):
                new_g.delete_edge(
                    vertex, j
                )  # Remove the edge if it exists in original graph
            else:
                new_g.add_edge(
                    vertex, j
                )  # Add the edge if it doesn't exist in original graph

        return FYP_Graph(new_g)
    
    def switch_graph_edges(self, list):
        new_g = self.obj.copy()
        vertices_to_switch = []

        # get the vertices to switch, V
        for i in list:
            vertices_to_switch.append(self.obj.vertices()[i])

        for j in self.obj.vertex_iterator():
            if j in vertices_to_switch:
                continue

            # now j is not in the list of vertices to switch
            # meaning that j is in S/V
            for i in vertices_to_switch:
                if self.obj.has_edge(i, j):
                    new_g.delete_edge(i, j)
                else:
                    new_g.add_edge(i, j)
        

        return FYP_Graph(new_g)


    def delete_vertex(self, i):
        """
        Delete the first vertex and return a new graph without the vertex.

        Args:
            i: index of vertex to delete

        Returns:
            FYP_Graph: A new graph object without the first vertex.
        """
        new_graph = self.obj.copy()
        new_graph.delete_vertex(new_graph.vertices()[i])
        return FYP_Graph(new_graph)
    
    def delete_vertices(self, list):
        new_graph = self.obj.copy()
        vertices_to_delete = []
        for i in list:
            vertices_to_delete.append(self.obj.vertices()[i])
            
        for i in vertices_to_delete:
            new_graph.delete_vertex(i)

        return FYP_Graph(new_graph)

        

    def new_line_graph(self):
        """
        Creates a line graph of the current graph (copy) and return it as a new FYP_Graph object

        Args:
            None

        Returns:
            FYP_Graph: A new graph object that is the line graph of the current graph
        """
        new_graph = self.obj.copy()
        line_graph = new_graph.line_graph()

        return FYP_Graph(line_graph)

    def new_complement(self):
        """
        Creates a complement graph of the current graph (copy) and return it as a new FYP_Graph object

        Args:
            None

        Returns:
            FYP_Graph: A new graph object that is the complement graph of the current graph
        """
        new_graph = self.obj.copy()
        complement_graph = new_graph.complement()

        return FYP_Graph(complement_graph)

    # Add more graph operations here


class OA_Graph(FYP_Graph):
    def __init__(self, m, n):
        """
        Initialize an Orthogonal Array Block Graph.
        """
        obj = self._create_orthogonal_array_block_graph(m, n)
        super().__init__(obj)

    def _create_orthogonal_array_block_graph(self, m, n):
        """
        Create an Orthogonal Array Block Graph using parameters m and n.
        """
        return graphs.OrthogonalArrayBlockGraph(m, n)


class Triangle_Graph(FYP_Graph):
    def __init__(self, n):
        """
        Initialize the line graph of a Complete Graph.
        """
        obj = self._create_complete_graph_line_graph(n)
        super().__init__(obj)

    def _create_complete_graph_line_graph(self, n):
        """
        Create the line graph of a Complete Graph.
        """
        complete_graph = graphs.CompleteGraph(n).line_graph()
        return complete_graph


# Add on more classes of graph here (follow the above 2 conventions)


# Helper functions


def cr(graph):
    """
    returns:
        - dimension of graph
        - fibres rank
        - coherent rank
        - array of configurations of the adjacency matrix
    """
    no_v = len(graph.vertices())
    with open("input.txt", "w") as f:
        f.write(f"order {no_v}\n")
        f.write(f"dimension {no_v}\n")

        # Write each row of the adjacency matrix on a single line
        for row in graph.adjacency_matrix():
            f.write(" ".join(map(str, row)) + "\n")

    return subprocess.run(
        ["./weisfeiler"],
        input=open("input.txt", "r").read(),
        text=True,
        capture_output=True,
        shell=True,
    ).stdout[:-1]


def get_intervals(lst):
    """Return intervals of consecutive identical elements in a sorted list."""
    intervals = []
    n = len(lst)
    start = 0

    for i in range(1, n + 1):
        if i == n or lst[i] != lst[start]:
            intervals.append((start, i - 1))
            start = i

    return intervals


def type_matrix(results):
    """
    Computes the type matrix from the Weisfeiler–Lehman (WL) refinement results.

    This function performs the following steps:
      1) Parses and reshapes the adjacency/configuration matrix from the WL results.
      2) Applies a first permutation (P1) based on the initial diagonal labels.
      3) Constructs an initial type matrix (m x m) using block structures determined by identical diagonal labels.
      4) Reorders blocks by constructing a block-level permutation (P2) to sort the type matrix diagonal.
      5) Applies the block-level permutation to reorder the original adjacency matrix.
      6) Re-condenses the permuted adjacency matrix to produce the final type matrix, ensuring a sorted diagonal.

    Parameters
    ----------
    results : list
        Weisfeiler–Lehman refinement results, containing:
          - results[0] : str
              Dimension (number of vertices) of the adjacency/configuration matrix.
          - results[1] : str
              Space-separated diagonal labels (one per vertex).
          - results[3] : str
              Space-separated flattened adjacency/configuration matrix (n^2 entries).

    Returns
    -------
    final_res : list of lists
        Final (m x m) type matrix with reordered blocks and a sorted diagonal.
    P1 : sage.matrix.matrix_integer_dense.Matrix_integer_dense
        First permutation matrix applied to group the initial diagonal labels.
    P2 : sage.matrix.matrix_integer_dense.Matrix_integer_dense
        Second permutation matrix applied to sort the reordered diagonal labels.
    """

    # ---------------------------
    # 0) Parse the initial matrix
    # ---------------------------

    dimension = int(results[0])
    matrix_data = list(map(int, filter(None, results[3].split(" "))))
    config_matrix = np.array(matrix_data).reshape((dimension, dimension))
    config_matrix = Matrix(config_matrix)

    # ----------------------------------------------
    # 1) First permutation (P1) by original diagonal
    # ----------------------------------------------
    diagonal_labels = list(map(int, filter(None, results[1].split(" ")[:-1])))

    # Sort diagonal labels in ascending order -> build P1
    sorted_diag, perm = zip(
        *sorted((val, idx) for idx, val in enumerate(diagonal_labels))
    )

    # Convert to 1-based (for Sage) and make the permutation matrix
    perm_1based = [i + 1 for i in perm]
    P1 = Permutation(perm_1based).to_matrix()

    # Apply P1^T * config_matrix * P1
    sorted_matrix = P1.T * config_matrix * P1

    # -----------------------------------------
    # 2) Build an initial (m x m) type matrix
    # -----------------------------------------
    intervals = get_intervals(sorted_diag)
    res = []
    for i_start, i_end in intervals:
        row_block = []
        for j_start, j_end in intervals:
            unique_vals = {
                sorted_matrix[r, c]
                for r in range(i_start, i_end + 1)
                for c in range(j_start, j_end + 1)
            }
            row_block.append(len(unique_vals))
        res.append(row_block)

    m = len(res)  # dimension of type matrix
    diag_tm = [res[i][i] for i in range(m)]

    pairs = [(val, idx) for idx, val in enumerate(diag_tm)]
    pairs.sort(key=lambda x: x[0])

    for idx, interval in enumerate(intervals):
        for i in range(interval[0], interval[1] + 1):
            diagonal_labels[i] = pairs[idx][1]

    # Sort diagonal labels in ascending order -> build P2
    sorted_diag, perm = zip(
        *sorted((val, idx) for idx, val in enumerate(diagonal_labels))
    )

    # Convert to 1-based (for Sage) and make the permutation matrix
    perm2 = [i + 1 for i in perm]
    P2 = Permutation(perm2).to_matrix()

    # Apply P2^T * sorted_matrix * P2
    sorted_matrix = P2.T * sorted_matrix * P2

    intervals = get_intervals(
        sorted_diag
    )  # (start, end) blocks for identical diagonal labels
    res = []
    for i_start, i_end in intervals:
        row_block = []
        for j_start, j_end in intervals:
            unique_vals = {
                sorted_matrix[r, c]
                for r in range(i_start, i_end + 1)
                for c in range(j_start, j_end + 1)
            }
            row_block.append(len(unique_vals))
        res.append(row_block)

    return res, P1, P2


def visualise_matrix(sub_matrix):
    """
    Visualize a graph represented by a submatrix of an adjacency matrix.

    This function performs the following steps:
      1) Converts the given submatrix into a contiguous NumPy array for compatibility.
      2) Converts the NumPy array into a SageMath matrix.
      3) Constructs a SageMath graph object from the Sage matrix.
      4) Visualizes the resulting graph using SageMath's built-in visualization tools.

    Parameters
    ----------
    sub_matrix : array-like or matrix
        A submatrix (or full adjacency matrix) representing the graph structure.
        The entries must be numeric (e.g., 0/1 for unweighted graphs, weights for weighted graphs).

    Returns
    -------
    None
        Displays a visual representation of the graph constructed from the adjacency matrix.
    """
    contiguous_matrix = np.ascontiguousarray(sub_matrix)

    # Convert to a Sage matrix
    adj_matrix = Matrix(contiguous_matrix)

    # Create the graph from the adjacency matrix
    G = Graph(adj_matrix)

    # Visualize the graph
    G.show()
