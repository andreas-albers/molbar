import networkx as nx
import numpy as np
import numbers
from molbar.io.filereader import FileReader
from molbar.molgraph.nodes.nodes import Nodes, node_data
from molbar.molgraph.edges.edges import Edges
from molbar.molgraph.edges.coordination import Coordination
from molbar.molgraph.edges.bond_order import BondOrder
from molbar.molgraph.cycles import Cycle
from molbar.molgraph.nodes.vsepr import VSEPR, ideal_soaps, ideal_geometries
from molbar.molgraph.edges.dihedral import Dihedral
from molbar.molgraph.repulsion import Repulsion
from molbar.molgraph.nodes.priorities import Priorities


class MolGraph(
    nx.Graph,
    FileReader,
    Nodes,
    Edges,
    Coordination,
    BondOrder,
    Priorities,
    Cycle,
    Repulsion,
    VSEPR,
    Dihedral,
):
    """
    MolGraph is a subclass of networkx.Graph. It inherits all methods from networkx.Graph and adds additional methods for the construction of molecular graphs.
    """

    def __init__(self, total_charge=0, data=None, **attr):
        """
        Constructs MolGraph object.
        Args:
            data (optional): Additional data needed for the construction fo the Networkx graph. Defaults to None.
        """
        super().__init__(data, **attr)
        assert isinstance(total_charge, numbers.Integral), "Total charge must be an integer."
        self.total_charge = total_charge
        self.constraints = {}
        self.unified_distance_matrix = np.array([])
        self.cycle_nodes = []
        self.visible_repulsion_nodes = []
        self.fragment_data = {}
        self.is_2D = False

    def from_file(self, filepath: str) -> None:
        """
        Reads a file and constructs a MolGraph object from coordinates.

        Args:
            filepath (str): Path to file containing coordinates. Either .xyz, .coord, .sdf or .mol file.

        Raises:
            XYZNotFound: If the file does not exist.
            NotXYZFormat: If the input geometry is not in .xyz format.
        """

        self.filepath = filepath
        # Load tabulated chemical data from data file such as electronegativity, atomic radius, etc.
        self.element_data = node_data
        self.ideal_soaps = ideal_soaps
        self.ideal_geometries = ideal_geometries
        # Read file and return number of atoms, coordinates and elements
        n_atoms, coordinates, elements = self.read_file()
        nodes = list(range(n_atoms))
        for node, coordinate, element in zip(nodes, coordinates, elements):
            # Get specific element data for each atom from tabulated element_data
            element_data = self._get_data_for_element(element)
            # Add node to graph
            self.add_node(
                node,
                coordinates=coordinate,
                elements=element,
                visible=True,
                **element_data,
            )

    def from_coordinates(self, coordinates: np.ndarray, elements: np.ndarray, edges=None, bond_orders=None, formal_charges=None) -> None:
        """
        Constructs a MolGraph object from coordinates given by np.ndarray and elements given by np.ndarray.

        Args:
            coordinates (np.ndarray): Coordinates of the molecule of shape (n_atoms, 3).
            elements (np.ndarray): Elements of the molecule of shape (n_atoms,).

        Raises:
            ValueError: If coordinates and elements do not have the same length.
        """

        if len(coordinates) != len(elements):
            raise ValueError("geometry and atoms must have the same length")
        nodes = list(range(len(coordinates)))
        # Load tabulated chemical data from data file such as electronegativity, atomic radius, etc.
        self.element_data = node_data
        # Load ideal SOAP values from data file.
        self.ideal_soaps = ideal_soaps
        self.ideal_geometries = ideal_geometries
        allowed_elements = list(self.element_data["atomic_numbers"].keys())
        if all([(isinstance(element, numbers.Integral)) for element in elements]):
            elements = [allowed_elements[element - 1] for element in elements]
        for node, coordinate, element in zip(nodes, coordinates, elements):
            # Get specific element data for each atom from tabulated element_data
            element_data = self._get_data_for_element(element)
            # Add node to graph
            self.add_node(
                node,
                coordinates=coordinate,
                elements=element,
                visible=True,
                **element_data,
            )
        self.custom_edges = edges
        self.custom_bond_orders = bond_orders
        if formal_charges is not None:
            self.formal_charges = True
            self.add_node_data("formal_charge", formal_charges)
        else:
            self.formal_charges = False

    def return_n_atoms(self, include_all=False) -> int:
        """
        Returns the number of atoms in the molecule.

        Args:
            include_all (bool, optional): If all nodes and not only the visible nodes should be considered. The default value is False

        Returns:
            int: _description_
        """
        if include_all:
            return self.number_of_nodes()
        return len(
            [node for node, data in self.nodes(data=True) if data.get("visible")]
        )

    def return_cn_matrix(self, include_all=False) -> np.ndarray:
        """
        Returns the coordination matrix/adjacency matrix of the molecule.

        Args:
            include_all (bool, optional): If all nodes and not only the visible nodes should be considered. The default value is False

        Returns:
            np.ndarray: Coordination matrix/adjacency matrix of the molecule.
        """
        if include_all:
            return nx.to_numpy_array(self)
        else:
            nodes = self.return_visible_nodes()
            cn_matrix = np.zeros((len(nodes), len(nodes)))
            for node1, node2, data in self.edges(data=True):
                if data["visible"]:
                    if (node1 in nodes) and (node2 in nodes):
                        cn_matrix[nodes.index(node1), nodes.index(node2)] = 1
                        cn_matrix[nodes.index(node2), nodes.index(node1)] = 1
            return cn_matrix

    def get_degree(self, include_all=False) -> np.ndarray:
        """
        Returns the degree of each node in the molecule.

        Args:
            include_all (bool, optional): If all nodes and not only the visible nodes should be considered. The default value is False

        Returns:
            np.ndarray: Degree of each node in the molecule.
        """
        if include_all:
            return np.array([degree for node, degree in self.degree()])
        subgraph = nx.Graph()
        subgraph.add_nodes_from(self.return_visible_nodes())
        subgraph.add_edges_from(self.return_visible_edges(nodes_visible=True))
        return np.array([degree for node, degree in subgraph.degree()])

    def connected_components(self) -> list:
        """
        Returns the connected nodes of the molecule. Nodes are connected if there is a visible edge between them. Further, nodes need to be visible.

        Returns:
            list: List of connected nodes.

        """
        visible_edges = self.return_visible_edges(nodes_visible=True)
        subgraph = nx.Graph()
        subgraph.add_nodes_from(self.return_visible_nodes())
        subgraph.add_edges_from(visible_edges)
        connected_components = list(
            [list(fragment) for fragment in nx.connected_components(subgraph)]
        )
        return connected_components

    def rigid_fragmentation(self, include_all=False) -> None:
        """

        Fragmentates the molecule into rigid parts by considering cycles and multiple bonds as rigid fragments.
        Adds attribute fragment_id to the nodes that cluster nodes into fragments.

        Args:

            include_all (bool, optional): If all nodes and not only the visible nodes should be considered. The default value is False

        """
        if include_all:
            self.set_nodes_visible(include_all=True)
        self.set_edges_visible(include_all=True)
        degrees = self.get_degree()
        visible_nodes = np.array(self.return_visible_nodes())
        unvisible_nodes = visible_nodes[np.where(degrees == 1)[0]]
        self.set_nodes_unvisible(unvisible_nodes=unvisible_nodes)
        self.set_edges_unvisible(attributes=["rigid", "cycle_id"], values=[False, None])
        self.set_edges_unvisible(
            attributes=["stick", "cycle_id"],
            values=[True, None],
            default_all_visible=False,
        )
        fragments = self.connected_components()
        self.set_edges_visible(include_all=True)
        fragment_id = 0
        sorted_fragments = sorted(fragments, key=lambda x: len(x), reverse=True)

        for fragment in sorted_fragments:
            self.set_nodes_visible(visible_nodes=fragment)
            self.add_node_data(
                attribute="fragment_id", new_data=[fragment_id] * len(fragment)
            )
            fragment_id += 1
        self.set_nodes_visible(include_all=True)

    def set_all_visible(self):
        """
        Sets all nodes and edges to visible.
        """
        self.set_nodes_visible(include_all=True)
        self.set_edges_visible(include_all=True)
