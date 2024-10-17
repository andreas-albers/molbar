import numpy as np
import networkx as nx
from molbar.exceptions.error import FileNotFound
from collections.abc import Iterable


class Edges:

    def define_edges(self, include_all=False):
        """
        Defines the edges of the molecule based on the coordinates of the atoms.

        Args:
            include_all (bool, optional): If True, all nodes are included. Defaults to False so that only visible nodes are included.
        """
        # Use custom edges if provided.
        if self.custom_edges is not None:
            rij0 = self.return_node_data(attribute="rcov", include_all=include_all)
            elements = self.return_node_data(attribute="elements", include_all=include_all)
            rij0[np.where(elements == "H")[0]] += 0.3
            for edge in self.custom_edges:
                dij = rij0[edge[0]] + rij0[edge[1]]
                self.add_edge(edge[0], edge[1], bo=1, visible=True, dij=dij, rigid=False)
        # else, get the edges from the molecule coordinates.
        else:
            # edges, rij0 = self._get_edges_d3(include_all=include_all)
            edges, rij0 = self._get_edges_gfnff(include_all=include_all)
            # Add the edges to the graph.
            for edge in edges:
                dij = rij0[edge[0]] + rij0[edge[1]]
                self.add_edge(edge[0], edge[1], bo=1, visible=True, dij=dij, rigid=False)
        self._add_terminal_attribute()
        self._add_stick_attribute()

    def return_edges(
        self, attributes=None, values=None, nodes_visible=False, include_all=False
    ) -> list:
        """
        Returns the edges of the molecule with certain attributes and values.

        Args:
            attributes (Union(list, str), optional): list of attributes the edges should be checked for. Defaults to None.
            values (Union(list, str), optional): list of values for the attributes the edges should have. Defaults to None.
            nodes_visible (bool, optional): If only visible nodes should be considered. Defaults to False.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.

        Raises:
            ValueError: If attributes and values are not provided as str or list/np.array or if they have different lengths.

        Returns:
            list: list of tuples containing edges with the attributes and values.
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        elif (isinstance(attributes, Iterable) and isinstance(values, Iterable)) and (
            len(attributes) != len(values)
        ):
            raise ValueError("attributes and values must have the same length")
        elif not attributes and not values and not include_all:
            raise ValueError(
                "You must provide either (attributes and values) as str (one attribute) or as lists/np.array (several attributes) or visible_nodes as a list."
            )
        if include_all:
            return self.edges
        else:
            visible_nodes = self.return_visible_nodes()
            if nodes_visible:
                visible_edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if all([node in visible_nodes for node in (node1, node2)])
                    and (data.get("visible"))
                ]
            else:
                visible_edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if (data.get("visible"))
                ]
            return [
                (node1, node2)
                for node1, node2, data in self.edges(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
                and ((node1, node2) in visible_edges)
            ]

    def return_edge_data(
        self, attribute: str, nodes_visible=False, include_all=False
    ) -> np.ndarray:
        """
        Returns the data of the edges.

        Args:
            attribute (str): Data of the edges to be returned.
            nodes_visible (bool, optional): If only visible nodes should be considered. Defaults to False.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.

        Returns:
            np.ndarray: Data of the edges.
        """
        if include_all:
            edge_data = [
                data.get(attribute) for node1, node2, data in self.edges(data=True)
            ]
        else:
            if nodes_visible:
                visible_nodes = self.return_visible_nodes()
                edge_data = [
                    data.get(attribute)
                    for node1, node2, data in self.edges(data=True)
                    if all([node in visible_nodes for node in (node1, node2)])
                    and (data.get("visible"))
                ]
            else:
                edge_data = [
                    data.get(attribute)
                    for node1, node2, data in self.edges(data=True)
                    if (data.get("visible"))
                ]
        try:
            return np.array(edge_data)
        except ValueError:
            return edge_data

    def return_visible_edges(self, nodes_visible=False, include_all=False) -> list:
        """
        Returns the visible edges of the molecule.

        Args:
            nodes_visible (bool, optional): If also nodes and addition to the edges should be visible. Defaults to False.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.

        Returns:
            list: list of visible edges.
        """
        if include_all:
            return self.edges
        else:
            if nodes_visible:
                visible_nodes = self.return_visible_nodes()
                edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if all([node in visible_nodes for node in (node1, node2)])
                    and (data.get("visible"))
                ]
            else:
                edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if (data.get("visible"))
                ]
        return edges

    def return_unvisible_edges(self, nodes_visible=False, include_all=False) -> list:
        """
        Returns the unvisible edges of the molecule.

        Args:
            nodes_visible (bool, optional): If also nodes and addition to the edges should be visible. Defaults to False.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.

        Returns:
            list: list of unvisible edges.
        """
        if include_all:
            return self.edges
        else:
            if nodes_visible:
                visible_nodes = self.return_visible_nodes()
                edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if all([node in visible_nodes for node in (node1, node2)])
                    and (not data.get("visible"))
                ]
            else:
                edges = [
                    (node1, node2)
                    for node1, node2, data in self.edges(data=True)
                    if (not data.get("visible"))
                ]
        return edges

    def set_edges_visible(
        self, attributes=None, values=None, visible_edges=None, include_all=False
    ) -> None:
        """
        Sets the edges of the molecule with certain attributes and values visible.

        Args:
            attributes (Union(list, str), optional): list of attributes the edges should be checked for. Defaults to None.
            values (Union(list, str), optional): list of values for the attributes the edges should have. Defaults to None.
            visible_edges (list, optional): list of edges that should be set visible as alternative to provide attributes and values. Defaults to None.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.

        Raises:
            ValueError: if attributes and values are not provided as str or list/np.array or if they have different lengths.
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        elif (isinstance(attributes, Iterable) and isinstance(values, Iterable)) and (
            len(attributes) != len(values)
        ):
            raise ValueError("attributes and values must have the same length")
        elif isinstance(visible_edges, Iterable) and not isinstance(visible_edges, str):
            pass
        elif not attributes and not values and not include_all:
            raise ValueError(
                "You must provide either (attributes and values) as str (one attribute) or as lists/np.array (several attributes) or visible_nodes as a list."
            )

        if include_all:
            for edge in self.edges:
                self.add_edge(edge[0], edge[1], visible=True)
            return
        if visible_edges is not None:
            for edge in self.edges():
                self.add_edge(edge[0], edge[1], visible=False)
            for edge in visible_edges:
                self.add_edge(edge[0], edge[1], visible=True)
        elif attributes and values:
            visible_edges = [
                (node1, node2)
                for node1, node2, data in self.edges(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
            ]
            for edge in self.edges():
                self.add_edge(edge[0], edge[1], visible=False)
            for edge in visible_edges:
                self.add_edge(edge[0], edge[1], visible=True)
        else:
            raise ValueError(
                "Either (attribute and value) or visible_nodes must be provided."
            )

    def set_edges_unvisible(
        self, attributes=None, values=None, include_all=False, default_all_visible=True
    ):
        """
        Sets the edges of the molecule with certain attributes and values unvisible.

        Args:
            attributes (str, optional): Sets the edges with this attribute unvisible if they are equal to the provided value. Defaults to None.
            values (str, optional): Sets the edges with this value for the attribute unvisible. Defaults to None.
            include_all (bool, optional): if all nodes should be considered. Defaults to False.
            default_all_visible (bool, optional): wether to set all edges to visible, before setting the edges with the provided attributes and values unvisible. Defaults to True.
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        if include_all:
            for edge in self.edges():
                self.add_edge(edge[0], edge[1], visible=False)
            return
        if default_all_visible:
            self.add_edges_from(self.edges, visible=True)
        unvisible_edges = [
            (node1, node2)
            for node1, node2, data in self.edges(data=True)
            if all(
                [
                    data.get(attribute) == value
                    for attribute, value in zip(attributes, values)
                ]
            )
        ]
        for edge in unvisible_edges:
            self.add_edge(edge[0], edge[1], visible=False)

    def get_graph_distance_matrix(self, include_all=False) -> dict:
        """
        Returns the distance matrix of the molecule as a dictionary#

        Args:
            include_all (bool, optional): If all nodes should be considered. Defaults to False.
        """
        if include_all:
            return nx.floyd_warshall_numpy(self)
        subgraph = nx.Graph()
        visible_nodes = self.return_visible_nodes()
        subgraph.add_nodes_from(visible_nodes)
        subgraph.add_edges_from(self.return_visible_edges())
        matrix = nx.floyd_warshall_numpy(subgraph)
        return matrix

    def get_all_shortest_paths(self, include_all=False) -> dict:
        """
        Returns the shortest graph paths between all nodes of the molecule based on the Dijkstra algorithm from networkx.

        Args:
            include_all (bool, optional): If all nodes should be considered. Defaults to False.

        Returns:
            dict: Dictionary containing the shortest paths between all nodes.
        """
        if include_all:
            visible_nodes = self.return_visible_nodes(include_all=include_all)
            return {
                source: {
                    target: nx.shortest_path(self, source, target)
                    for target in visible_nodes
                    if source != target
                }
                for source in visible_nodes
            }
        else:
            visible_nodes = self.return_visible_nodes(include_all=include_all)
            subgraph = nx.Graph()
            subgraph.add_nodes_from(visible_nodes)
            subgraph.add_edges_from(self.return_visible_edges())
            return {
                source: {
                    target: nx.shortest_path(subgraph, source, target)
                    for target in visible_nodes
                    if source != target
                }
                for source in visible_nodes
            }

    def get_shortest_path(self, source: int, target: int, include_all=False) -> list:
        """
        Returns the shortest path between two nodes of the molecule based on the Dijkstra algorithm from networkx.
        Args:
            source (int): source node.
            target (int): target node.
            include_all (bool, optional): If all nodes should be considered in the path. Defaults to False.
        Returns:
            list: Shortest path between source and target.
        """
        if include_all:
            return nx.shortest_path(self, source, target)
        subgraph = nx.Graph()
        visible_nodes = self.return_visible_nodes()
        subgraph.add_nodes_from(visible_nodes)
        subgraph.add_edges_from(self.return_visible_edges())
        return nx.shortest_path(subgraph, source, target)

    def _add_terminal_attribute(self) -> None:
        """
        Adds the terminal attribute to the edges of the graph.
        """
        # Get the degree of each visible node in the order of visible_nodes.
        degrees = self.get_degree()
        visible_nodes = self.return_visible_nodes()
        for edge in self.edges:
            degree1 = degrees[visible_nodes.index(edge[0])]
            degree2 = degrees[visible_nodes.index(edge[1])]
            # If one of the nodes has degree 1, the edge is terminal.
            if degree1 == 1 or degree2 == 1:
                self.add_edge(edge[0], edge[1], terminal=True, rigid=True)

    def _add_stick_attribute(self) -> None:
        """
        Adds the stick attribute to the edges of the graph.
        """
        # Get the degree of each visible node in the order of visible_nodes.
        degrees = self.get_degree()
        elements = self.return_node_data(attribute="elements")
        visible_nodes = self.return_visible_nodes()
        nodes_with_condition = [
            node
            for idx, node in enumerate(visible_nodes)
            if (degrees[idx] == 2) and (elements[idx] == "C")
        ]
        self.set_nodes_visible(visible_nodes=nodes_with_condition)
        connected_components = self.connected_components()
        for component in connected_components:
            if len(component) >= 2:
                self.set_nodes_visible(visible_nodes=component)
                self.add_node_data(attribute="stick", new_data=[True] * len(component))
                for edge in self.return_visible_edges(nodes_visible=True):
                    self.add_edge(edge[0], edge[1], stick=True)
        self.set_nodes_visible(visible_nodes=visible_nodes)

    def _add_bond_constraints(self) -> None:
        """
        Adds bond constraints to the graph.
        """
        if self.constraints.get("constraints") is not None:
            extra_constraints = self.constraints["constraints"].get("bonds")
            if extra_constraints is not None:
                for extra_constraint in extra_constraints:
                    nodes = extra_constraint["atoms"]
                    self.add_edge(
                        nodes[0] - 1,
                        nodes[1] - 1,
                        bo=1,
                        visible=True,
                        dij=extra_constraint["value"],
                    )
