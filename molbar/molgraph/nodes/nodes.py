import os
import json
import numpy as np
import networkx as nx
from collections.abc import Iterable
from scipy.spatial.distance import cdist

with open(
    os.path.join(
        os.path.dirname(__file__), "..", "..", "data", "element_data.json"
    )
) as json_file:
    node_data = json.load(json_file)

class Nodes:

    def return_nodes(self, attributes=None, values=None, include_all=False) -> list:
        """
        Returns the nodes of the graph that have the specified attributes and values.

        Args:
            attributes (Union(list, str), optional): Attributes to filter the nodes. Defaults to None.
            values (_type_, optional): Values the attributes of the nodes should have to filter the nodes. Defaults to None.
            include_all (bool, optional): If True, all nodes are returned. Defaults to False.

        Raises:
            ValueError: If attributes and values are provided as lists, they must have the same length.

        Returns:
            list: List of nodes that have the specified attributes and values.
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        elif (isinstance(attributes, Iterable) and isinstance(values, Iterable)) and (
            not isinstance(attributes, str) and not isinstance(values, str)
        ):
            raise ValueError("attributes and values must have the same length")
        elif not include_all:
            raise ValueError(
                "You must provide either (attributes and values) as str (one attribute) or as lists/np.array (several attributes) or visible_nodes as a list."
            )
        if include_all:
            return [
                node
                for node, data in self.nodes(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
            ]
        else:
            return [
                node
                for node, data in self.nodes(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
                and (data.get("visible"))
            ]

    def return_node_data(self, attribute: str, include_all=False) -> np.array:
        """
        Returns the data of the nodes of the graph that have the specified attributes and values.

        Args:
            attribute (str): Attribute to filter the nodes.
            include_all (bool, optional): If True, all nodes are returned. Defaults to False.

        Returns:
            np.array: Array of data of the nodes that have the specified attributes and values.
        """
        if include_all:
            data = [data.get(attribute) for node, data in self.nodes(data=True)]
        else:
            data = [
                data.get(attribute)
                for node, data in self.nodes(data=True)
                if data.get("visible")
            ]
        try:
            return np.array(data)
        except ValueError:
            return data

    def add_node_data(self, attribute: str, new_data: list) -> None:
        """
        Adds new data to the nodes of the graph that have the specified attributes and values.

        Args:
            attribute (str): New attribute to add to the nodes.
            new_data (list): Data to add to the nodes.

        Raises:
            AssertionError: If the length of the new data is not the same as the number of visible nodes.

        Returns:
            None
        """
        visible_nodes = [
            (node, data) for node, data in self.nodes(data=True) if data.get("visible")
        ]
        assert len(visible_nodes) == len(
            new_data
        ), "new data must have the same length as the number of visible nodes"
        for node_data, new_data_point in zip(visible_nodes, new_data):
            if node_data[1].get("visible"):
                node_data[1][attribute] = new_data_point
                self.add_node(node_data[0], **node_data[1])

    def _get_data_for_element(self, element: str):
        """
        Returns the data for the specified element.

        Args:
            element (str): Element to get the data for.

        Returns:
            dict: Dictionary with the data for the specified element. {"atomic_numbers": atomic_number, "rcov":rcov, "en_pauling":en_pauling, "cn_fak":cn_fak, "per_fak":per_fak, "masses":masses, "is_metal": is_metal, "atomic_valences": valences}
        """

        atomic_number = self.element_data["atomic_numbers"].get(element)
        rcov = self.element_data["rcov"][atomic_number - 1]
        en_pauling = self.element_data["en_pauling2"][atomic_number - 1]
        cn_fak = self.element_data["cn_fak"][atomic_number - 1]
        per_fak = self.element_data["per_fak"][atomic_number - 1]
        masses = self.element_data["masses"][atomic_number - 1]
        is_metal = self.element_data["is_metal"][atomic_number - 1]
        valences = self.element_data["atomic_valences"][atomic_number - 1]
        normcn = self.element_data["normcn"][atomic_number - 1]
        rnorm = self.element_data["rnorm"][atomic_number - 1]
        valence_electrons = self.element_data["valence_electrons"][atomic_number - 1]
        return {
            "valence_electrons": valence_electrons,
            "atomic_numbers": atomic_number,
            "rcov": rcov,
            "en_pauling": en_pauling,
            "cn_fak": cn_fak,
            "per_fak": per_fak,
            "masses": masses,
            "is_metal": is_metal,
            "atomic_valences": valences,
            "normcn": normcn,
            "rnorm": rnorm,
        }

    def return_distance_matrix(self, include_all=False):
        """
        Returns the distance matrix of the graph.

        Args:
            include_all (bool, optional): If True, all nodes are returned. Defaults to False.

        Returns:
            np.array: Cartesian distance matrix of the nodes in the graph.
        """

        if include_all:
            return cdist(
                self.return_node_data(attribute="coordinates", include_all=True),
                self.return_node_data(attribute="coordinates", include_all=True),
            )
        return cdist(
            self.return_node_data(attribute="coordinates"),
            self.return_node_data(attribute="coordinates"),
        )

    def set_nodes_visible(
        self, attributes=None, values=None, visible_nodes=None, include_all=False
    ):
        """
        Sets the nodes of the graph that have the specified attributes and values to visible.

        Args:
            attributes (_type_, optional): _description_. Defaults to None.
            values (_type_, optional): _description_. Defaults to None.
            visible_nodes (_type_, optional): _description_. Defaults to None.
            include_all (bool, optional): _description_. Defaults to False.

        Raises:
            ValueError: _description_
            ValueError: _description_
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        elif (isinstance(attributes, Iterable) and isinstance(values, Iterable)) and (
            not isinstance(attributes, str) and not isinstance(values, str)
        ):
            assert len(attributes) == len(
                values
            ), "attributes and values must have the same length"
        elif isinstance(visible_nodes, Iterable) and not isinstance(visible_nodes, str):
            pass
        elif isinstance(visible_nodes, int):
            visible_nodes = [visible_nodes]
        elif not include_all:
            raise ValueError(
                "You must provide either (attributes and values) as str (one attribute) or as lists/np.array (several attributes) or visible_nodes as a list."
            )
        if include_all:
            new_attributes = {node: {"visible": True} for node in self.nodes()}
            nx.set_node_attributes(self, new_attributes)
            return
        new_attributes = {node: {"visible": False} for node in self.nodes()}
        nx.set_node_attributes(self, new_attributes)
        if visible_nodes is not None:
            new_attributes = {node: {"visible": True} for node in visible_nodes}
            nx.set_node_attributes(self, new_attributes)
        elif attributes and values:
            visible_nodes = [
                node
                for node, data in self.nodes(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
            ]
            # New attributes to set
            new_attributes = {node: {"visible": True} for node in visible_nodes}
            # Set new attributes
            nx.set_node_attributes(self, new_attributes)
        else:
            raise ValueError(
                "Either (attribute and value) or visible_nodes must be provided."
            )

    def set_nodes_unvisible(
        self, attributes=None, values=None, unvisible_nodes=None, include_all=False
    ):
        """
        Sets the nodes of the graph that have the specified attributes and values to unvisible.

        Args:
            attribute (str): Attribute to filter the nodes.
            value (_type_, optional): Value the attribute of the nodes should have to filter the nodes. Defaults to None.
            include_all (bool, optional): If True, all nodes are set unvisble. Defaults to False.
        """
        if isinstance(attributes, str):
            attributes = [attributes]
            values = [values]
        elif (isinstance(attributes, Iterable) and isinstance(values, Iterable)) and (
            not isinstance(attributes, str) and not isinstance(values, str)
        ):
            assert len(attributes) == len(
                values
            ), "attributes and values must have the same length"
        elif isinstance(unvisible_nodes, Iterable) and not isinstance(
            unvisible_nodes, str
        ):
            pass
        elif isinstance(unvisible_nodes, int):
            unvisible_nodes = [unvisible_nodes]

        elif not include_all:
            raise ValueError(
                "You must provide either (attributes and values) as str (one attribute) or as lists/np.array (several attributes) or visible_nodes as a list."
            )
        if include_all:
            new_attributes = {node: {"visible": False} for node in self.nodes()}
            nx.set_node_attributes(self, new_attributes)
            return
        if unvisible_nodes is not None:
            new_attributes = {node: {"visible": False} for node in unvisible_nodes}
            nx.set_node_attributes(self, new_attributes)
        elif attributes and values:
            unvisible_nodes = [
                node
                for node, data in self.nodes(data=True)
                if all(
                    [
                        data.get(attribute) == value
                        for attribute, value in zip(attributes, values)
                    ]
                )
            ]
            new_attributes = {node: {"visible": False} for node in unvisible_nodes}
            nx.set_node_attributes(self, new_attributes)
        else:
            raise ValueError(
                "Either (attribute and value) or visible_nodes must be provided."
            )

    def return_adjacent_nodes(
        self, core_nodes=None, include_all=False, is_adjacent_visible=True
    ):
        """
        Returns the adjacent nodes of the specified nodes.

        Args:
            core_nodes (Union(list, np.int64)): Nodes to get the adjacent nodes for.
            include_all (bool, optional): If True, all nodes are considered as core nodes. Defaults to False.
            is_adjacent_visible (bool, optional): If True, only visible nodes are considered to be adjacent to core nodes. Defaults to True.

        Returns:
            list: List of adjacent nodes of the specified nodes.
        """
        if isinstance(core_nodes, int) or isinstance(core_nodes, np.int64):
            core_nodes = [core_nodes]
        if include_all:
            adjacent_nodes = [self.neighbors(node) for node in self.nodes()]
        else:
            nodes_data = self.nodes(data=True)
            if is_adjacent_visible:
                adjacent_nodes = [
                    adjacent_node
                    for core_node in core_nodes
                    if nodes_data[core_node].get("visible")
                    for adjacent_node in self.neighbors(core_node)
                    if (nodes_data[adjacent_node].get("visible"))
                    and (adjacent_node not in core_nodes)
                ]
            else:
                adjacent_nodes = [
                    adjacent_node
                    for core_node in core_nodes
                    if nodes_data[core_node].get("visible")
                    for adjacent_node in self.neighbors(core_node)
                    if (adjacent_node not in core_nodes)
                ]
        return list(set(adjacent_nodes))

    def return_visible_nodes(self, return_data=False, include_all=False) -> list:
        """
        Returns the visible nodes of the graph.

        Args:
            data (bool, optional): If True, the data of the visible nodes is returned. Defaults to False.
            include_all (bool, optional): If True, all nodes are considered as visible. Defaults to False.

        Returns:
            list: List of visible nodes of the graph.
        """
        if include_all:
            return self.nodes(data=return_data)
        else:
            if return_data:
                return [
                    (node, data)
                    for node, data in self.nodes(data=True)
                    if data.get("visible")
                ]
            else:
                return [
                    node for node, data in self.nodes(data=True) if data.get("visible")
                ]
