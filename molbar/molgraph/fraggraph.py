import numpy as np
import networkx as nx
from molbar.helper.rank import Rank
from molbar.molgraph.molgraph import MolGraph


class FragGraph:

    def __init__(self, molgraph: MolGraph) -> None:
        """
        Constructs a FragGraph object.

        Args:
            molgraph (MolGraph): A MolGraph object.
        """

        self.molgraph = molgraph

        self.fragments = self.get_fragments()

        self.n_fragments = self.return_n_fragments()

    def get_fragments(self) -> list:
        """
        Returns the fragments of the molecule.

        Returns:
            list: lists of nodes of the fragments.
        """

        visible_nodes = self.molgraph.return_visible_nodes()

        fragment_id_per_node = self.molgraph.return_node_data(attribute="fragment_id")

        self.fragment_ids = list(self.molgraph.fragments_data.keys())

        self.fragments = [
            [
                node
                for node in visible_nodes
                if fragment_id_per_node[node] == fragment_id
            ]
            for fragment_id in self.fragment_ids
        ]

        return self.fragments

    def return_n_fragments(self):
        """
        Returns the number of fragments.

        Returns:
            int: number of fragments
        """

        return len(self.fragments)

    def get_abs_conf_indices(self) -> np.ndarray:
        """
        Returns the absolute configuration indices of the fragments in the order of self.fragment_ids.

        Returns:
            np.ndarray: Absolute configuration indices of the fragments.
        """

        abs_conf_indices = np.array(
            [
                self.molgraph.fragments_data[fragment_id]["absconf_index"]
                for fragment_id in self.fragment_ids
            ]
        )

        return abs_conf_indices

    def return_fragment_priorities(self) -> np.ndarray:
        """
        Returns the priorities of the fragments in the order of self.fragment_ids.

        Returns:
            np.ndarray: Priorities of the fragments.
        """

        visible_nodes = self.molgraph.return_visible_nodes()

        priorities = np.round(self.molgraph.return_node_data(attribute="priorities"))

        fragments_priorities = {
            id: sorted(
                [priorities[visible_nodes.index(node)] for node in fragment],
                reverse=True,
            )
            for id, fragment in zip(self.fragment_ids, self.fragments)
        }

        merged_fragment_priorities = self._get_fragment_priorities(fragments_priorities)

        return np.array(
            [
                merged_fragment_priorities[fragment_id]
                for fragment_id in self.fragment_ids
            ]
        )

    def _get_fragment_priorities(self, priorities) -> list:
        """
        Determines priorities of fragments.

        Returns:
            list: Priorities of fragments.
        """

        ids = list(priorities.keys())

        priorities = list(priorities.values())

        # Find the maximum number of atoms in any fragment
        max_number_of_atoms_in_fragment = 0

        for j in range(len(priorities)):
            max_number_of_atoms_in_fragment = max(
                max_number_of_atoms_in_fragment, len(priorities[j])
            )

        # Create a list of filler values to pad the fragment priorities to the same length
        filler = [0] * max_number_of_atoms_in_fragment

        # Pad the priorities of atoms in each fragment and calculate the fragment priorities
        equalized_priorities_of_atoms_in_fragments = [
            sublist[:max_number_of_atoms_in_fragment] + filler[len(sublist) :]
            for sublist in priorities
        ]

        return {
            id: prio
            for id, prio in zip(
                ids, Rank(equalized_priorities_of_atoms_in_fragments).priorities
            )
        }

    def get_fragment_distance_matrix(self) -> np.ndarray:
        """
        Returns the fragment graph distance matrix.

        Returns:
            np.ndarray: Distance matrix of the molecule in the order of self.fragment_ids.
        """

        adjacent_nodes = [
            self.molgraph.return_adjacent_nodes(core_nodes=fragment)
            for fragment in self.fragments
        ]

        fragment_graph = self._set_up_fragment_graph(adjacent_nodes)

        return self._get_fragment_distance_matrix(fragment_graph)

    def _get_fragment_distance_matrix(self, fragment_graph: nx.Graph) -> np.ndarray:

        distance_matrix = np.zeros((self.n_fragments, self.n_fragments))

        for ith_fragment_id in self.fragment_ids:

            for jth_fragment_id in self.fragment_ids:

                if jth_fragment_id > ith_fragment_id:

                    try:

                        distance = nx.shortest_path_length(
                            fragment_graph,
                            source=ith_fragment_id,
                            target=jth_fragment_id,
                        )

                    except nx.exception.NetworkXNoPath:

                        distance = np.inf

                    distance_matrix[ith_fragment_id, jth_fragment_id] = distance

                    distance_matrix[jth_fragment_id, ith_fragment_id] = distance

        return distance_matrix

    def _set_up_fragment_graph(self, adjacent_nodes: list) -> nx.Graph:
        """
        Sets up the fragment graph.

        Args:
            adjacent_nodes (list): adjacent nodes of the fragments in the order of self.fragment_ids.

        Returns:
            nx.Graph: nx.Graph object of the fragment graph.
        """

        fragment_graph = nx.Graph()

        fragment_graph.add_nodes_from(self.fragment_ids)

        for ith_id in range(self.n_fragments):

            ith_fragment = self.fragments[ith_id]

            for jth_id in range(self.n_fragments):

                if jth_id > ith_id:

                    if (
                        len(set(ith_fragment).intersection(set(adjacent_nodes[jth_id])))
                        > 0
                    ):

                        fragment_graph.add_edge(ith_id, jth_id)

        return fragment_graph
