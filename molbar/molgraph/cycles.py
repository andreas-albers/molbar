import networkx as nx
import numpy as np
from collections import defaultdict


class Cycle:

    def detect_cycles(self, include_all=False) -> None:
        """
        Detects cycles in the molecule and adds them to the nodes and edges.

        Args:
            include_all (bool, optional): If all nodes and not only the visible no  des should be considered. The default value is False.

        Returns:
            None
        """

        first_visible_nodes = self.return_visible_nodes()

        if include_all:

            cycle_basis = nx.minimum_cycle_basis(self)

        else:

            self.set_nodes_visible(attributes="is_metal", values=False)

            visible_nodes = self.return_visible_nodes()

            cycle_basis = self._get_minimum_cycle_basis(self.subgraph(visible_nodes))

        cycles = [cycle for cycle in cycle_basis if len(cycle) <= 8]

        self._find_connected_cycles(cycles)

        if not include_all:

            cycles = self._add_metal_cycles(cycles)

        self._add_cycles_to_nodes(cycles)

        self._add_cycles_to_edges()

        self.set_nodes_visible(visible_nodes=first_visible_nodes)

        self.add_node_data(
            attribute="smallest_cycle_size", new_data=[0] * len(first_visible_nodes)
        )

        self._add_smallest_cycle_size_to_nodes(cycles)

        self.set_nodes_visible(include_all=True)

        self.cycle_nodes = cycles

    def _find_connected_cycles(self, cycles: list) -> list:
        adjacent_nodes_per_cycle = []
        cycles_nodes = set(node for cycle in cycles for node in cycle)
        self.set_nodes_visible(visible_nodes=cycles_nodes)
        for cycle in cycles:
            adjacent_nodes = self.return_adjacent_nodes(cycle)
            adjacent_nodes_per_cycle.append(adjacent_nodes)
        connected_cycles = []
        for i, cycle in enumerate(cycles):
            for j, adjacent_nodes in enumerate(adjacent_nodes_per_cycle):
                if j > i:
                    if len(set(cycle).intersection(set(adjacent_nodes))) > 1:
                        connected_cycles.append([i, j])
        for connected_cycle in connected_cycles:
            first_cycle = set(cycles[connected_cycle[0]])
            second_cycle = set(cycles[connected_cycle[1]])
            both_cycles = first_cycle.union(second_cycle)
            self.set_nodes_visible(visible_nodes=both_cycles)
            edges = self.return_visible_edges(nodes_visible=True)
            for edge in edges:
                if (
                    (edge[0] in first_cycle)
                    and (edge[1] in second_cycle)
                    or (edge[0] in second_cycle)
                    and (edge[1] in first_cycle)
                ):
                    self.add_edge(edge[0], edge[1], rigid=True)

    def _add_cycles_to_nodes(self, cycles: list) -> None:
        """
        Adds the cycles to the nodes data structure as cycle_ids.

        Args:
            cycles (list): list of cycles

        Returns:
            None
        """

        for i, cycle in enumerate(cycles):

            for node in cycle:

                if "cycle_id" in self.nodes[node]:

                    self.nodes[node]["cycle_id"].append(i)

                else:

                    self.add_node(node, cycle_id=[i])

    def _add_cycles_to_edges(self) -> None:
        """
        Adds the cycles to the edges data structure as cycle_ids.

        Returns:
            None
        """

        for edge in self.edges():

            self.set_nodes_visible(visible_nodes=edge)

            cycle_ids = self.return_node_data("cycle_id")

            if (cycle_ids[0] is not None) and (cycle_ids[1] is not None):

                intersection = set(cycle_ids[0]).intersection(set(cycle_ids[1]))

                if intersection:

                    self.add_edge(edge[0], edge[1], cycle_id=list(intersection))

    def _get_minimum_cycle_basis(self, subgraph: nx.Graph) -> list:

        cycle_basis = nx.minimum_cycle_basis(subgraph)

        cycles = self._refine_cycle_basis(cycle_basis)

        return cycles

    def _refine_cycle_basis(self, cycles: list) -> list:

        cycle_nodes = set(node for cycle in cycles for node in cycle)

        self.minimal_cycle_sizes = self._get_minimal_cycle_sizes(cycles)

        self.set_nodes_visible(visible_nodes=cycle_nodes)

        cycles = [
            cycle
            for starting_node in cycle_nodes
            for cycle in self._find_cycles_from_starting_node(starting_node)
        ]

        cycles = [
            list(unique_cycle)
            for unique_cycle in set(tuple(sorted(cycle)) for cycle in cycles)
        ]

        # cycles = self._remove_epoxide_rings(cycles)

        minimal_cycle_sizes = [self.minimal_cycle_sizes[node] for node in cycle_nodes]

        self.add_node_data(
            attribute="smallest_cycle_size", new_data=minimal_cycle_sizes
        )

        return cycles

    def _find_cycles_from_starting_node(self, starting_node: int) -> list:

        smaller_cycles = set()

        start_adjacent_nodes = self.return_adjacent_nodes(
            [starting_node], is_adjacent_visible=True
        )

        queue = set(
            tuple([starting_node, adjacent_atom])
            for adjacent_atom in start_adjacent_nodes
        )

        while queue:

            current_cycle = list(queue.pop())

            latest_node = current_cycle[-1]

            for adjacent_node in self.return_adjacent_nodes(
                [latest_node], is_adjacent_visible=True
            ):

                if (adjacent_node == starting_node) & (len(current_cycle) > 2):

                    for new_cycle_node in current_cycle:

                        if self.minimal_cycle_sizes[new_cycle_node] > len(
                            current_cycle
                        ):

                            self.minimal_cycle_sizes[new_cycle_node] = len(
                                current_cycle
                            )

                    if tuple(sorted(current_cycle)) not in [
                        tuple(sorted(ring)) for ring in smaller_cycles
                    ]:

                        smaller_cycles.add(tuple(current_cycle))

                # If the adjacent atom is already part of the current ring, skip it
                elif adjacent_node in current_cycle:

                    continue

                # If the ring size exceeds 12 atoms, skip it. And if the ring size is larger than the smallest ring size of the startung atom, skip it. (we want the smallest ring)
                elif (len(current_cycle) > self.minimal_cycle_sizes[starting_node]) | (
                    len(current_cycle) >= 8
                ):

                    continue

                else:
                    # Expand the current ring with the adjacent atom and add it to the queue
                    expanded_ring = current_cycle + [adjacent_node]

                    if tuple(sorted(expanded_ring)) not in [
                        tuple(sorted(ring)) for ring in smaller_cycles
                    ]:

                        queue.add(tuple(expanded_ring))

        return [
            cycle
            for cycle in smaller_cycles
            if (len(cycle) == self.minimal_cycle_sizes[starting_node])
        ]

    def _get_minimal_cycle_sizes(self, cycles):

        # Create a dictionary with the smallest ring size for each ring atom
        minimum_ring_size_per_atom = defaultdict(lambda: float("inf"))

        # Loop over all rings in the cycle basis
        for cycle in cycles:
            # Loop over all atoms in the ring
            for atom in cycle:
                # Update the smallest ring size for each atom in the ring, if the current ring is smaller than the stored smallest ring size
                minimum_ring_size_per_atom[atom] = min(
                    len(cycle), minimum_ring_size_per_atom[atom]
                )
        # Return the dictionary
        return dict(minimum_ring_size_per_atom)

    def _add_metal_cycles(self, cycles: list) -> None:
        """

        Detects additonal metal cycles.
        Metal cycles are cycles in the graph that contains one metal node and one node that has not been part of a cycle before.

        Args:
            cycles (list): list of cycles

        Returns:
            cycles(list): updated list of cycles
        """

        self.set_nodes_visible(include_all=True)

        visible_nodes = self.return_visible_nodes()

        metal_nodes = self.return_nodes(attributes="is_metal", values=True)

        if len(metal_nodes) != 0:

            visible_nodes += metal_nodes

            self.set_nodes_visible(visible_nodes=visible_nodes)

            metal_cycle_basis = self._get_minimum_cycle_basis(
                self.subgraph(visible_nodes)
            )

            metal_cycles = [
                cycle
                for cycle in metal_cycle_basis
                if (len(cycle) <= 8)
                and (len(cycle) > 3)
                and any([metal_node in cycle for metal_node in metal_nodes])
            ]

            for cycle in metal_cycles:

                visible_nodes = [node for node in cycle if node not in metal_nodes]

                self.set_nodes_visible(visible_nodes=visible_nodes)

                if any(
                    data is None for data in self.return_node_data(attribute="cycle_id")
                ):

                    cycles.append(cycle)

        return cycles

    def _add_smallest_cycle_size_to_nodes(self, cycles: list) -> None:
        """
        Adds the size of the smallest cycle a node is part of to the nodes data structure as smallest_cycle_size attribute.

        Args:
            cycles (list): list cycles

        Returns:
            None
        """

        cycle_nodes = set(node for cycle in cycles for node in cycle)

        minimal_cycle_sizes = {}

        for node in cycle_nodes:

            min_size = np.inf

            for cycle in cycles:

                if (node in cycle) and (len(cycle) < min_size):

                    min_size = len(cycle)

            minimal_cycle_sizes[node] = min_size

        self.set_nodes_visible(visible_nodes=cycle_nodes)

        visible_nodes = self.return_visible_nodes()

        minimal_cycle_sizes = [minimal_cycle_sizes[node] for node in visible_nodes]

        self.add_node_data(
            attribute="smallest_cycle_size", new_data=minimal_cycle_sizes
        )

    def _remove_epoxide_rings(self, cycles: list) -> list:

        filtered_cycles = []

        for cycle in cycles:

            if len(cycle) == 3:

                cycle_nodes = set(cycle)

                degree = self.degree()

                if any(degree[node] == 2 for node in cycle_nodes):

                    continue

            filtered_cycles.append(cycle)

        return filtered_cycles
