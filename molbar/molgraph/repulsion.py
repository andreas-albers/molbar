import numpy as np
from typing import Union


class Repulsion:

    def return_repulsion_pairs(self) -> Union[list, list]:
        """
        Returns the repulsion pairs and the corresponding charges for visible node pairs.

        Returns:
            list: List of repulsion pairs.
            list: List of charges for the repulsion pairs.
        """
        visible_nodes = self.return_visible_nodes()

        repulsion = []

        charges = []

        for i in visible_nodes:

            for j in visible_nodes:

                if j > i:

                    if (
                        j not in self.visible_repulsion_nodes
                        or i not in self.visible_repulsion_nodes
                    ):

                        continue

                    ith_new_index = self.visible_repulsion_nodes.index(i)

                    jth_new_index = self.visible_repulsion_nodes.index(j)

                    charge = self.repulsion_pairs[ith_new_index, jth_new_index]

                    if charge == 0.0:

                        continue

                    repulsion.append([i, j])

                    charges.append(charge)

        return repulsion, charges

    def _get_repulsion_pairs(self):
        """
        Generates the repulsion pairs for the molecule.
        Sets the repulsion pairs and the visible repulsion nodes as attributes of the MolGraph object.
        """

        visible_nodes = self.return_visible_nodes()

        if self.constraints.get("repulsion_charge") != None:

            self.charge = self.constraints.get("repulsion_charge")

        else:

            self.charge = 1e2

        if self.constraints.get("unique_repulsion") == True:

            repulsion = self._get_unique_repulsion_pairs(visible_nodes)

        else:

            repulsion = []

        repulsion_pairs = (
            np.ones((len(visible_nodes), len(visible_nodes))) * self.charge
        )

        for combination in repulsion:

            repulsion_pairs[combination[0], combination[1]] = 5 * self.charge

            repulsion_pairs[combination[1], combination[0]] = 5 * self.charge

        self.repulsion_pairs = repulsion_pairs

        self.visible_repulsion_nodes = visible_nodes

    def _get_unique_repulsion_pairs(self, visible_nodes: list) -> list:
        """
        Generates the unique repulsion pairs for the molecule.
        Searches for single bonds that are not part of a cycle and adds the corresponding nodes to the repulsion pairs.
        Gets the indices of the adjacent nodes to that bond with the highest priority and adds them to the repulsion pairs with higher charge.

        Args:
            visible_nodes (list): List of visible nodes of the molecule.

        Returns:
            list: List of unique repulsion pairs.
        """

        repulsion = []

        visible_edges = self.return_visible_edges(nodes_visible=True)

        for edge in visible_edges:

            if any(
                [
                    True if edge[0] in cycle and edge[1] in cycle else False
                    for cycle in self.cycle_nodes
                ]
            ):

                continue

            if self.edges[edge]["bo"] == 1.0:

                adjacent_nodes_0 = self.return_adjacent_nodes(core_nodes=edge[0])

                adjacent_nodes_0.remove(edge[1])

                adjacent_nodes_1 = self.return_adjacent_nodes(core_nodes=edge[1])

                adjacent_nodes_1.remove(edge[0])

                if len(adjacent_nodes_0) <= 1 or len(adjacent_nodes_1) <= 1:

                    continue

                single_bond_nodes = set(
                    list(edge) + adjacent_nodes_0 + adjacent_nodes_1
                )

                self.set_nodes_visible(visible_nodes=single_bond_nodes)

                single_bond_nodes = self.return_visible_nodes()

                adjacent_indices_0 = [
                    single_bond_nodes.index(node) for node in adjacent_nodes_0
                ]

                adjacent_indices_1 = [
                    single_bond_nodes.index(node) for node in adjacent_nodes_1
                ]

                priorities = self.return_node_data(attribute="priorities")

                adjacent_priorities_0 = [
                    priorities[node] for node in adjacent_indices_0
                ]

                adjacent_priorities_1 = [
                    priorities[node] for node in adjacent_indices_1
                ]

                first_index = first_unique_or_max(adjacent_priorities_0)

                second_index = first_unique_or_max(adjacent_priorities_1)

                repulsion.append(
                    sorted(
                        [
                            single_bond_nodes[first_index],
                            single_bond_nodes[second_index],
                        ]
                    )
                )

                self.set_nodes_visible(visible_nodes=visible_nodes)

        return repulsion


def first_unique_or_max(lst) -> int:
    """
    Returns the index of the first unique value in a list or the index of the maximum value if no unique value is found.

    Args:
        lst (list): list object

    Returns:
        int: index of unique element of list
    """
    # Create a dictionary to store the frequency of each element
    frequency = {}

    # Count the frequency of each element in the list
    for item in lst:
        if item in frequency:
            frequency[item] += 1
        else:
            frequency[item] = 1

    # Find the first unique value or the maximum value
    result = None

    for item in sorted(lst, reverse=True):
        if frequency[item] == 1:
            result = item
            break  # Break the loop when the first unique value is found

    # If no unique value is found, return the maximum value
    if result is None:
        result = sorted(lst, reverse=True)[0]

    return lst.index(result)
