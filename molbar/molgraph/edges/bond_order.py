import numpy as np
import itertools
import molbar.fortranlib.bo_matrix as bo


class BondOrder:

    def __init__(self) -> None:
        pass

    def assign_bond_orders(self):
        """
        Adds bond orders as edge weights to the graph.

        Returns:
            None
        """
        if self.custom_bond_orders is not None:
            if 'double' in self.custom_bond_orders:
                for edge in self.custom_bond_orders['double']:
                    self.add_edge(
                        edge[0],
                        edge[1],
                        bo=2,
                        rigid=True,
                    )
            if 'tripple' in self.custom_bond_orders:
                for edge in self.custom_bond_orders['tripple']:
                    self.add_edge(
                        edge[0],
                        edge[1],
                        bo=3,
                        rigid=True,
                    )

        else:
            initial_visible_nodes = self.return_visible_nodes()

            # Set metal nodes invisible as they are not considered for bond order detection (set to bond order of 1 automically)
            self.set_nodes_visible(attributes="is_metal", values=False)

            elements = self.return_node_data(attribute="elements")

            visible_nodes = self.return_visible_nodes()

            # Get degree (e.g. coordination number) of all visible nodes
            degree = self.get_degree()

            # Get valences of all visible nodes as tabulated for each element. 4 for carbon, 2 for oxygen, etc.
            valences = self.return_node_data(attribute="atomic_valences")

            # Check if atom is possible unsaturated based on coordination number and valence and set attribute "possible_unsaturated" to True if so. If True, atom might be involved in multiple bonds.
            is_unsaturated = np.array(
                [
                    (
                        True
                        if (node_data[1].get("visible") == True)
                        and (degree[i] < valences[i][0])
                        and (
                            node_data[1].get("stick") is None
                            or (
                                node_data[1].get("stick") == True
                                and node_data[1].get("cycle_id") is not None
                            )
                        )
                        else False
                    )
                    for i, node_data in enumerate(
                        self.return_visible_nodes(return_data=True)
                    )
                ]
            )

            is_unsaturated = self._saturate_groups(
                is_unsaturated, visible_nodes, elements, degree, valences
            )

            # Add attribute "possible_unsaturated" to all nodes
            self.add_node_data(attribute="possible_unsaturated", new_data=is_unsaturated)

            # Set all nodes invisible that are not possible unsaturated
            self.set_nodes_visible(
                attributes=["possible_unsaturated", "is_metal"], values=[True, False]
            )

            # Get fragments of connected possible unsaturated atoms
            unsaturated_fragments = self.connected_components()

            # Iterate over all fragments
            for core_nodes in unsaturated_fragments:
                # Only bond order detection for fragments with more than one node, otherwise there is no other possibility that atom is unsaturated.
                if len(core_nodes) > 1:

                    # Determine bond orders for one unsaturated fragment
                    self._get_bond_order_for_fragment(core_nodes, visible_nodes)

            self.set_nodes_visible(visible_nodes=initial_visible_nodes)

        self._add_valences_to_nodes()

        self._flexify_nitroso()

        self._flexify_nitrogen_valence_5()

    def _get_bond_order_for_fragment(
        self, core_nodes: list, visible_nodes: list
    ) -> None:
        """
        Assigns bond orders to a one unsaturated fragment.

        Args:
            core_nodes (list): core nodes of one fragment
            visible_nodes (list): visible nodes of the molecule
        """
        self.set_nodes_visible(visible_nodes=visible_nodes)
        # Get adjacent nodes of the unsaturated fragment
        adjacent_nodes = self.return_adjacent_nodes(core_nodes)
        # print(np.array(core_nodes)+1)
        # Only consider the fragment nodes for further bond order detection
        self.set_nodes_visible(visible_nodes=core_nodes + adjacent_nodes)

        n_atoms = self.return_n_atoms()

        visible_nodes = self.return_visible_nodes()
        # print(np.array(core_nodes)+1, np.array(visible_nodes)+1)
        # Determine which atoms are core atoms, so possible unsaturated. Adjacent atoms are saturated by definition.
        are_core = [1 if node in core_nodes else 0 for node in visible_nodes]

        # Get valences of all atoms as tabulated for each element. 4 for carbon, 2 for oxygen, etc in the order of visible_nodes.
        atomic_valences = self.return_node_data(attribute="atomic_valences")

        # Get degree (e.g. coordination number) of all atoms in the order of visible_nodes
        degrees = self.get_degree()

        # If atom is core atom, then all valences are possible. If atom is adjacent atom, then only the coordination number is possible.
        atomic_valences = [
            list(valence) if is_core == 1 else [degree]
            for valence, degree, is_core in zip(atomic_valences, degrees, are_core)
        ]

        # Get atomic numbers of all atoms in the order of visible_nodes (e.g. 6 for carbon, 8 for oxygen, etc.)
        atomic_numbers = self.return_node_data(attribute="atomic_numbers")

        # Get elements of all atoms in the order of visible_nodes (e.g. "C" for carbon, "O" for oxygen, etc.)
        elements = self.return_node_data(attribute="elements")

        # Get all possible valence combinations for all atoms in the order of visible_nodes
        all_valence_combinations = self._get_combinations_of_valences(
            elements, atomic_valences, degrees, are_core, visible_nodes
        )

        # Get connectivity matrix of the fragment
        cn_matrix = self.return_cn_matrix()
        # Get best bond order matrix for the fragment
        fragment_bo_matrix = bo.bond_order_detection.get_best_bo_matrix(
            n_atoms,
            cn_matrix,
            degrees,
            atomic_numbers,
            all_valence_combinations,
            are_core,
        )
        # print(fragment_bo_matrix)
        self._add_weights(fragment_bo_matrix, elements, degrees, visible_nodes)

        self.set_nodes_visible(visible_nodes=visible_nodes)

        # self._add_amide_group()

    def _add_valences_to_nodes(self):

        visible_nodes = self.return_visible_nodes()

        visible_edges = self.return_visible_edges(nodes_visible=True)

        valences = np.zeros(len(visible_nodes), dtype=int)

        for edge in visible_edges:

            bo = self.edges[edge]["bo"]

            valences[visible_nodes.index(edge[0])] += bo

            valences[visible_nodes.index(edge[1])] += bo

        self.add_node_data(attribute="valences", new_data=valences)

        valences = self.return_node_data(attribute="valences")

    def _add_weights(self, bo_matrix, elements, degrees, visible_nodes):

        for i in range(bo_matrix.shape[0]):

            for j in range(bo_matrix.shape[1]):

                if (j > i) and (bo_matrix[i, j] > 1):

                    ith_element = elements[i]

                    jth_element = elements[j]

                    if (ith_element == "C") & (jth_element == "C"):

                        if degrees[i] == 3 and degrees[j] == 2:

                            self.add_edge(
                                visible_nodes[i], visible_nodes[j], bo=2, rigid=True
                            )

                        elif degrees[i] == 2 and degrees[j] == 3:

                            self.add_edge(
                                visible_nodes[i], visible_nodes[j], bo=2, rigid=True
                            )

                        elif degrees[i] == 2 and degrees[j] == 2:

                            self.add_edge(
                                visible_nodes[i], visible_nodes[j], bo=2, rigid=True
                            )

                        else:

                            self.add_edge(
                                visible_nodes[i],
                                visible_nodes[j],
                                bo=bo_matrix[i, j],
                                rigid=True,
                            )

                    else:
                        # Add bond order as edge weight to the graph
                        self.add_edge(
                            visible_nodes[i],
                            visible_nodes[j],
                            bo=bo_matrix[i, j],
                            rigid=True,
                        )

    def _flexify_nitroso(self):

        self.set_nodes_visible(include_all=True)

        visible_edges = self.return_visible_edges(nodes_visible=True)

        visible_nodes = self.return_visible_nodes()

        visible_elements = self.return_node_data(attribute="elements")

        degrees = self.get_degree()

        for edge in visible_edges:

            elements_in_edge = [
                visible_elements[visible_nodes.index(node)] for node in edge
            ]

            if "N" in elements_in_edge and "O" in elements_in_edge:

                oxygen_index = edge[elements_in_edge.index("O")]

                nitrogen_index = edge[elements_in_edge.index("N")]

                adjacent_to_nitrogen = self.return_adjacent_nodes(nitrogen_index)

                oxygen_degree = degrees[visible_nodes.index(oxygen_index)]

                if oxygen_degree == 1:

                    for atom in adjacent_to_nitrogen:

                        if atom != oxygen_index:

                            self.add_edge(nitrogen_index, atom, rigid=False)

    def _flexify_nitrogen_valence_5(self):

        self.set_nodes_visible(include_all=True)

        visible_nodes = self.return_visible_nodes()

        visible_elements = self.return_node_data(attribute="elements")

        valences = self.return_node_data(attribute="valences")

        for node, element, valence in zip(visible_nodes, visible_elements, valences):

            if element == "N" and valence >= 4:

                adjacent_nodes = self.return_adjacent_nodes(node)

                for atom in adjacent_nodes:

                    self.add_edge(node, atom, rigid=False)

    def _add_amide_group(self):

        self.set_nodes_visible(include_all=True)

        visible_edges = self.return_visible_edges(nodes_visible=True)

        visible_nodes = self.return_visible_nodes()

        visible_elements = self.return_node_data(attribute="elements")

        for edge in visible_edges:

            elements_in_edge = [
                visible_elements[visible_nodes.index(node)] for node in edge
            ]

            if "N" in elements_in_edge and "C" in elements_in_edge:

                carbon_index = edge[elements_in_edge.index("C")]

                nitrogen_index = edge[elements_in_edge.index("N")]

                adjacent_to_carbon = self.return_adjacent_nodes(carbon_index)

                if len(adjacent_to_carbon) == 3:

                    self.set_nodes_visible(visible_nodes=adjacent_to_carbon)

                    adjacent_nodes = self.return_visible_nodes()

                    elements_adjacent_to_carbon = self.return_node_data(
                        attribute="elements"
                    )

                    self.set_nodes_visible(include_all=True)

                    unique_elements, counts = np.unique(
                        elements_adjacent_to_carbon, return_counts=True
                    )

                    if ("O" in elements_adjacent_to_carbon) and (
                        counts[np.where(unique_elements == "O")[0][0]] == 1
                    ):

                        oxygen_index = adjacent_nodes[
                            np.where(elements_adjacent_to_carbon == "O")[0][0]
                        ]
                        bond_order_C0 = self.edges[carbon_index, oxygen_index]["bo"]

                        if bond_order_C0 == 2:

                            self.add_edge(
                                carbon_index, nitrogen_index, bo=2, rigid=False
                            )

                    elif ("S" in elements_adjacent_to_carbon) and (
                        counts[np.where(unique_elements == "S")[0][0]] == 1
                    ):

                        sulphur_index = adjacent_nodes[
                            np.where(elements_adjacent_to_carbon == "S")[0][0]
                        ]

                        bond_order_C0 = self.edges[carbon_index, sulphur_index]["bo"]

                        if bond_order_C0 == 2:

                            self.add_edge(
                                carbon_index, nitrogen_index, bo=2, rigid=False
                            )

    def _get_combinations_of_valences(
        self,
        atoms,
        atomic_valences: list,
        degrees: list,
        are_core: list,
        visible_nodes: list,
    ) -> list:
        """
        Checks for each which valences are possible (valence must be greater than actual coordination number) and creates itertools.product
        representing all atomic valence combinations.

        Returns:
            itertools.product: all atomic valence combinations
        """
        all_possible_valences = []

        for n, is_core in enumerate(are_core):

            if is_core == True:

                atom = atoms[n]

                degree = degrees[n]
                # For each atom, gets all possible valences

                possible_valences = atomic_valences[n]

                node = visible_nodes[n]

                if atom == "S":

                    if degree <= 2:

                        possible_valence_for_atom = [2]

                    elif degree == 3:

                        possible_valence_for_atom = [3]

                    elif (degree == 4) | (degree == 5) | (degree == 6):

                        possible_valence_for_atom = [6]

                    else:

                        possible_valence_for_atom = [degree]

                elif atom == "P":

                    if degree <= 3:

                        possible_valence_for_atom = [3]

                    else:

                        possible_valence_for_atom = [5]

                elif atom == "N":

                    if degree == 4:

                        possible_valence_for_atom = [4]

                    elif degree == 3:

                        adjacent_nodes = self.return_adjacent_nodes(core_nodes=node)

                        degree = self.get_degree(include_all=True)

                        elements = [
                            self.nodes(data=True)[node]["elements"]
                            for node in adjacent_nodes
                        ]

                        if not "O" in elements:
                            possible_valence_for_atom = [4, 3]
                        else:
                            oxygen_counter = 0
                            for node, element in zip(adjacent_nodes, elements):
                                if element == "O" and degree[node] == 1:
                                    oxygen_counter += 1

                            if oxygen_counter == 1:
                                possible_valence_for_atom = [5, 4]
                            elif oxygen_counter == 2:
                                possible_valence_for_atom = [6, 5]
                            elif oxygen_counter == 3:
                                possible_valence_for_atom = [6]
                            else:
                                possible_valence_for_atom = [4, 3]

                    elif degree == 2:

                        possible_valence_for_atom = [3]

                    elif degree == 1:

                        possible_valence_for_atom = [3]

                    else:

                        possible_valence_for_atom = [degree]

                elif atom == "Si":

                    if degree == 2:

                        possible_valence_for_atom = [2]

                    else:

                        possible_valence_for_atom = [4]

                # If atom has no valence of 0 and coordination number is less than maximum possible valence, then all select all possible valences.
                elif (0 not in possible_valences) & (degree <= max(possible_valences)):

                    possible_valence_for_atom = [
                        valence for valence in atomic_valences[n] if valence >= degree
                    ]
                else:
                    # If atom has valence of 0 or coordination number is greater than maximum possible valence, then only select coordination number as possible valence.
                    # This refers to being only singly bonded.
                    possible_valence_for_atom = [degree]

                all_possible_valences.append(possible_valence_for_atom)

            else:

                # For each adjacent atom, aka saturated atom, only select coordination number as possible valence
                # (we cut the graph, and its coordination number is only one while already beeing saturated)

                all_possible_valences.append([degrees[n]])

        # Creates itertools.product of all possible valences
        return np.asfortranarray(
            [list for list in itertools.product(*all_possible_valences)]
        ).T

    def _saturate_groups(
        self,
        is_unsaturated: list,
        visible_nodes: list,
        elements: list,
        degree: list,
        valences: list,
    ):
        """
        Saturates carbonyl groups and carbon atoms with a degree of 2 and their adjacent atoms to simplify bond order detection.

        Args:
            is_unsaturated (list): list of all atoms that are possibly unsaturated
            visible_nodes (list): visible nodes of the molecule
            elements (list): elements of all atoms
            degree (list): degree of all atoms
            valences (list): valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """

        unsaturated_nodes = np.argwhere(is_unsaturated == True).flatten()

        for unsaturated_node in unsaturated_nodes:
            is_unsaturated = self._saturate_carbonyl_group(
                unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            )
            if not is_unsaturated[unsaturated_node]:
                continue

            is_unsaturated = self._saturate_alkyne(
                unsaturated_node,
                elements,
                degree,
                visible_nodes,
                is_unsaturated,
                valences,
            )
            if not is_unsaturated[unsaturated_node]:
                continue

            is_unsaturated = self._saturate_terminal_alkyne(
                unsaturated_node,
                elements,
                degree,
                visible_nodes,
                is_unsaturated,
                valences,
            )
            if not is_unsaturated[unsaturated_node]:
                continue

            is_unsaturated = self._saturate_cyanide(
                unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            )

            if not is_unsaturated[unsaturated_node]:
                continue
            is_unsaturated = self._saturate_diazo(
                unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            )
            if not is_unsaturated[unsaturated_node]:
                continue

            is_unsaturated = self._saturate_isocyanate(
                unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            )
            if not is_unsaturated[unsaturated_node]:
                continue

            # is_unsaturated = self._saturate_amine_oxide(
            #    unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            # )
            # if not is_unsaturated[unsaturated_node]:
            #    continue

            is_unsaturated = self._saturate_allene(
                unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            )
            if not is_unsaturated[unsaturated_node]:
                continue
            # is_unsaturated = self._saturate_nitroxyl(
            #    unsaturated_node, elements, degree, visible_nodes, is_unsaturated
            # )
            # if not is_unsaturated[unsaturated_node]:
            #    continue

        return is_unsaturated

    def _saturate_carbonyl_group(
        self, unsaturated_node, elements, degree, visible_nodes, is_unsaturated
    ):
        """
        Assigns a double bond to a carbonyl group before bond order detection to simplify bond order detection.

        Args:
            unsaturated_node (int): Atom index of the unsaturated atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """

        if elements[unsaturated_node] == "C" and degree[unsaturated_node] == 3:
            carbon_node = visible_nodes[unsaturated_node]
            adjacent_nodes = self.return_adjacent_nodes(core_nodes=carbon_node)
            adjacent_oxygen_index = [
                i
                for i, node in enumerate(adjacent_nodes)
                if (elements[visible_nodes.index(node)] == "O")
                or (elements[visible_nodes.index(node)] == "S")
            ]

            if adjacent_oxygen_index:
                for index in adjacent_oxygen_index:
                    oxygen_node = adjacent_nodes[index]
                    oxygen_degree = degree[visible_nodes.index(oxygen_node)]

                    if oxygen_degree == 1:
                        self.add_edge(carbon_node, oxygen_node, bo=2, rigid=True)
                        is_unsaturated[unsaturated_node] = False
                        is_unsaturated[visible_nodes.index(oxygen_node)] = False

        return is_unsaturated

    def _saturate_terminal_alkyne(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
        valences: list,
    ):
        """
        Saturates a carbon atom with a degree of 2 and its adjacent atoms.

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "C" and degree[unsaturated_node] == 1:

            carbon_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=carbon_node)

            for adjacent_node in adjacent_nodes:
                try:
                    adjacent_index = visible_nodes.index(adjacent_node)
                    adjacent_degree = degree[adjacent_index]
                    adjacent_element = elements[adjacent_index]
                except ValueError:
                    continue

                if adjacent_element == "C" and adjacent_degree == 2:

                    is_unsaturated[unsaturated_node] = False

                    is_unsaturated[adjacent_index] = False

                    self.add_edge(carbon_node, adjacent_node, bo=3, rigid=True)

        return is_unsaturated

    def _saturate_alkyne(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
        valences: list,
    ):
        """
        Saturates a carbon atom with a degree of 2 and its adjacent atoms.

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "C" and degree[unsaturated_node] == 2:

            carbon_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=carbon_node)

            adjacent_degrees = []

            adjacent_elements = []

            for adjacent_node in adjacent_nodes:

                try:

                    adjacent_index = visible_nodes.index(adjacent_node)

                    adjacent_degree = degree[adjacent_index]

                    adjacent_element = elements[adjacent_index]

                    adjacent_degrees.append(adjacent_degree)

                    adjacent_elements.append(adjacent_element)

                except ValueError:

                    return is_unsaturated

            if "C" in adjacent_elements and "H" in adjacent_elements:

                adjacent_carbon_index = adjacent_elements.index("C")

                adjacent_carbon_degree = adjacent_degrees[adjacent_carbon_index]

                adjacent_hydrogen_index = adjacent_elements.index("H")

                adjacent_hydrogen_degree = adjacent_degrees[adjacent_hydrogen_index]

                if (adjacent_carbon_degree != 2) or (adjacent_hydrogen_degree != 1):

                    return is_unsaturated

                global_carbon_index = adjacent_nodes[adjacent_carbon_index]

                self.add_edge(carbon_node, global_carbon_index, bo=2, rigid=True)

                is_unsaturated[visible_nodes.index(global_carbon_index)] = False

        return is_unsaturated

    def _saturate_amine_oxide(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):

        if elements[unsaturated_node] == "O" and degree[unsaturated_node] == 1:

            central_node = visible_nodes[unsaturated_node]

            nitrogen_node = self.return_adjacent_nodes(core_nodes=central_node)[0]

            nitrogen_index = visible_nodes.index(nitrogen_node)

            adjacent_element = elements[nitrogen_index]

            if adjacent_element == "N":

                is_unsaturated[nitrogen_index] = False

                is_unsaturated[unsaturated_node] = False

                self.add_edge(central_node, nitrogen_node, bo=2, rigid=False)

                adjacent_nodes = self.return_adjacent_nodes(core_nodes=nitrogen_node)

                for adjacent_node in adjacent_nodes:

                    adjacent_index = visible_nodes.index(adjacent_node)

                    adjacent_degree = degree[adjacent_index]

                    adjacent_element = elements[adjacent_index]

                    # if adjacent_element == "C" and adjacent_degree == 3:

                    #    self.add_edge(nitrogen_node, adjacent_node, bo=2, rigid=False)

        return is_unsaturated

    def _saturate_(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a nitrogen atom with a degree of 1 and its adjacent atoms.

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "N" and degree[unsaturated_node] == 1:

            central_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=central_node)

            for adjacent_node in adjacent_nodes:

                first_adjacent_index = visible_nodes.index(adjacent_node)

                first_adjacent_degree = degree[first_adjacent_index]

                first_element = elements[first_adjacent_index]

                if first_element == "N" and first_adjacent_degree == 2:

                    is_unsaturated[unsaturated_node] = False

                    is_unsaturated[first_adjacent_index] = False

                    self.add_edge(
                        central_node,
                        visible_nodes[first_adjacent_index],
                        bo=2,
                        rigid=True,
                    )

        return is_unsaturated

    def _saturate_diazo(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a nitrogen atom with a degree of 1 and its adjacent atoms.

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "N" and degree[unsaturated_node] == 1:

            central_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=central_node)

            for adjacent_node in adjacent_nodes:

                first_adjacent_index = visible_nodes.index(adjacent_node)

                first_adjacent_degree = degree[first_adjacent_index]

                first_element = elements[first_adjacent_index]

                if first_element == "N" and first_adjacent_degree == 2:

                    is_unsaturated[unsaturated_node] = False

                    is_unsaturated[first_adjacent_index] = False

                    self.add_edge(
                        central_node,
                        visible_nodes[first_adjacent_index],
                        bo=2,
                        rigid=True,
                    )

        return is_unsaturated

    def _saturate_isocyanate(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a isocyanate group

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """

        if elements[unsaturated_node] == "C" and degree[unsaturated_node] == 2:

            carbon_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=carbon_node)

            adjacent_degrees = []

            adjacent_elements = []

            for adjacent_node in adjacent_nodes:

                try:

                    adjacent_index = visible_nodes.index(adjacent_node)

                    adjacent_degree = degree[adjacent_index]

                    adjacent_element = elements[adjacent_index]

                    adjacent_degrees.append(adjacent_degree)

                    adjacent_elements.append(adjacent_element)

                except ValueError:

                    return is_unsaturated

            if "N" in adjacent_elements and "O" in adjacent_elements:

                adjacent_oxygen_index = adjacent_elements.index("O")

                adjacent_oxygen_degree = adjacent_degrees[adjacent_oxygen_index]

                adjacent_nitrogen_index = adjacent_elements.index("N")

                adjacent_nitrogen_degree = adjacent_degrees[adjacent_nitrogen_index]

                if adjacent_oxygen_degree != 1 and adjacent_nitrogen_degree != 2:

                    return is_unsaturated

                global_oxygen_index = adjacent_nodes[adjacent_oxygen_index]

                global_nitrogen_index = adjacent_nodes[adjacent_nitrogen_index]

                self.add_edge(carbon_node, global_oxygen_index, bo=2, rigid=True)

                is_unsaturated[visible_nodes.index(global_oxygen_index)] = False

                self.add_edge(carbon_node, global_nitrogen_index, bo=2, rigid=True)

                is_unsaturated[visible_nodes.index(global_nitrogen_index)] = False

        return is_unsaturated

    def _saturate_nitroxyl(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a nitroxyl group

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """

        if elements[unsaturated_node] == "N" and degree[unsaturated_node] == 3:

            nitrogen_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=nitrogen_node)

            adjacent_degrees = []

            adjacent_elements = []

            for adjacent_node in adjacent_nodes:

                try:

                    adjacent_index = visible_nodes.index(adjacent_node)

                    adjacent_degree = degree[adjacent_index]

                    adjacent_element = elements[adjacent_index]

                    adjacent_degrees.append(adjacent_degree)
                    adjacent_elements.append(adjacent_element)

                except ValueError:

                    return is_unsaturated

            if adjacent_elements.count("O") == 1:

                adjacent_oxygen_index = adjacent_elements.index("O")

                adjacent_oxygen_degree = adjacent_degrees[adjacent_oxygen_index]

                if adjacent_oxygen_degree != 1:

                    return is_unsaturated

                global_oxygen_index = adjacent_nodes[adjacent_oxygen_index]

                self.add_edge(nitrogen_node, global_oxygen_index, bo=2, rigid=True)

                local_oxygen_index = visible_nodes.index(global_oxygen_index)

                is_unsaturated[unsaturated_node] = False

                is_unsaturated[local_oxygen_index] = False

                for adjacent_node in adjacent_nodes:

                    adjacent_index = visible_nodes.index(adjacent_node)

                    adjacent_degree = degree[adjacent_index]

                    adjacent_element = elements[adjacent_index]

                    if adjacent_element == "C" and adjacent_degree == 3:

                        is_unsaturated[unsaturated_node] = False

                        is_unsaturated[adjacent_index] = False

                        self.add_edge(nitrogen_node, adjacent_node, bo=2, rigid=True)
        return is_unsaturated

    def _saturate_cyanide(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a cyanide group

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "N" and degree[unsaturated_node] == 1:

            nitrogen_node = visible_nodes[unsaturated_node]

            adjacent_node = self.return_adjacent_nodes(core_nodes=nitrogen_node)[0]

            try:

                adjacent_index = visible_nodes.index(adjacent_node)

                adjacent_degree = degree[adjacent_index]

                adjacent_element = elements[adjacent_index]

            except ValueError:

                return is_unsaturated

            if adjacent_element == "C" and adjacent_degree == 2:

                is_unsaturated[unsaturated_node] = False

                is_unsaturated[adjacent_index] = False

                self.add_edge(nitrogen_node, adjacent_node, bo=2, rigid=True)

        return is_unsaturated

    def _saturate_allene(
        self,
        unsaturated_node: int,
        elements: list,
        degree: list,
        visible_nodes: list,
        is_unsaturated: list,
    ):
        """
        Saturates a allene group.

        Args:
            unsaturated_node (int): Atom index of the unsaturated carbon atom
            elements (list): Elements of all atoms
            degree (list): Degree of all atoms
            visible_nodes (list): Visible nodes of the molecule
            is_unsaturated (list): List of all atoms that are possibly unsaturated
            valences (list): Valences of all atoms

        Returns:
            is_unsaturated: Refined list of all atoms that are possibly unsaturated
        """
        if elements[unsaturated_node] == "C" and degree[unsaturated_node] == 2:

            central_node = visible_nodes[unsaturated_node]

            adjacent_nodes = self.return_adjacent_nodes(core_nodes=central_node)

            first_adjacent_index = visible_nodes.index(adjacent_nodes[0])

            first_adjacent_degree = degree[first_adjacent_index]

            first_element = elements[first_adjacent_index]

            second_adjacent_index = visible_nodes.index(adjacent_nodes[1])

            second_adjacent_degree = degree[second_adjacent_index]

            second_element = elements[second_adjacent_index]

            if (first_element == "C" and first_adjacent_degree == 3) and (
                second_element == "C" and second_adjacent_degree == 3
            ):

                is_unsaturated[unsaturated_node] = False

                is_unsaturated[first_adjacent_index] = False

                is_unsaturated[second_adjacent_index] = False

                self.add_edge(
                    central_node, visible_nodes[first_adjacent_index], bo=2, rigid=True
                )

                self.add_edge(
                    central_node, visible_nodes[second_adjacent_index], bo=2, rigid=True
                )

        return is_unsaturated
