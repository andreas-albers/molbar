import numpy as np
from molbar.helper.vector import Vector
import warnings


class Dihedral:

    def __init__(self) -> None:

        pass

    def _add_allene_dihedral_constraints(self) -> None:

        visible_nodes = self.return_visible_nodes()

        elements = self.return_node_data(attribute="elements")

        vsepr_classes = self.return_node_data(attribute="vsepr_classes")

        allowed_classes = set(["trigonal_planar", "bent"])

        for central_node in visible_nodes:

            if (
                elements[visible_nodes.index(central_node)] == "C"
                and vsepr_classes[visible_nodes.index(central_node)] == "linear_cn2"
            ):

                adjacent_nodes = self.return_adjacent_nodes(
                    core_nodes=[central_node], is_adjacent_visible=True
                )

                first_vsepr = vsepr_classes[visible_nodes.index(adjacent_nodes[0])]

                second_vsepr = vsepr_classes[visible_nodes.index(adjacent_nodes[1])]

                first_element = elements[visible_nodes.index(adjacent_nodes[0])]

                second_element = elements[visible_nodes.index(adjacent_nodes[1])]

                adjacent_vsepr_classes = set([first_vsepr, second_vsepr])

                if (
                    adjacent_vsepr_classes.issubset(allowed_classes)
                    and first_element == "C"
                    and second_element == "C"
                ):

                    self._add_allene_dihedral(
                        central_node, adjacent_nodes, visible_nodes
                    )
                else:

                    self.add_edge(central_node, adjacent_nodes[0], rigid=False)

                    self.add_edge(central_node, adjacent_nodes[1], rigid=False)

    def _add_allene_dihedral(
        self, central_node, core_nodes: np.ndarray, visible_nodes: np.ndarray
    ) -> None:

        first_adjacent_nodes = [
            node
            for node in self.return_adjacent_nodes(
                core_nodes=[core_nodes[0]], is_adjacent_visible=True
            )
            if node != central_node
        ]

        second_adjacent_nodes = [
            node
            for node in self.return_adjacent_nodes(
                core_nodes=[core_nodes[1]], is_adjacent_visible=True
            )
            if node != central_node
        ]

        self.set_nodes_visible(
            visible_nodes=list(core_nodes)
            + first_adjacent_nodes
            + second_adjacent_nodes
        )

        local_visible_nodes = self.return_visible_nodes()

        local_coordinates = self.return_node_data(attribute="coordinates")

        self.set_nodes_visible(visible_nodes=visible_nodes)

        if len(set(local_visible_nodes)) != len(local_visible_nodes):

            return

        edge_constraints = []

        for first_adjacent_node in first_adjacent_nodes:

            for second_adjacent_node in second_adjacent_nodes:

                dihedral = [
                    first_adjacent_node,
                    core_nodes[0],
                    core_nodes[1],
                    second_adjacent_node,
                ]

                local_nodes = [local_visible_nodes.index(node) for node in dihedral]

                if len(local_nodes) != len(set(local_nodes)):

                    raise ValueError("Dihedral constraint is not unique.")

                dihedral_coordinates = local_coordinates[local_nodes]

                real_dihedral = Vector(dihedral_coordinates).calculate_dihedral()

                deviations = {
                    90.0: np.abs(90.0 - real_dihedral),
                    -90.0: np.abs(-90.0 - real_dihedral),
                }

                ideal_dihedral = min(deviations, key=deviations.get)

                abs_deviation = np.abs(ideal_dihedral - real_dihedral)

                if abs_deviation > 45.0:

                    return

                # Order the nodes in the dihedral constraint so that the double bond nodes are ordered in ascending order.
                if dihedral[1] < dihedral[2]:

                    new_constraint = {
                        "nodes": [dihedral[0], dihedral[3]],
                        "dihedral": ideal_dihedral * np.pi / 180.0,
                    }

                else:

                    new_constraint = {
                        "nodes": [dihedral[3], dihedral[0]],
                        "dihedral": ideal_dihedral * np.pi / 180.0,
                    }

                edge_constraints.append(new_constraint)

        self.add_edge(
            dihedral[1], dihedral[2], dihedral_constraints=edge_constraints, rigid=True
        )

    def _add_db_dihedral_constraints(self) -> None:
        """
        Adds dihedral constraints for assigned double bonds if double bond nodes have allowed VSEPR classes (bent, trigonal_planar, t_shaped).

        """

        allowed_vsepr_classes = ["bent", "trigonal_planar", "t_shaped"]

        visible_nodes = self.return_visible_nodes()

        visible_edges = self.return_visible_edges()

        # Get all double bonds
        double_bonds = self.return_edges(attributes=["bo", "rigid"], values=[2, True])

        triple_bonds = self.return_edges(attributes="bo", values=3)

        cycle_ids = self.return_edge_data(attribute="cycle_id")

        for edge in triple_bonds:

            self.add_edge(edge[0], edge[1], rigid=True)

        # Get all vsepr classes
        vsepr_classes = self.return_node_data(
            attribute="vsepr_classes", include_all=True
        )

        # Get all extra constraints defined by the input file
        if self.constraints.get("constraints") != None:

            extra_constraints = self.constraints["constraints"].get("dihedrals")

        # If no extra constraints are defined, set extra_constraints to an empty dictionary
        else:

            extra_constraints = {}

        # Iterate over all double bonds
        for edge in double_bonds:

            # If the double bond is already constrained by an extra constraint, skip it. Prevents double constraints.
            # It is possible to do it like this, because the extra constraints dictionary is defined in a way that
            # the edge node with the smaller node id is the key and the value is a list of all constraints for this central node.

            # Check if the ith_node and the jth_node are part of the same cycle. If they are part of the same cycle, the dihedral constraint is not added.

            if extra_constraints and (extra_constraints.get(edge[0]) != None):

                continue

            elif cycle_ids[visible_edges.index(edge)] is not None:

                continue

            # Check if the vsepr classes of the double bond nodes are allowed
            elif (vsepr_classes[edge[0]] in allowed_vsepr_classes) and (
                vsepr_classes[edge[1]] in allowed_vsepr_classes
            ):

                self._add_db_dihedral(edge)

                self.set_nodes_visible(visible_nodes=visible_nodes)

            elif (
                vsepr_classes[edge[0]] in allowed_vsepr_classes
                and vsepr_classes[edge[1]] == "linear_cn2"
            ) or (
                vsepr_classes[edge[1]] in allowed_vsepr_classes
                and vsepr_classes[edge[0]] == "linear_cn2"
            ):

                self.add_edge(edge[0], edge[1], rigid=True)

                continue

            else:
                self.add_edge(edge[0], edge[1], rigid=False)

        # Add extra constraints

        self._add_db_extra_constraint(extra_constraints)

    def _add_db_extra_constraint(self, extra_constraints: dict) -> None:
        """
        Adds extra constraints defined by the input file for double bonds.
        extra_constraints is a dictionary with the central node as key and a list of all constraints for this central node as value.
        {node1: [{"nodes": [node1, node2], "dihedral": dihedral}, ...], node2:...
        """

        # Iterate over all extra constraints

        for central_node, constraints in extra_constraints.items():

            # Iterate over all constraints for the central node
            for constraint in constraints:

                dihedral = [
                    constraint["nodes"][0],
                    central_node,
                    constraint["nodes"][1],
                    constraint["nodes"][2],
                ]

                edge_constraints = [
                    {
                        "nodes": [dihedral[0], dihedral[3]],
                        "dihedral": constraint["dihedral"],
                    }
                ]

                rcovij = (
                    self.nodes(data=True)[central_node]["rcov"]
                    + self.nodes(data=True)[constraint["nodes"][1]]["rcov"]
                )

                # Add constraint to the graph, bond order is set to 2, so that is considered in the fragmentation scheme to be rigid.
                self.add_edge(
                    dihedral[1],
                    dihedral[2],
                    dihedral_constraints=edge_constraints,
                    bo=1,
                    dij=rcovij,
                    rigid=True,
                )

    def _add_db_dihedral(self, edge: list) -> None:
        """
        Adds four dihedral constraints for a single double bond.

        Args:
            edge (list): List of two nodes that are connected by a double bond.
        """

        # Get all adjacent nodes of the double bond nodes
        ith_adjacent_nodes = [
            node
            for node in self.return_adjacent_nodes(
                int(edge[0]), is_adjacent_visible=True
            )
            if node not in edge
        ]

        jth_adjacent_nodes = [
            node
            for node in self.return_adjacent_nodes(
                int(edge[1]), is_adjacent_visible=True
            )
            if node not in edge
        ]

        # Only consider the double bond atoms and their adjacent atoms in the following.
        self.set_nodes_visible(
            visible_nodes=list(edge) + ith_adjacent_nodes + jth_adjacent_nodes
        )

        visible_nodes = self.return_visible_nodes()

        # Get the coordinates of the double bond atoms and their adjacent atoms
        coordinates = self.return_node_data(attribute="coordinates")

        is_metal = self.return_node_data(attribute="is_metal")

        if any(is_metal):
            self.add_edge(edge[0], edge[1], rigid=False)
            return

        # to_be_excluded = self._exclude_NH_groups_from_db_dihedral(
        #    edge, ith_adjacent_nodes, jth_adjacent_nodes, visible_nodes
        # )

        # if to_be_excluded:

        #    return

        # Iterate over all adjacent nodes of the first double bond node, so that each double bond is constrained by 4 dihedrals
        for ith_node in ith_adjacent_nodes:

            for jth_node in jth_adjacent_nodes:

                if jth_node > ith_node:

                    dihedral = [ith_node, edge[0], edge[1], jth_node]

                else:

                    dihedral = [jth_node, edge[1], edge[0], ith_node]

                self._add_single_db_dihedral(dihedral, coordinates, visible_nodes)

                self.set_nodes_visible(
                    visible_nodes=list(edge) + ith_adjacent_nodes + jth_adjacent_nodes
                )

    def _exclude_NH_groups_from_db_dihedral(
        self,
        edge: list,
        ith_adjacent_nodes: np.ndarray,
        jth_adjacent_nodes: np.ndarray,
        visible_nodes: np.ndarray,
    ) -> bool:
        """
        Checks whether the double bond is part of a NH group. If it is part of a NH group, the double bond is excluded from the dihedral constraint.

        Args:
            edge (list): List of two nodes that are connected by a double bond.
            ith_adjacent_nodes (np.ndarray): adjacent nodes to first double bond node
            jth_adjacent_nodes (np.ndarray): adjacent nodes to second double bond node
            visible_nodes (np.ndarray): List of all visible nodes.

        Returns:
            bool: True if the double bond is part of a NH group, False otherwise.
        """

        elements = self.return_node_data(attribute="elements")

        ith_element = elements[visible_nodes.index(edge[0])]

        jth_element = elements[visible_nodes.index(edge[1])]

        if (
            len(ith_adjacent_nodes) == 1
            and elements[visible_nodes.index(ith_adjacent_nodes[0])] == "H"
            and ith_element == "N"
        ):

            return True

        elif (
            len(jth_adjacent_nodes) == 1
            and elements[visible_nodes.index(jth_adjacent_nodes[0])] == "H"
            and jth_element == "N"
        ):

            return True

        return False

    def _add_single_db_dihedral(
        self, dihedral: list, coordinates: np.ndarray, visible_nodes: np.ndarray
    ) -> None:
        """
        Adds a single dihedral constraint for a double bond.

        Args:
            dihedral (list): List of four nodes that are connected by a dihedral.
            coordinates (np.ndarray): Coordinates of all 4 nodes in the dihedral.
            visible_nodes (list): List of all visible nodes.
        """

        local_nodes = [visible_nodes.index(node) for node in dihedral]
        # Calculate the real dihedral angle
        real_dihedral = Vector(coordinates[local_nodes]).calculate_dihedral()

        # Calculate the deviation from the all three possible ideal dihedral angles (-180°, 0, 180°)
        deviations = {
            -180.0: np.abs(-180.0 - real_dihedral),
            0.0: np.abs(0.0 - real_dihedral),
            180.0: np.abs(180.0 - real_dihedral),
        }

        # Get the ideal dihedral angle based on the smallest deviation.
        ideal_dihedral = min(deviations, key=deviations.get)
        # Calculate the absolute deviation
        abs_deviation = np.abs(ideal_dihedral - real_dihedral)

        # Get existing dihedral constraints for the edge.
        edge_constraints = self.get_edge_data(
            dihedral[1], dihedral[2], default=None
        ).get("dihedral_constraints")

        # Order the nodes in the dihedral constraint so that the double bond nodes are ordered in ascending order.
        if dihedral[1] < dihedral[2]:

            new_constraint = {
                "nodes": [dihedral[0], dihedral[3]],
                "dihedral": ideal_dihedral * np.pi / 180.0,
            }

        else:

            new_constraint = {
                "nodes": [dihedral[3], dihedral[0]],
                "dihedral": ideal_dihedral * np.pi / 180.0,
            }

        # If the absolute deviation is larger than 20°, the dihedral constraint is added to the graph and the bond order is set to 1.
        if abs_deviation > 90.0:

            self.add_edge(dihedral[1], dihedral[2], rigid=False)

            return

        # If there are no existing dihedral constraints, add the new constraint to the graph.
        if edge_constraints is None:

            self.add_edge(
                dihedral[1],
                dihedral[2],
                dihedral_constraints=[new_constraint],
                rigid=True,
            )

            self.set_nodes_visible(visible_nodes=[dihedral[1], dihedral[2]])

            self.add_node_data(attribute="in_accepted_db", new_data=[True] * 2)

        else:

            edge_constraints.append(new_constraint)

            self.add_edge(
                dihedral[1],
                dihedral[2],
                dihedral_constraints=edge_constraints,
                rigid=True,
            )

    def _add_cycle_dihedral_constraints(self) -> None:
        """
        Adds dihedral constraints for all cycles in the molecule.
        """

        # Get the cycle ids all edges are part of.
        all_cycle_ids = self.return_edge_data(
            attribute="cycle_id", nodes_visible=True, include_all=False
        )

        # Get all edges that are visible
        edges = self.return_visible_edges(nodes_visible=True)

        # Get all unique cycle ids
        cycles = list(
            set(
                [
                    cycle_id
                    for cycle_ids in all_cycle_ids
                    if cycle_ids is not None
                    for cycle_id in cycle_ids
                ]
            )
        )

        # Get all edges that are part of a cycle
        edges_in_cycles = [
            [
                edge
                for i, edge in enumerate(edges)
                if (all_cycle_ids[i] is not None) and (cycle in all_cycle_ids[i])
            ]
            for cycle in cycles
        ]

        # Get all nodes that are part of a cycle
        nodes_in_cycle = [
            set([node for edge in edges_in_cycle for node in edge])
            for edges_in_cycle in edges_in_cycles
        ]

        # Iterate over all cycles via their edges and nodes.
        for edges_in_cycle, cycle_nodes in zip(edges_in_cycles, nodes_in_cycle):

            if len(cycle_nodes) > 3:

                self._add_cycle_dihedral(edges_in_cycle, cycle_nodes)

    def _add_cycle_dihedral(self, edges_in_cycle: list, nodes_in_cycle: list) -> None:
        """
        Adds dihedral constraints for a single cycle.

        Args:
            edges_in_cycle (list): List of edges that are part of the cycle.
            nodes_in_cycle (list): List of nodes that are part of the cycle.
        """

        # Only consider the nodes in the cycle in the following.
        self.set_nodes_visible(visible_nodes=list(nodes_in_cycle))

        # Only consider the edges in the cycle in the following except the first one.
        self.set_edges_visible(visible_edges=edges_in_cycle[1:])

        # Get the shortest path between the first and the second node in the cycle.
        starting_edge = edges_in_cycle[0]

        # Get the path the defines the cycle.
        shortest_path = self.get_shortest_path(starting_edge[0], starting_edge[1])

        # Iterate over all paths in the cycle so that every node is the first node in the path once.
        for _ in range((len(nodes_in_cycle))):

            # Shift the shortest path by one node to get the next path.
            shortest_path = shortest_path[1:] + [shortest_path[0]]

            # Add the dihedral constraint for the cycle.
            self._add_single_cycle_dihedral(shortest_path)

    def _add_single_cycle_dihedral(self, path: list) -> None:
        """
        Adds a single dihedral constraint based on a path of 4 nodes in a cycle.

        Args:
            path (list): List of 4 nodes that are connected by a dihedral.
        """

        # Get the cycles the first edge is part of.
        first_cycles = self.get_edge_data(path[0], path[1], default=None).get(
            "cycle_id"
        )

        # Get the cycles the second edge is part of.
        second_cycles = self.get_edge_data(path[2], path[3], default=None).get(
            "cycle_id"
        )

        ## If the first and the second edge are part of the same cycles, the dihedral constraint is added to the graph.
        # if len(first_cycles) == len(second_cycles):

        #    if set(first_cycles) != set(second_cycles):

        ##        return
        # else:

        #    return

        # If the bond order of the core edge is larger than 1 it is checked whether the core edge is already constrained by a double bond dihedral constraint.
        if self.get_edge_data(path[1], path[2], default=None).get("bo") > 1:

            self.set_nodes_visible(visible_nodes=[path[1], path[2]])

            if all(self.return_node_data(attribute="in_accepted_db", include_all=True)):

                return

        # Get existing dihedral constraints for the edge.
        edge_constraints = self.get_edge_data(path[1], path[2], default=None).get(
            "dihedral_constraints"
        )

        # Order the nodes in the dihedral constraint so that the first and the second node are ordered in ascending order.
        if path[1] > path[2]:

            new_constraint = {"nodes": [path[3], path[0]], "dihedral": 0.0}

        else:

            new_constraint = {"nodes": [path[0], path[3]], "dihedral": 0.0}

        # If there are no existing dihedral constraints, add the new constraint to the graph.
        if edge_constraints is None:

            self.add_edge(path[1], path[2], dihedral_constraints=[new_constraint])

        else:

            edge_constraints.append(new_constraint)

            self.add_edge(path[1], path[2], dihedral_constraints=edge_constraints)
