import os
import json
import math
import numpy as np
# import warnings as warning
from ase import Atoms
# from molbar.io.filereader import FileReader
from dscribe.descriptors import SOAP
from molbar.helper.vector import Vector
from scipy.spatial.transform import Rotation
from itertools import permutations

with open(
    os.path.join(os.path.dirname(__file__), "..", "..", "data", "soap.json")
) as json_file:
    ideal_soaps = json.load(json_file)


with open(
    os.path.join(
        os.path.dirname(__file__), "..", "..", "data", "ideal_geometries.json"
    )
) as json_file:

    ideal_geometries = json.load(json_file)

class VSEPR:

    def classify_nodes_geometry(
        self, is_adjacent_visible=False, include_all=False
    ) -> None:
        """
        Classifies the geometry of all nodes of the molecule according to VSEPR theory.

        Args:
            is_adjacent_visible (bool, optional): Whether to consider only visible adjacent nodes for the classification. Defaults to False.
            include_all (bool, optional): Whether to classify all nodes in the classification. Otherwise only visible nodes are classified. Defaults to False.

        Returns:
            None
        """

        if self.is_2D:

            self.classify_nodes_geometry_from2D()

            return

        if include_all:

            nodes = [node for node, data in self.nodes(data=True)]

        else:

            nodes = [
                node
                for node, data in self.nodes(data=True)
                if data.get("visible") == True
            ]

        # Classify all nodes
        vsepr_classes = [
            self._classify_node_geometry(
                node, nodes, is_adjacent_visible=is_adjacent_visible
            )
            for node in nodes
        ]

        # Add classification to nodes
        self.add_node_data(attribute="vsepr_classes", new_data=vsepr_classes)

    def _classify_node_geometry(
        self, central_node: int, initial_visible_nodes: list, is_adjacent_visible=False
    ) -> str:
        """
        Classifies the geometry of one node of the molecule according to VSEPR theory.

        Args:
            node (int): node id
            nodes (list): list of all visible nodes
            is_adjacent_visible (bool, optional): Whether to consider only visible adjacent nodes for the classification. Defaults to False.

        Returns:
            str: VSEPR class of the node.
        """
        # Get adjacent nodes to central node
        adjacent_nodes = self.return_adjacent_nodes(
            central_node, is_adjacent_visible=is_adjacent_visible
        )

        self.set_nodes_visible(visible_nodes=[central_node])

        central_cycle_ids = self.return_node_data(attribute="cycle_id")

        # Check if node is in cycle, has influence on VSEPR class, e.g. usally only linear geometry is allowed, but is allowed to be bent if in cycle
        if self.return_node_data(attribute="cycle_id")[0] is not None:

            is_in_cycle = True

            n_cycles = len(central_cycle_ids[0])

            smallest_cycle_size = self.return_node_data(
                attribute="smallest_cycle_size"
            )[0]

        else:

            is_in_cycle = False

            n_cycles = 0

            smallest_cycle_size = 0

        nuclear_charge = self.return_node_data(attribute="atomic_numbers")[0]

        valence = self.return_node_data(attribute="valences")[0]

        valence_electrons = self.return_node_data(attribute="valence_electrons")[0]

        formal_charge = self.return_node_data(attribute="formal_charge")[0]

        self.set_nodes_visible(visible_nodes=[central_node] + adjacent_nodes)

        visible_nodes = self.return_visible_nodes()
        # Get index of central node in visible nodes
        core_index = visible_nodes.index(central_node)

        # Get coordinates of all visible nodes
        coordinates = self.return_node_data(attribute="coordinates")

        elements = self.return_node_data(attribute="elements")

        # Check if one of the adjacent nodes or central node is a metal
        is_metal = any(self.return_node_data(attribute="is_metal"))

        n_atoms = self.return_n_atoms()

        self.set_nodes_visible(visible_nodes=initial_visible_nodes)

        # Rescale coordinates to unit vectors, so that there is no distance dependency
        coordinates = np.array(
            [Vector(coordinates[[core_index, i]]).unit_vector() for i in range(n_atoms)]
        )

        # Calculate SOAP for central node
        soap = self._get_soap(coordinates, n_atoms, core_index)

        non_core_indices = [i for i in range(n_atoms) if i != core_index]

        if n_atoms == 4 or n_atoms == 5:
            try:
                volumes = []
                for i in non_core_indices:
                    for j in non_core_indices:
                        for k in non_core_indices:
                            if k > j and j > i:
                                volumes.append(
                                    np.abs(np.linalg.det(coordinates[[i, j, k]]))
                                )
                volume = np.mean(volumes)
            except np.linalg.LinAlgError:
                volume = 0.0
        else:
            volume = 0.0

        if self.formal_charges:
            vsepr_names = {
                'AX1': ['linear_cn1'],
                'AX1E0': ['linear_cn1'],
                'AX1E1': ['linear_cn1'],
                'AX1E2': ['linear_cn1'],
                'AX1E3': ['linear_cn1'],
                'AX1E4': ['linear_cn1'],
                'AX1E5': ['linear_cn1'],
                'AX1E6': ['linear_cn1'],
                'AX2': ['linear_cn2', 'bent', 'bent_tetrahedral'],
                'AX2E0': ['linear_cn2'],
                'AX2E1': ['bent'],
                'AX2E2': ['bent_tetrahedral'],
                'AX3': ['trigonal_planar', 'trigonal_bent', 't_shaped'],
                'AX3E0': ['trigonal_planar'],
                'AX3E1': ['trigonal_bent'],
                'AX3E2': ['trigonal_planar', 't_shaped'], ###
                'AX4': ['tetrahedral', 'trigonal_pyramidal', 'seesaw', 'seesaw_rectangular', 'square_planar'],
                'AX4E0': ['tetrahedral'],
                'AX4E1': ['trigonal_pyramidal', 'seesaw'], ###
                'AX4E2': ['square_planar', 'seesaw_rectangular'], ###
                'AX5': ['trigonal_bipyramidal', 'square_pyramidal'],
                'AX5E0': ['trigonal_bipyramidal'],
                'AX5E1': ['square_pyramidal'],
                'AX6': ['octahedral', 'pentagonal_pyramidal'],
                'AX6E0': ['octahedral'],
                'AX6E1': ['pentagpnal_pyramidal'],
                'AX7': ['pentagpnal_bipyramidal'],
                'AX7E0':['pentagpnal_bipyramidal'],
                'AX8': ['square_antiprismatic'],
                'AX8E0':['square_antiprismatic']
            }

            lone_pairs = math.floor(0.5*(valence_electrons-formal_charge-valence))
            lone_pairs = int((abs(lone_pairs)+lone_pairs)/2) # in case a negative value is obtained, lone_pairs is set to zero

            # ToDo: Check for aromatic rings: lone pair of hetero atoms in aromatic rings is part of pi-system and does not contribute to VSEPR polyeder.
            #       This leads to wrong VSEPR structures if ring_hetero_planarisation=False.
            #       Current Workaround: formal charge of aromatic hetero atoms is raised by 2 by the user, to remove the lone pair from vsepr classification

            if len(vsepr_names[f'AX{n_atoms-1}E{lone_pairs}']) == 1:
                vsepr_class = vsepr_names[f'AX{n_atoms-1}E{lone_pairs}'][0]
            else:
                possible_classes = vsepr_names[f'AX{n_atoms-1}E{lone_pairs}']
                vsepr_class = self._classify_angle(coordinates, core_index, non_core_indices, possible_classes)

            central_element = elements[core_index]

            # apply chirality settings:
            if self.stereo_settings['amine_inversion']:
                if (central_element == "N" and vsepr_class == "trigonal_bent"):
                    vsepr_class = "trigonal_planar"
            if self.stereo_settings['ring_hetero_planarisation']:
                if (is_in_cycle == True and vsepr_class == "trigonal_bent"):
                    vsepr_class = "trigonal_planar"
            if self.stereo_settings['imine_inversion']:
                if (central_element == "N" and is_in_cycle == False and vsepr_class == "bent"):
                    vsepr_class = "linear_cn2"
            if not self.stereo_settings['ring_planarisation']:
                raise ValueError("ring_planarisation == False is not jet supported.")

        else:
            # Classify node according to VSEPR theory by comparing SOAP to ideal SOAPs
            vsepr_class = self._classifiy_soap(
                soap,
                elements[core_index],
                nuclear_charge,
                valence,
                n_atoms,
                list(elements),
                is_metal,
                is_in_cycle,
                n_cycles,
                smallest_cycle_size,
                volume,
            )

        return vsepr_class

    def classify_nodes_geometry_from2D(self):
        visible_nodes = self.return_visible_nodes()
        elements = self.return_node_data(attribute="elements")
        degrees = self.get_degree()

        vsepr_classes = []

        # for uncharged species:
        for node, element, degree in zip(visible_nodes, elements, degrees):
            self.set_nodes_visible(visible_nodes=[node])

            # Check if node is in cycle, has influence on VSEPR class, e.g. usally only linear geometry is allowed, but is allowed to be bent if in cycle
            if self.return_node_data(attribute="cycle_id")[0] is not None:
                is_in_cycle = True

            else:
                is_in_cycle = False

            ### 0 binding partners ####

            if degree == 0:
                vsepr_classes.append("dot")

            ### 1 bidning partner ###

            elif degree == 1:
                vsepr_classes.append("linear_cn1")

            ### 2 binding partners ###

            # If node is in cycle, it is allowed to be bent, otherwise it is linear, bent contains no important information
            elif (degree == 2) and (is_in_cycle == False):
                vsepr_classes.append("linear_cn2")

            elif (degree == 2) and (is_in_cycle == True):
                vsepr_classes.append("bent")

            ### 3 binding partners ###

            elif (degree == 3) and (element == "P"):
                vsepr_classes.append("trigonal_bent")

            elif (degree == 3) and (element != "P"):
                vsepr_classes.append("trigonal_planar")

            ### 4 binding partners ###

            elif degree == 4:
                vsepr_classes.append("tetrahedral")

            ### 5 binding partners ###

            elif degree == 5:
                vsepr_classes.append("trigonal_bipyramidal")

            ### 6 binding partners ####

            elif degree == 6:
                vsepr_classes.append("octahedral")

            ### more binding partners ####

            else:
                vsepr_classes.append("unknown")

        self.set_nodes_visible(visible_nodes=visible_nodes)

        self.add_node_data(attribute="vsepr_classes", new_data=vsepr_classes)

    def _classifiy_soap(
        self,
        soap: list,
        central_element: str,
        nuclear_charge: int,
        valence: int,
        n_atoms: int,
        all_elements: list,
        adjacent_metal=False,
        is_in_cycle=False,
        n_cycles=0,
        smallest_cycle_size=0,
        volume=0.0,
    ):
        """
        Classifies the geometry of one node of the molecule according to VSEPR theory based on the SOAP of the node.

        Args:
            soap (list): SOAP of the node
            central_element (str): element of the central node
            nuclear_charge (int): nuclear charge of the central node
            n_atoms (int): number of nodes in local geometry
            adjacent_metal (bool, optional): Whether one of the adjacent nodes is a metal. Defaults to False.
            is_in_cycle (bool, optional): Whether the node is in a cycle. Defaults to False.


        Returns:
            str: VSEPR class of the node.
        """
        if n_atoms == 1:

            return "dot"

        elif n_atoms == 2:

            return "linear_cn1"

        ################## n_adj=3 ##################

        # If node is in cycle, it is allowed to be bent, otherwise it is linear, bent contains no important information
        elif (
            (n_atoms == 3)
            and central_element not in {"N", "C", "B"}
            and (is_in_cycle == False)
            and (valence != 2)
        ):
            # non_core_coordinates = coordinates[[x for x in range(3) if x != core_index], :]
            # angle = self._get_bond_angle(coordinates[core_index, :], non_core_coordinates[0], non_core_coordinates[1])
            # if angle > 160:
            #     return "linear_cn2"
            # else:
            #     return "bent"
            return self._get_best_fit(soap, ["bent", "linear_cn2"])

        elif (n_atoms == 3) and (central_element == "B"):

            return "linear_cn2"

        elif (
            (n_atoms == 3)
            and (central_element != "N" and central_element != "C")
            and (is_in_cycle == False)
            and (valence == 2)
        ):

            return "linear_cn2"

        elif (
            (n_atoms == 3)
            and (central_element == "C")
            and (is_in_cycle == False)
            and (valence == 4)
        ):

            return "linear_cn2"

        elif (
            (n_atoms == 3)
            and (central_element == "N")
            and (is_in_cycle == False)
            and (valence == 4)
        ):

            return "linear_cn2"

        elif (
            (n_atoms == 3)
            and (central_element == "N")
            and (is_in_cycle == False)
            and (valence != 4)
            and ("H" not in all_elements)
        ):

            return "bent"

        elif (
            (n_atoms == 3)
            and (central_element == "N")
            and (is_in_cycle == False)
            and (valence != 4)
            and ("H" in all_elements)
        ):

            return "bent"

        elif (n_atoms == 3) and (is_in_cycle == True):

            return "bent"

        ################## n_adj=4 ##################

        # If element is not C, N or O, it is allowed to be trigonal planar, trigonal pyramidal or t-shaped, otherwise it assigned via hacks.
        elif (n_atoms == 4) and (central_element == "Cl") and (central_element == "Xe"):

            return "t_shaped"

        elif (n_atoms == 4) and (nuclear_charge >= 8):

            if volume < 0.1:

                return "trigonal_planar"

            elif volume > 0.1:

                return "trigonal_bent"

        elif (
            (n_atoms == 4)
            and (central_element == "N")
            and (is_in_cycle)
            and (valence == 3)
            and ("H" not in all_elements)
            and (smallest_cycle_size == 3)
        ):

            return "trigonal_bent"

        elif (n_atoms == 4) and (central_element == "N"):

            return "trigonal_planar"

        elif (n_atoms == 4) and (nuclear_charge < 8):

            return "trigonal_planar"

        ################## n_adj=5 ##################

        elif (
            (n_atoms == 5)
            and (central_element == "C")
            and all_elements.count("C") == 1
            and all_elements.count("O") == 2
            and all_elements.count("H") == 1
            and (adjacent_metal == True)
        ):

            return "formate"

        elif (
            (n_atoms == 5)
            and (central_element == "C")
            and all_elements.count("C") >= 3
            and (adjacent_metal == True)
            and (is_in_cycle == True)
        ):

            return "ferrocene"

        elif n_atoms == 5:

            if volume < 0.25:

                return "square_planar"

            elif volume > 0.25:

                return "tetrahedral"

        ################## n_adj=6 ##################

        elif n_atoms == 6:

            return self._get_best_fit(
                soap, ["trigonal_bipyramidal", "square_pyramidal"]
            )

        ################## n_adj>6 ##################

        elif n_atoms == 7:

            return "octahedral"

        elif n_atoms == 8:

            return "pentagonal_bipyramidal_formate"

        elif n_atoms == 9:

            return "square_antiprismatic"

        else:

            return "unknown"

    def _classify_angle(self, coords, core_index, non_core_indices, possible_classes):
        ideal_angles = {
            'trigonal_planar': [120, 120, 120],
            't_shaped': [90, 90, 180],
            'trigonal_pyramidal': [90, 90, 90, 120, 120, 120],
            'seesaw': [90, 90, 90, 90, 120, 180],
            'square_planar': [90, 90, 90, 90, 180, 180],
            'seesaw_rectangular': [90, 90, 90, 90, 90, 180]
        }

        valence_angles = []
        for idx, i in enumerate(non_core_indices):
            for j in non_core_indices[idx:]:
                if i != j:
                    valence_angles.append(self._get_bond_angle(coords[core_index,:], coords[i,:], coords[j,:]))
        valence_angles.sort()
        
        sum_of_deviations = {}
        for vsepr_class in possible_classes:
            sum_of_deviations[vsepr_class] = np.sum(np.abs(np.array(valence_angles) - np.array(ideal_angles[vsepr_class])))

        sum_of_deviations = {k: v for k, v in sorted(sum_of_deviations.items(), key=lambda item: item[1])}
        best_fit = next(iter(sum_of_deviations.items()))
        return best_fit[0]

    def _get_bond_angle(self, s, a, b):
        """
        Calculate bond angle between atoms a-s-b.

        Parameters
        ----------
        s, a, b: np.array
            Coordinates of atoms s, a, and b.

        Return
        ------
        angle: float
        """
        u = np.array(a) - np.array(s)
        v = np.array(b) - np.array(s)
        cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        if cos > 1:
            cos = 1
        if cos < -1:
            cos = -1
        angle = np.degrees(np.arccos(cos))
        return angle

    def _get_best_fit(self, soap, possible_classes):

        # Calculate deviations from ideal local geometries
        deviations = {
            key: np.linalg.norm((self.ideal_soaps[key] - soap))
            for key in possible_classes
        }

        # Return geometry class with smallest deviation
        chosen_class = min(deviations, key=deviations.get)

        return chosen_class

    def _get_soap(self, coordinates, n_atoms, core_index):

        elements = ["H"] * n_atoms

        elements[core_index] = "C"

        soap = SOAP(
            species=["C", "H"],
            periodic=False,
            r_cut=3.0,
            n_max=2,
            l_max=1,
        )

        ase_atoms = Atoms("".join(elements), positions=coordinates)

        return soap.create(ase_atoms, centers=[core_index])

    def constrain_nodes_angles(self, is_adjacent_visible=False, include_all=False):
        if include_all:

            visible_nodes = self.nodes()

        else:

            visible_nodes = self.return_visible_nodes()

        self.set_nodes_visible(visible_nodes=visible_nodes)

        vsepr_classes = self.return_node_data(
            attribute="vsepr_classes", include_all=include_all
        )

        angle_data = []

        if self.constraints.get("constraints") != None:

            extra_constraints = self.constraints["constraints"].get("angles")

        else:

            extra_constraints = {}

        # Get new vsepr classes as they might have changed _constrain_node_angle

        self.set_nodes_visible(visible_nodes=visible_nodes)

        vsepr_classes = self.return_node_data(
            attribute="vsepr_classes", include_all=include_all
        )

        for node, vsepr_class in zip(visible_nodes, vsepr_classes):

            node_angle = self._constrain_node_angle(
                node,
                visible_nodes,
                vsepr_class,
                is_adjacent_visible=is_adjacent_visible,
                include_all=include_all,
            )

            if extra_constraints and extra_constraints.get(node) != None:

                for extra_constraint in extra_constraints[node]:

                    node_angle += [extra_constraint]

                filtered_angles = {}

                for data_dict in node_angle:
                    nodes_key = tuple(data_dict["nodes"])
                    filtered_angles[nodes_key] = data_dict

                # Extract the filtered list of dictionaries containing the latest occurrences
                node_angle = list(filtered_angles.values())

            angle_data.append(node_angle)

        self.set_nodes_visible(visible_nodes=visible_nodes)

        self.add_node_data(attribute="angle_constraints", new_data=angle_data)

        dihedral_data = [
            self._constrain_node_dihedral(
                node,
                visible_nodes,
                vsepr_class,
                is_adjacent_visible=is_adjacent_visible,
                include_all=include_all,
            )
            for node, vsepr_class in zip(visible_nodes, vsepr_classes)
        ]

        self.set_nodes_visible(visible_nodes=visible_nodes)

        self.add_node_data(attribute="dihedral_constraints", new_data=dihedral_data)


    def _constrain_node_angle(
        self, node, nodes, vsepr_class, is_adjacent_visible=False, include_all=False
    ):

        if (
            vsepr_class == "unknown"
            or vsepr_class == "dot"
            or vsepr_class == "linear_cn1"
        ):

            return []

        if vsepr_class in [
            "bent",
            "bent_tetrahedral",
            "trigonal_planar",
            "trigonal_bent",
            "tetrahedral",
        ]:

            return self._constrain_simple_node_angles(
                node, nodes, vsepr_class, is_adjacent_visible, include_all
            )

        self.set_nodes_visible(visible_nodes=nodes)

        adjacent_atoms = self.return_adjacent_nodes(
            node, is_adjacent_visible=is_adjacent_visible, include_all=include_all
        )

        self.set_nodes_visible(visible_nodes=[node] + adjacent_atoms)

        coordinates = self.return_node_data(attribute="coordinates")

        n_atoms = self.return_n_atoms()

        visible_nodes = self.return_visible_nodes()

        core_index = visible_nodes.index(node)

        if vsepr_class not in ["bent", "bent_tetrahedral", "linear_cn2"]:

            best_order = self._get_best_order_for_match(
                n_atoms, coordinates, core_index, vsepr_class
            )

        else:

            best_order = [core_index] + [
                i for i in range(0, n_atoms) if i != core_index
            ]

        ideal_angles = self.ideal_geometries[vsepr_class]["angles"]

        final_angles = []

        for ideal_angle in ideal_angles:

            angle_nodes = [ideal_angle["nodes"][0], ideal_angle["nodes"][2]]

            angle_nodes_global = [
                visible_nodes[best_order[node]] for node in angle_nodes
            ]

            final_angles.append(
                {
                    "nodes": angle_nodes_global,
                    "angle": ideal_angle["angle"] * np.pi / 180.0,
                }
            )

        return final_angles

    def _get_best_order_for_match(self, n_atoms, coordinates, core_index, vsepr_class):

        coordinates = np.array(
            [
                Vector(coordinates[[core_index, i]]).unit_vector()
                for i in range(0, n_atoms)
            ]
        )

        # ideal_file_path = os.path.join(
        #     os.path.dirname(__file__),
        #     "..",
        #     "..",
        #     "data",
        #     "ideal_geometries",
        #     f"{vsepr_class}.xyz",
        # )

        # _, ideal_coordinates, _ = FileReader(ideal_file_path).read_file()

        ideal_coordinates = np.array(self.ideal_geometries[vsepr_class]["coords"])

        order = [i for i in range(0, n_atoms) if i != core_index]

        shuffled = list(permutations(order))

        rmsds = []

        for permutation in shuffled:

            coordinates_permuted = [coordinates[core_index]] + [
                coordinates[i] for i in permutation
            ]

            _, rmsd = Rotation.align_vectors(coordinates_permuted, ideal_coordinates)

            rmsds.append(rmsd)

        best_order_index = np.argmin(rmsds)

        return [core_index] + list(shuffled[best_order_index])

    def _constrain_simple_node_angles(
        self, node, nodes, vsepr_class, is_adjacent_visible, include_all
    ):

        self.set_nodes_visible(visible_nodes=node)

        cycle_size = self.return_node_data("smallest_cycle_size")[0]

        cycle_ids_of_node = self.return_node_data(attribute="cycle_id")[0]

        central_node_in_cycle = cycle_ids_of_node is not None

        if central_node_in_cycle:

            n_cycles = len(cycle_ids_of_node)

        else:

            n_cycles = 0

        self.set_nodes_visible(visible_nodes=nodes)

        adjacent_atoms = self.return_adjacent_nodes(
            node, is_adjacent_visible=is_adjacent_visible, include_all=include_all
        )

        self.set_nodes_visible(visible_nodes=[node] + adjacent_atoms)

        visible_nodes = self.return_visible_nodes()

        cycle_ids_of_nodes = self.return_node_data(attribute="cycle_id")

        core_index = visible_nodes.index(node)

        if central_node_in_cycle:

            central_cycle_id = cycle_ids_of_node[0]

            local_non_cycle_nodes = [
                i
                for i, indices in enumerate(cycle_ids_of_nodes)
                if i != core_index
                and (indices is None or central_cycle_id not in indices)
            ]

        else:

            local_non_cycle_nodes = []

        if (
            central_node_in_cycle
            and n_cycles == 1
            and cycle_size == 3
            and vsepr_class == "tetrahedral"
        ):

            final_angles = []

            global_non_cycle_nodes = [
                visible_nodes[i] for i in local_non_cycle_nodes
            ]

            global_cycle_nodes = [
                i
                for i in visible_nodes
                if i not in global_non_cycle_nodes and i != node
            ]

            for global_cycle_node in global_cycle_nodes:

                for global_non_cycle_node in global_non_cycle_nodes:

                    final_angles.append(
                        {
                            "nodes": [global_non_cycle_node, global_cycle_node],
                            "angle": 120.0 * np.pi / 180.0,
                        }
                    )

            final_angles.append(
                {"nodes": global_cycle_nodes, "angle": 60.0 * np.pi / 180.0}
            )

            final_angles.append(
                {"nodes": global_non_cycle_nodes, "angle": 109.0 * np.pi / 180.0}
            )

            return final_angles

        elif (
            central_node_in_cycle
            and n_cycles == 1
            and cycle_size == 3
            and vsepr_class in ["bent", "bent_tetrahedral"]
        ):

            return [{"nodes": adjacent_atoms, "angle": 60.0 * np.pi / 180.0}]

        elif (
            central_node_in_cycle
            and n_cycles == 1
            and vsepr_class == "trigonal_planar"
            and len(local_non_cycle_nodes) == 1
        ):

            final_angles = []

            if cycle_size == 3:

                major_angle = 150.0

                minor_angle = 60.0

            elif cycle_size == 4:

                major_angle = 135.0

                minor_angle = 90.0

            elif cycle_size == 5:

                major_angle = 126.0

                minor_angle = 108.0

            else:

                major_angle = 120.0

                minor_angle = 120.0

            global_non_cycle_node = [visible_nodes[i] for i in local_non_cycle_nodes][0]

            global_cycle_nodes = [
                i for i in visible_nodes if i != global_non_cycle_node and i != node
            ]

            for cycle_node in global_cycle_nodes:

                final_angles.append(
                    {
                        "nodes": [global_non_cycle_node, cycle_node],
                        "angle": major_angle * np.pi / 180.0,
                    }
                )

            final_angles.append(
                {"nodes": global_cycle_nodes, "angle": minor_angle * np.pi / 180.0}
            )

            return final_angles

        coordinates = self.return_node_data(attribute="coordinates")

        n_atoms = self.return_n_atoms()

        local_angle_nodes = [
            [i, core_index, j]
            for i in range(n_atoms)
            if i != core_index
            for j in range(n_atoms)
            if j != core_index and j > i
        ]

        global_angle_nodes = [
            [visible_nodes[angle_nodes[0]], visible_nodes[angle_nodes[2]]]
            for angle_nodes in local_angle_nodes
        ]

        real_angles = [
            Vector(coordinates[angle]).calculate_angle() for angle in local_angle_nodes
        ]

        ideal_angles = list(
            set(
                [
                    angle["angle"]
                    for angle in self.ideal_geometries[vsepr_class]["angles"]
                ]
            )
        )

        final_angles = []

        for angle_nodes, real_angle in zip(global_angle_nodes, real_angles):

            deviations = [
                np.abs(real_angle - ideal_angle) for ideal_angle in ideal_angles
            ]

            argmin_dev = np.argmin(deviations)

            final_angles.append(
                {
                    "nodes": angle_nodes,
                    "angle": ideal_angles[argmin_dev] * np.pi / 180.0,
                }
            )

        return final_angles

    def _constrain_node_dihedral(
        self, node, nodes, vsepr_class, is_adjacent_visible=False, include_all=False
    ):

        if vsepr_class not in [
            "tetrahedral",
            "trigonal_bent",
            # "seesaw",
            # "seesaw_rectangular",
            "trigonal_planar",
            "square_planar",
            # "formate",
        ]:

            return []

        self.set_nodes_visible(visible_nodes=nodes)

        adjacent_nodes = self.return_adjacent_nodes(
            node, is_adjacent_visible=is_adjacent_visible, include_all=include_all
        )

        self.set_nodes_visible(visible_nodes=[node] + adjacent_nodes)

        visible_nodes = self.return_visible_nodes()

        core_index = visible_nodes.index(node)

        if vsepr_class == "tetrahedral":
            # return self._constrain_pyramidal_dihedral(core_index, visible_nodes)
            return self._constrain_squezzed_tetrahedral_dihedral(
                core_index, adjacent_nodes, visible_nodes
            )

        elif vsepr_class == "trigonal_bent":
            return self._constrain_pyramidal_dihedral(core_index, visible_nodes)

        else:
            return self._constrain_planar_dihedral(core_index, visible_nodes)

    def _constrain_squezzed_tetrahedral_dihedral(
        self,
        core_node_index: int,
        adjacent_nodes: np.ndarray,
        visible_nodes: np.ndarray,
    ):
        """
        Constrain the dihedral angles of a node with tetrahedral geometry if node is part of two cycles and the tetrahedral geometry is far from ideal.
        Prevents flip of chirality.

        Returns:
            list: list of dictionaries containing the nodes and the ideal dihedral angle
        """
        cycle_ids = self.return_node_data(attribute="cycle_id")
        core_cycle_id = cycle_ids[core_node_index]
        # If tetrahedral node is not part of two cycles return no dihedral constraints as otherwise it is constrained enough
        if core_cycle_id is None or len(core_cycle_id) < 2:
            return []
        adjacent_cycle_ids = [cycle_ids[visible_nodes.index(i)] for i in adjacent_nodes]
        non_cycle_nodes = [
            adjacent_node
            for cycle_id, adjacent_node in zip(adjacent_cycle_ids, adjacent_nodes)
            if cycle_id is None or len(set(cycle_id).intersection(core_cycle_id)) == 0
        ]
        n_non_cycle_nodes = len(non_cycle_nodes)
        ## If no adjacent is no part of a cycle return no dihedral constraints as otherwise it is constrained enough.
        if n_non_cycle_nodes == 0:
            cycle_size = self.return_node_data(attribute="smallest_cycle_size")[
                core_node_index
            ]
            if cycle_size == 3:
                return []
            else:
                return self._constrain_pyramidal_dihedral(
                    core_node_index, visible_nodes
                )
        elif n_non_cycle_nodes > 1:
            return []
        non_cycle_node = non_cycle_nodes[0]
        non_cycle_node_index = visible_nodes.index(non_cycle_node)
        cycle_node_indices = [
            visible_nodes.index(i) for i in adjacent_nodes if i != non_cycle_node
        ]
        core_node = visible_nodes[core_node_index]
        self.nodes[core_node]["vsepr_classes"] = "seesaw_rectangular"
        core_node_index = visible_nodes.index(core_node)
        cycle_ids = self.return_node_data(attribute="cycle_id")
        rect_seesaw = False
        try:
            orthogonal_node_index = [
                cycle_node_index
                for cycle_node_index in cycle_node_indices
                if set(core_cycle_id).issubset(set(cycle_ids[cycle_node_index]))
            ][0]
            linear_dependent_nodes = [
                node
                for node in cycle_node_indices
                if (node != orthogonal_node_index) and (node != core_node_index)
            ]
            angle_data = [
                {
                    "nodes": [
                        visible_nodes[linear_dependent_nodes[0]],
                        visible_nodes[linear_dependent_nodes[1]],
                    ],
                    "angle": np.pi,
                },
                {
                    "nodes": [
                        visible_nodes[linear_dependent_nodes[0]],
                        visible_nodes[orthogonal_node_index],
                    ],
                    "angle": np.pi / 2.0,
                },
                {
                    "nodes": [
                        visible_nodes[linear_dependent_nodes[1]],
                        visible_nodes[orthogonal_node_index],
                    ],
                    "angle": np.pi / 2.0,
                },
                {
                    "nodes": [visible_nodes[linear_dependent_nodes[0]], non_cycle_node],
                    "angle": np.pi / 2.0,
                },
                {
                    "nodes": [visible_nodes[linear_dependent_nodes[1]], non_cycle_node],
                    "angle": np.pi / 2.0,
                },
                {
                    "nodes": [visible_nodes[orthogonal_node_index], non_cycle_node],
                    "angle": np.pi / 2.0,
                },
            ]
            self.nodes[core_node]["angle_constraints"] = angle_data

            dihedral_nodes = [
                [
                    linear_dependent_nodes[0],
                    core_node_index,
                    orthogonal_node_index,
                    linear_dependent_nodes[1],
                ],
                [
                    linear_dependent_nodes[1],
                    core_node_index,
                    orthogonal_node_index,
                    linear_dependent_nodes[0],
                ],
                [
                    orthogonal_node_index,
                    core_node_index,
                    linear_dependent_nodes[0],
                    visible_nodes.index(non_cycle_nodes[0]),
                ],
                [
                    orthogonal_node_index,
                    core_node_index,
                    linear_dependent_nodes[1],
                    visible_nodes.index(non_cycle_nodes[0]),
                ],
            ]
            rect_seesaw = True

        except IndexError:
            dihedral_nodes = [
                [
                    cycle_node_indices[0],
                    core_node_index,
                    cycle_node_indices[1],
                    non_cycle_node_index,
                ],
                [
                    cycle_node_indices[0],
                    core_node_index,
                    cycle_node_indices[2],
                    non_cycle_node_index,
                ],
                [
                    cycle_node_indices[1],
                    core_node_index,
                    cycle_node_indices[0],
                    non_cycle_node_index,
                ],
                [
                    cycle_node_indices[1],
                    core_node_index,
                    cycle_node_indices[2],
                    non_cycle_node_index,
                ],
                [
                    cycle_node_indices[2],
                    core_node_index,
                    cycle_node_indices[0],
                    non_cycle_node_index,
                ],
                [
                    cycle_node_indices[2],
                    core_node_index,
                    cycle_node_indices[1],
                    non_cycle_node_index,
                ],
            ]

        coordinates = self.return_node_data(attribute="coordinates")
        real_dihedrals = [
            Vector(coordinates[dihedral]).calculate_dihedral()
            for dihedral in dihedral_nodes
        ]
        assigned_dihedrals = []
        for i, real_dihedral in enumerate(real_dihedrals):
            deviations = []
            if i < 2 and rect_seesaw:
                ideal_dihedrals = [180.0, -180.0]
            elif i >= 2 and rect_seesaw:
                ideal_dihedrals = [90.0, -90.0]
            else:
                ideal_dihedrals = [120.0, -120.0]
            for ideal_dihedral in ideal_dihedrals:
                deviation = np.abs(ideal_dihedral - real_dihedral)
                if deviation < 180.0:
                    deviations.append(deviation)
                else:
                    deviations.append(np.inf)
            smallest_deviation = np.argmin(deviations)
            if deviations[smallest_deviation] != np.inf:
                assigned_dihedrals.append(ideal_dihedrals[smallest_deviation])
            else:
                assigned_dihedrals.append(np.inf)

        final_dihedrals = [
            {
                "nodes": [
                    visible_nodes[dihedral[0]],
                    visible_nodes[dihedral[2]],
                    visible_nodes[dihedral[3]],
                ],
                "angle": assigned_angle * np.pi / 180.0,
            }
            for dihedral, assigned_angle in zip(dihedral_nodes, assigned_dihedrals)
            if assigned_angle != np.inf
        ]

        return final_dihedrals

    def _constrain_pyramidal_dihedral(self, core_index, visible_nodes):
        cycle_ids = self.return_node_data(attribute="cycle_id")
        core_cycle_id = cycle_ids[core_index]
        if core_cycle_id is None:
            return []
        cycle_size = self.return_node_data(attribute="smallest_cycle_size")[core_index]
        if cycle_size == 3:
            cycle_ids = self.return_node_data(attribute="cycle_id")
            core_cycle_id = cycle_ids[core_index]
            core_node = visible_nodes[core_index]
            adjacent_nodes = self.return_adjacent_nodes(core_node)
            adjacent_cycle_ids = [
                cycle_ids[visible_nodes.index(i)] for i in adjacent_nodes
            ]
            non_cycle_nodes = [
                adjacent_node
                for cycle_id, adjacent_node in zip(adjacent_cycle_ids, adjacent_nodes)
                if cycle_id is None
                or len(set(cycle_id).intersection(core_cycle_id)) == 0
            ]
            if len(non_cycle_nodes) == 1:
                non_cycle_node = non_cycle_nodes[0]
                cycle_nodes = [
                    node for node in adjacent_nodes if node != non_cycle_node
                ]
                cycle_node_indices = [visible_nodes.index(node) for node in cycle_nodes]
                non_cycle_node_index = visible_nodes.index(non_cycle_node)
                angle_data = [
                    {
                        "nodes": [cycle_nodes[0], non_cycle_node],
                        "angle": 110.0 * np.pi / 180.0,
                    },
                    {
                        "nodes": [cycle_nodes[1], non_cycle_node],
                        "angle": 110.0 * np.pi / 180.0,
                    },
                ]
                self.nodes[core_node]["angle_constraints"] = angle_data
                self.nodes[core_node]["vsepr_classes"] = "trigonal_bent_squeezed"
                dihedral_nodes = [
                    [
                        cycle_node_indices[0],
                        core_index,
                        cycle_node_indices[1],
                        non_cycle_node_index,
                    ],
                    [
                        cycle_node_indices[1],
                        core_index,
                        cycle_node_indices[0],
                        non_cycle_node_index,
                    ],
                ]
                coordinates = self.return_node_data(attribute="coordinates")
                real_dihedrals = [
                    Vector(coordinates[dihedral]).calculate_dihedral()
                    for dihedral in dihedral_nodes
                ]
                assigned_dihedrals = []
                ideal_dihedrals = [110.0, -110.0]
                for real_dihedral in real_dihedrals:
                    deviations = [
                        np.abs(ideal_dihedral - real_dihedral)
                        for ideal_dihedral in ideal_dihedrals
                    ]
                    smallest_deviation = np.argmin(deviations)
                    assigned_dihedrals.append(ideal_dihedrals[smallest_deviation])

                final_dihedrals = [
                    {
                        "nodes": [
                            visible_nodes[dihedral[0]],
                            visible_nodes[dihedral[2]],
                            visible_nodes[dihedral[3]],
                        ],
                        "angle": assigned_angle * np.pi / 180.0,
                    }
                    for dihedral, assigned_angle in zip(
                        dihedral_nodes, assigned_dihedrals
                    )
                ]
                return final_dihedrals

        ideal_dihedrals = [120.0, -120.0]

        n_atoms = self.return_n_atoms()

        coordinates = self.return_node_data(attribute="coordinates")

        dihedral_nodes = self._get_dihedral_nodes(coordinates, n_atoms, core_index)

        real_dihedrals = [
            Vector(coordinates[dihedral]).calculate_dihedral()
            for dihedral in dihedral_nodes
        ]

        assigned_dihedrals = []

        for real_dihedral in real_dihedrals:

            deviations = []

            for ideal_dihedral in ideal_dihedrals:

                deviation = np.abs(ideal_dihedral - real_dihedral)

                if deviation < 90.0:

                    deviations.append(deviation)

                else:

                    deviations.append(np.inf)

            smallest_deviation = np.argmin(deviations)

            if deviations[smallest_deviation] != np.inf:

                assigned_dihedrals.append(ideal_dihedrals[smallest_deviation])

            else:

                assigned_dihedrals.append(np.inf)

        return [
            {
                "nodes": [
                    visible_nodes[dihedral[0]],
                    visible_nodes[dihedral[2]],
                    visible_nodes[dihedral[3]],
                ],
                "angle": assigned_angle * np.pi / 180.0,
            }
            for dihedral, assigned_angle in zip(dihedral_nodes, assigned_dihedrals)
            if assigned_angle != np.inf
        ]

    def _constrain_planar_dihedral(self, core_index, visible_nodes):

        ideal_dihedrals = [180.0, -180.0]

        n_atoms = self.return_n_atoms()

        coordinates = self.return_node_data(attribute="coordinates")

        dihedral_nodes = self._get_dihedral_nodes(coordinates, n_atoms, core_index)

        real_dihedrals = [
            Vector(coordinates[dihedral]).calculate_dihedral()
            for dihedral in dihedral_nodes
        ]

        assigned_dihedrals = []

        for real_dihedral in real_dihedrals:

            deviations = []

            for ideal_dihedral in ideal_dihedrals:

                deviation = np.abs(ideal_dihedral - real_dihedral)

                # Dummy value, alway flatten, no matter the deviation
                if deviation < 180.0:

                    deviations.append(deviation)

                else:

                    deviations.append(np.inf)

            smallest_deviation = np.argmin(deviations)

            if deviations[smallest_deviation] != np.inf:

                assigned_dihedrals.append(ideal_dihedrals[smallest_deviation])

            else:

                assigned_dihedrals.append(np.inf)

        return [
            {
                "nodes": [
                    visible_nodes[dihedral[0]],
                    visible_nodes[dihedral[2]],
                    visible_nodes[dihedral[3]],
                ],
                "angle": assigned_angle * np.pi / 180.0,
            }
            for dihedral, assigned_angle in zip(dihedral_nodes, assigned_dihedrals)
            if assigned_angle != np.inf
        ]

    def _get_dihedral_nodes(self, coordinates, n_atoms, core_index):

        # dihedral_nodes =  [[i, core_index, j, k] for i in range(n_atoms) for j in range(n_atoms) for k in range(n_atoms) if len(set([i, core_index, j, k])) == len([i, core_index, j, k]) and k > i and i != core_index and j != core_index and k != core_index]

        dihedral_nodes = [
            [i, core_index, j, k]
            for i in range(n_atoms)
            for j in range(n_atoms)
            for k in range(n_atoms)
            if len(set([i, core_index, j, k])) == len([i, core_index, j, k])
        ]

        final_dihedral_nodes = []

        for dihedral in dihedral_nodes:

            angle1 = Vector(
                coordinates[[dihedral[0], dihedral[1], dihedral[2]]]
            ).calculate_angle()

            angle2 = Vector(
                coordinates[[dihedral[2], dihedral[1], dihedral[3]]]
            ).calculate_angle()

            if (angle1 < 150.0 or angle1 > 210.0) and (
                angle2 < 150.0 or angle2 > 210.0
            ):

                final_dihedral_nodes.append(dihedral)

        return final_dihedral_nodes
