import yaml
import numpy as np
from molbar.exceptions.error import NotInputFormat


def _get_constraints(input_file: str) -> dict:
    """
    Returns the constraints from the input file.

    Args:
        input_file (str): file path of input file.

    Returns:

        dict: constraints from input file.
    """

    with open(input_file, "r") as file:
        input = yaml.safe_load(file)

    keys = list(input.keys())

    possible_keys = [
        "bond_order_assignment",
        "cycle_detection",
        "set_edges",
        "set_angles",
        "set_dihedrals",
        "set_repulsion",
        "repulsion_charge",
        "constraints",
    ]

    if all([True if key in possible_keys else False for key in keys]) == False:

        raise NotInputFormat(input_file)

    return input


def _transform_constraints(file_path: str, input: dict) -> dict:
    """
    Transforms the constraints from the input file to format that is used for MolGraph.

    Args:
        file_path (str): file path of input file.
        input (dict): constraints in MolBar format.

    Returns:

        dict: constraints in MolGraph format.
    """

    if input and input.get("constraints") != None:

        bond_constraints = input.get("constraints").get("bonds")

        if bond_constraints != None:

            input["constraints"]["bonds"] = _transform_bond_constraints(
                file_path, bond_constraints
            )

        angle_constraints = input.get("constraints").get("angles")

        if angle_constraints != None:

            input["constraints"]["angles"] = _transform_angle_constraints(
                file_path, angle_constraints
            )

        dihedral_constraints = input.get("constraints").get("dihedrals")

        if dihedral_constraints != None:

            input["constraints"]["dihedrals"] = _transform_dihedral_constraints(
                file_path, dihedral_constraints
            )

    if not input:

        input = {}

    return input


def _transform_bond_constraints(file, bond_constraints: list) -> dict:
    """
    Transforms the angle constraints from the input file to format that is used for MolGraph.
    Args:
        bond_constraints (list): List of bond constraints in input file

    Raises:
        NotYMLFormat: If input file is not in MolBar input YML format.

    Returns:
        dict: bond constraints in MolGraph format. {node1: [{"nodes": [node1], "bond_length": length}, ...], node2:... ]
    """

    transform_bond_constraints = {}

    for bond in bond_constraints:

        try:

            bond_atoms = sorted(bond["atoms"])

            new_constraint = [{"nodes": bond["atoms"][1], "bond_length": bond["value"]}]

            if bond_atoms[0] in transform_bond_constraints.keys():

                transform_bond_constraints[bond_atoms[0]] += new_constraint

            else:

                transform_bond_constraints[bond_atoms[0]] = new_constraint

        except KeyError:

            raise NotInputFormat(file)

    return transform_bond_constraints


def _transform_angle_constraints(file, angle_constraints: list) -> dict:
    """
    Transforms the angle constraints from the input file to format that is used for MolGraph.
    Args:
        angle_constraints (list): List of angle constraints in input file

    Raises:
        NotYMLFormat: If input file is not in MolBar input YML format.

    Returns:
        dict: Angle constraints in MolGraph format. {node1: [{"nodes": [node1, node2], "angle": angle}, ...], node2:... ]
    """

    transform_angle_constraints = {}

    for angle in angle_constraints:

        try:

            central_node = angle["atoms"][1] - 1

            adjacent_nodes = [angle["atoms"][0] - 1, angle["atoms"][2] - 1]

            if adjacent_nodes[0] > adjacent_nodes[1]:

                adjacent_nodes = [adjacent_nodes[1], adjacent_nodes[0]]

            new_constraint = [
                {"nodes": adjacent_nodes, "angle": angle["value"] * np.pi / 180.0}
            ]

            if central_node in transform_angle_constraints.keys():

                transform_angle_constraints[central_node] += new_constraint

            else:

                transform_angle_constraints[central_node] = new_constraint

        except KeyError:

            raise NotInputFormat(file)

    return transform_angle_constraints


def _transform_dihedral_constraints(file, dihedral_constraints: list) -> dict:
    """
    Transforms the dihedral constraints from the input file to format that is used for MolGraph.
    Args:
        file (str): File path of input file.
        dihedral_constraints (list): List of angle constraints in input file

    Raises:
        NotYMLFormat: If input file is not in MolBar input YML format.

    Returns:
        dict: Dihedral constraints in MolGraph format. {node1: [{"nodes": [node1, node2, node3], "dihedral": angle}, ...], node2:... ]
    """

    transform_dihedral_constraints = {}

    for dihedral in dihedral_constraints:

        try:

            dihedral_nodes = [
                dihedral["atoms"][0] - 1,
                dihedral["atoms"][1] - 1,
                dihedral["atoms"][2] - 1,
                dihedral["atoms"][3] - 1,
            ]

            if dihedral_nodes[1] > dihedral_nodes[2]:

                dihedral_nodes = [
                    dihedral_nodes[3],
                    dihedral_nodes[2],
                    dihedral_nodes[1],
                    dihedral_nodes[0],
                ]

            central_node = dihedral_nodes[1]

            adjacent_nodes = [dihedral_nodes[0], dihedral_nodes[2], dihedral_nodes[3]]

            new_constraint = [
                {"nodes": adjacent_nodes, "dihedral": dihedral["value"] * np.pi / 180.0}
            ]

            if central_node in transform_dihedral_constraints.keys():

                transform_dihedral_constraints[central_node] += new_constraint

            else:

                transform_dihedral_constraints[central_node] = new_constraint

        except KeyError:

            raise NotInputFormat(file)

    return transform_dihedral_constraints
