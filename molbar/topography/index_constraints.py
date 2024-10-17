def _get_angle_constraints(
    fragment: list, all_angle_constraints: list, index_conversion: dict
) -> tuple:
    """
    Converts angle constraints with global node indices to local node indices.

    Args:
        fragment (list): list of global node indices of one fragment
        all_angle_constraints (list): list of all angle constraints for each node in the fragment
        index_conversion (dict): dictionary that maps global node indices to local node indices

    Returns:
        tuple: list of angle constraints with local node indices, list of ideal angles
    """
    angles = []

    ideal_angles = []

    # Iterate over all nodes in the fragment and get angle constraints for each node
    for node, angle_constraints in zip(fragment, all_angle_constraints):

        # If no angle constraints for the node, continue
        if angle_constraints is None:

            continue

        # Iterate over all angle constraints for the node
        for angle in angle_constraints:
            # Reindex the angle constraints to local node indices
            angles.append(
                [
                    index_conversion[angle["nodes"][0]],
                    index_conversion[node],
                    index_conversion[angle["nodes"][1]],
                ]
            )

            ideal_angles.append(angle["angle"])

    return angles, ideal_angles


def _get_dihedral_constraints(
    visible_edges: list,
    visible_nodes: list,
    all_edge_dihedral_constraints: list,
    all_node_dihedral_constraints: list,
    index_conversion: dict,
) -> tuple:
    """
    Converts dihedral constraints with global node indices to local node indices.

    Args:
        visible_edges (list): list of visible edges, each edge is a tuple of two global node indices
        visible_nodes (list): list of visible nodes, each node is a global node index
        all_edge_dihedral_constraints (list): list of all dihedral constraints for each edge in the fragment
        all_node_dihedral_constraints (list): list of all dihedral constraints for each node in the fragment
        index_conversion (dict): dictionary that maps global node indices to local node indices

    Returns:
        tuple: list of dihedral constraints with local node indices, list of ideal dihedrals
    """

    dihedrals = []

    ideal_dihedrals = []

    # Get dihedral constraints for edges such as for double bond and cycles
    dihedrals, ideal_dihedrals = _get_dihedral_edge_constraints(
        visible_edges,
        all_edge_dihedral_constraints,
        index_conversion,
        dihedrals,
        ideal_dihedrals,
    )

    # Get dihedral constraints for nodes to constraint trigonal planar atoms
    dihedrals, ideal_dihedrals = _get_dihedral_node_constraints(
        visible_nodes,
        all_node_dihedral_constraints,
        index_conversion,
        dihedrals,
        ideal_dihedrals,
    )

    return dihedrals, ideal_dihedrals


def _get_dihedral_edge_constraints(
    visible_edges: list,
    all_edge_dihedral_constraints: list,
    index_conversion: dict,
    dihedrals: list,
    ideal_dihedrals: list,
) -> tuple:
    """
    Converts dihedral edge constraints with global node indices to local node indices.

    Args:
        visible_edges (list): list of visible edges, each edge is a tuple of two global node indices
        all_edge_dihedral_constraints (list): list of all dihedral constraints for each edge in the fragment
        index_conversion (dict): dictionary that maps global node indices to local node indices
        dihedrals (list): list of four local node indices for each dihedral constraint
        ideal_dihedrals (list): list of ideal dihedral values for each dihedral constraint

    Returns:
        tuple: list of dihedral constraints with local node indices, list of ideal dihedrals
    """

    # Iterate over all visible edges
    for edge, dihedral_constraints in zip(visible_edges, all_edge_dihedral_constraints):

        # If no dihedral constraints for the edge, continue
        if dihedral_constraints is None:

            continue
        # Iterate over all dihedral constraints for the edge
        for dihedral in dihedral_constraints:

            # Reindex the dihedral constraints to local node indices, but consider ascending order of the edge indices
            if edge[0] < edge[1]:

                new_dihedral = [
                    index_conversion[dihedral["nodes"][0]],
                    index_conversion[edge[0]],
                    index_conversion[edge[1]],
                    index_conversion[dihedral["nodes"][1]],
                ]

            else:

                new_dihedral = [
                    index_conversion[dihedral["nodes"][0]],
                    index_conversion[edge[1]],
                    index_conversion[edge[0]],
                    index_conversion[dihedral["nodes"][1]],
                ]

            dihedrals.append(new_dihedral)

            ideal_dihedrals.append(dihedral["dihedral"])

    return dihedrals, ideal_dihedrals


def _get_dihedral_node_constraints(
    visible_nodes: list,
    all_node_dihedral_constraints: list,
    index_conversion: dict,
    dihedrals: list,
    ideal_dihedrals: list,
) -> tuple:
    """
    Converts dihedral node constraints with global node indices to local node indices.

    Args:
        visible_nodes (list): list of visible nodes, each node is a global node index
        all_node_dihedral_constraints (list): list of all dihedral constraints for each node in the fragment
        index_conversion (dict): dictionary that maps global node indices to local node indices
        dihedrals (list): list of four local node indices for each dihedral constraint
        ideal_dihedrals (list): list of ideal dihedral values for each dihedral constraint

    Returns:
        tuple: list of dihedral constraints with local node indices, list of ideal dihedrals
    """

    # Iterate over all visible nodes
    for node, dihedral_constraints in zip(visible_nodes, all_node_dihedral_constraints):
        # If no dihedral constraints for the node, continue
        if dihedral_constraints is None:

            continue
        # Iterate over all dihedral constraints for the node
        for dihedral in dihedral_constraints:
            # Reindex the dihedral constraints to local node indices
            dihedrals.append(
                [
                    index_conversion[dihedral["nodes"][0]],
                    index_conversion[node],
                    index_conversion[dihedral["nodes"][1]],
                    index_conversion[dihedral["nodes"][2]],
                ]
            )

            ideal_dihedrals.append(dihedral["angle"])

    return dihedrals, ideal_dihedrals
