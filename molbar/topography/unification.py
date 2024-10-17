import numpy as np
import os
from molbar.helper.printer import Printer
from molbar.helper.debug import _mkdir_debug_path
from molbar.molgraph.molgraph import MolGraph
from molbar.topography.geomopt import GeomOpt
from molbar.topography.index_constraints import (
    _get_angle_constraints,
    _get_dihedral_constraints,
)
from scipy.spatial.distance import cdist


def unify_fragments(molgraph: MolGraph, write_trj=False, debug_path="debug"):
    """
    Unifies the coordinates of all fragments of a molecule.

    Args:
        molgraph (MolGraph): MolGraph object of the molecule.
        write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
        debug_path (str, optional): Path to write the trajectory to. Defaults to "debug".

    Returns:
        MolGraph: MolGraph object of the molecule with unified coordinates.
    """

    # Get fragment ids of all nodes
    fragment_ids = molgraph.return_node_data(attribute="fragment_id")

    # Get visible nodes
    visible_nodes = molgraph.return_visible_nodes()

    visible_edges = molgraph.return_visible_edges()

    # Get number of atoms
    n_atoms = molgraph.return_n_atoms()

    # Cluster nodes by fragment id
    fragment_nodes = {
        fragment_id: [
            node for node in visible_nodes if fragment_ids[node] == fragment_id
        ]
        for fragment_id in set(fragment_ids)
        if fragment_id is not None
    }

    fragments_data = {}

    if len(fragment_nodes.keys()) == 0:

        molgraph.unified_coulomb_matrix = _refine_cm_matrix(
            np.zeros((n_atoms, n_atoms)), visible_edges, molgraph
        )

        molgraph.fragments_data = fragments_data

        return molgraph

    # Initialize empty cm matrix
    cm_matrix = np.empty((n_atoms, n_atoms, len(fragment_nodes.keys())))

    # Fill cm matrix with nan values
    cm_matrix[:] = np.nan

    # If trajectory should be written, create directory
    if write_trj:

        _mkdir_debug_path(debug_path)

    fragment_ids = sorted(list(fragment_nodes.keys()))

    # Iterate over all fragments
    for fragment_id in fragment_ids:

        fragment = fragment_nodes[fragment_id]

        # Unify fragment based on constraints
        fragments_data, fragment_cm_matrix, visible_nodes = _unify_fragment(
            molgraph,
            fragments_data,
            fragment_id,
            fragment,
            write_trj=write_trj,
            debug_path=debug_path,
        )

        ith_nodes, jth_nodes = np.meshgrid(visible_nodes, visible_nodes)

        # Store cm matrix of fragment in total cm matrix in the corresponding fragment slice
        cm_matrix[ith_nodes, jth_nodes, fragment_id] = fragment_cm_matrix

    # Count the number how often one cm matrix element is in a fragment, that is, how often it is not nan
    count_element_in_fragments = np.count_nonzero(~np.isnan(cm_matrix), axis=2)

    # Set all number of elements that are 0 to 1 to avoid division by zero
    count_element_in_fragments[np.where(count_element_in_fragments == 0)] = 1

    # Set all nan values to 0 to avoid nan values in the average cm matrix
    cm_matrix[np.isnan(cm_matrix)] = 0

    # Average cm matrix as some elements might be assigned twice if two atoms are part of several fragments

    average_cm_matrix = np.sum(cm_matrix, axis=2) / count_element_in_fragments

    refined_cm_matrix = _refine_cm_matrix(average_cm_matrix, visible_edges, molgraph)

    molgraph.unified_coulomb_matrix = refined_cm_matrix

    molgraph.fragments_data = fragments_data

    return molgraph


def _refine_cm_matrix(
    cm_matrix: np.ndarray, visible_edges: np.ndarray, molgraph: MolGraph
) -> np.ndarray:
    """
    Refines the Coulomb matrix by setting 1/rij to 1/dij for adjacent atoms as slight deviations can occur during the unification process.

    Args:

        cm_matrix (np.ndarray): Coulomb matrix of the molecule.
        visible_edges (np.ndarray): Array of visible edges.
        molgraph (MolGraph): MolGraph object of the molecule.

    Returns:

        np.ndarray: Refined Coulomb matrix.
    """

    for edge in visible_edges:

        try:

            cm_matrix[edge[0], edge[1]] = 1 / molgraph.edges[edge]["dij"]

            cm_matrix[edge[1], edge[0]] = 1 / molgraph.edges[edge]["dij"]

        except KeyError:

            continue

    return cm_matrix


def _unify_fragment(
    molgraph: MolGraph,
    fragments_data: dict,
    fragment_id: int,
    fragment: list,
    write_trj=False,
    debug_path=None,
) -> tuple:
    """
    Unifies the coordinates of one fragment using a force field based on constraints.

    Args:
        molgraph (MolGraph): MolGraph object of the molecule.
        fragments_data (dict): dict containing the data of all fragments.
        fragment_id (int): fragment id of the fragment to unify.
        fragment (list): list of nodes of the fragment to unify.
        write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
        debug_path (str, optional): Path to write the trajectory to. Defaults to "debug".

    Returns:
        tuple: tuple containing: fragments_data (dict): dict containing the data of all fragments, coulomb_matrix (np.array): Coulomb matrix of the fragment, visible_nodes (list): list of visible nodes of the fragment.
    """
    # Reindex constraints for fragment
    (
        core_nodes,
        core_indices,
        bonds,
        ideal_bonds,
        angles,
        ideal_angles,
        dihedrals,
        ideal_dihedrals,
        repulsion,
        charges,
        visible_nodes,
    ) = _prepare_fragment_for_unification(molgraph, fragment)

    # Initialize GeomOpt object
    geomopt = GeomOpt(
        n_atoms=molgraph.return_n_atoms(),
        n_core=len(core_nodes),
        core_indices=core_indices,
        coordinates=molgraph.return_node_data(attribute="coordinates"),
        atoms=molgraph.return_node_data(attribute="elements"),
        bonds=bonds,
        ideal_bonds=ideal_bonds,
        angles=angles,
        ideal_angles=ideal_angles,
        dihedrals=dihedrals,
        ideal_dihedrals=ideal_dihedrals,
        repulsion=repulsion,
        charges=charges,
        write_trj=write_trj,
        get_index=True,
        fragment_id=fragment_id,
        debug_path=debug_path,
    )

    # Unify fragment
    geomopt.optimize()

    # Adjust priorities for pseudochirality, fragments have been sorted before so that possible chiral fragments
    # are treated first, that their priorities can be adjusted to account for pseudo chirality for a later otherwise
    # non-chiral fragment.
    # molgraph._adjust_priorities_for_pseudochirality(core_indices, geomopt.absconf_index)

    # Get Coulomb matrix
    rij = cdist(geomopt.optimized_geometry, geomopt.optimized_geometry)

    np.fill_diagonal(rij, 1)

    coulomb_matrix = 1 / rij

    np.fill_diagonal(coulomb_matrix, 0)

    topography_matrix = coulomb_matrix.copy()

    np.fill_diagonal(
        topography_matrix, molgraph.return_node_data(attribute="atomic_numbers")
    )

    # Store fragment data
    fragment_data = {
        "final_energy": geomopt.energy,
        "convergence": geomopt.convergence_type,
    }

    fragment_data["visible_nodes"] = visible_nodes

    fragment_data["nodes"] = [str(node + 1) for node in fragment]

    fragment_data["core_indices"] = core_indices

    fragment_data["final_geometry"] = geomopt.optimized_geometry

    fragment_data["coulomb_matrix"] = coulomb_matrix

    if write_trj and debug_path is not None:

        fragment_spectrum, _ = np.linalg.eigh(topography_matrix)
        fragment_data["topography"] = list(fragment_spectrum)

        Printer(
            molgraph.return_n_atoms(),
            geomopt.energy,
            geomopt.optimized_geometry,
            molgraph.return_node_data(attribute="elements"),
            os.path.join(debug_path, str(fragment_id + 1) + ".xyz"),
        ).print()

    fragment_data["elements"] = molgraph.return_node_data(attribute="elements")

    fragments_data[fragment_id] = fragment_data

    return fragments_data, coulomb_matrix, visible_nodes


def _prepare_fragment_for_unification(molgraph: MolGraph, fragment: list) -> tuple:
    """
    Prepares a fragment for unification by reindexing constraints.

    Args:
        molgraph (MolGraph): MolGraph object of the molecule.
        fragment (list): list of nodes of the fragment to unify.

    Returns:
        tuple: tuple containing: core_nodes (list): list of core nodes of the fragment, core_indices (list): list of core indices of the fragment, bonds (list): list of bonds of the fragment, ideal_bonds (list): list of ideal bonds of the fragment, angles (list): list of angles of the fragment, ideal_angles (list): list of ideal angles of the fragment, dihedrals (list): list of dihedrals of the fragment, ideal_dihedrals (list): list of ideal dihedrals of the fragment, repulsion (list): list of repulsion constraints of the fragment, charges (list): list of charges of the fragment, visible_nodes (list): list of visible nodes of the fragment.
    """

    molgraph.set_nodes_visible(include_all=True)

    # Get adjacent nodes of fragment
    adjacent_nodes = molgraph.return_adjacent_nodes(core_nodes=fragment)

    molgraph.set_nodes_visible(visible_nodes=fragment)

    # Get core nodes of fragment
    core_nodes = molgraph.return_visible_nodes()

    # Get all angle constraints for core nodes
    all_angle_constraints = molgraph.return_node_data(attribute="angle_constraints")

    # Get all dihedral constraints for core nodes
    all_edge_dihedral_constraints = molgraph.return_edge_data(
        attribute="dihedral_constraints", nodes_visible=True
    )

    all_node_dihedral_constraints = molgraph.return_node_data(
        attribute="dihedral_constraints"
    )

    visible_edges = molgraph.return_visible_edges(nodes_visible=True)

    molgraph.set_nodes_visible(visible_nodes=fragment + adjacent_nodes)

    visible_nodes = molgraph.return_visible_nodes()

    # Generate index conversion dict (global index -> local index)
    index_conversion = {node: index for index, node in enumerate(visible_nodes)}

    # Reindex angle constraints
    angles, ideal_angles = _get_angle_constraints(
        fragment, all_angle_constraints, index_conversion
    )

    # Reindex dihedral constraints
    dihedrals, ideal_dihedrals = _get_dihedral_constraints(
        visible_edges,
        core_nodes,
        all_edge_dihedral_constraints,
        all_node_dihedral_constraints,
        index_conversion,
    )

    # Reindex bond constraints
    bonds = [
        [index_conversion[bond[0]], index_conversion[bond[1]]]
        for bond in molgraph.return_visible_edges(nodes_visible=True)
    ]

    # Get ideal bond lengths
    ideal_bonds = molgraph.return_edge_data(attribute="dij", nodes_visible=True)

    # Filter out edges that are not constrained, C=C=C structures, outside C atoms are defined as edge
    bonds = [bond for bond, dij in zip(bonds, ideal_bonds) if dij]

    ideal_bonds = [dij for dij in ideal_bonds if dij]

    visible_nodes = molgraph.return_visible_nodes()

    repulsion, charges = molgraph.return_repulsion_pairs()

    # Reindex bond constraints
    repulsion = [
        [index_conversion[combination[0]], index_conversion[combination[1]]]
        for combination in repulsion
    ]

    core_nodes = [
        core_node
        for core_node in core_nodes
        if len(molgraph.return_adjacent_nodes(core_nodes=core_node)) > 1
    ]

    core_indices = [visible_nodes.index(node) for node in core_nodes]

    return (
        core_nodes,
        core_indices,
        bonds,
        ideal_bonds,
        angles,
        ideal_angles,
        dihedrals,
        ideal_dihedrals,
        repulsion,
        charges,
        visible_nodes,
    )
