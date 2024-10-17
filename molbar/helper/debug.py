import json
import os
import shutil
from molbar.molgraph.molgraph import MolGraph
from molbar.barcodes.barcodes import _calculate_barcode
from molbar.topography.index_constraints import (
    _get_angle_constraints,
    _get_dihedral_constraints,
)


def _mkdir_debug_path(debug_path: str) -> None:
    """
    Creates directory for debug files.

    Args:
        debug_path (str): debug directory path.

    Returns:
        None
    """

    if os.path.isdir(debug_path):

        shutil.rmtree(debug_path, ignore_errors=True)

    os.mkdir(debug_path)


def _get_debug_path(filepath: str) -> str:
    """
    Returns the file path of the debug file.

    Args:
        filepath (str): file path of input file.

    Returns:
        str: file path of debug file.
    """
    if not filepath:

        return False

    filename_with_extension = os.path.basename(filepath)

    parent_directory = os.path.dirname(filepath)

    filename, extension = os.path.splitext(filename_with_extension)

    debug_path = os.path.join(
        parent_directory, "data_molbar_" + filename + "_" + extension.replace(".", "")
    )

    return debug_path


def _prepare_debug_data(molbar: str, molgraph: MolGraph) -> None:
    """
    Writes debug.json file containing all information about the molecule.
    """
    data = {}

    molgraph.set_nodes_visible(include_all=True)

    data["MolBar"] = molbar

    try:
        data["topology_spectrum"] = list(
            [float(eps) for eps in molgraph.topology_spectrum]
        )
    except AttributeError:
        data["topology_spectrum"] = []

    try:
        data["heavy_atom_topology_spectrum"] = list(
            [float(eps) for eps in molgraph.heavy_atom_topology_spectrum]
        )
    except AttributeError:
        data["heavy_atom_topology_spectrum"] = []

    try:
        data["topography_spectrum"] = list(
            [float(eps) for eps in molgraph.topography_spectrum]
        )
    except AttributeError:
        data["topography_spectrum"] = []

    try:
        data["absolute_configuration_spectrum"] = list(
            [float(eps) for eps in molgraph.absconfiguration_spectrum]
        )
    except AttributeError:
        data["absolute_configuration_spectrum"] = []

    try:
        data["unified_coulomb_matrix"] = [
            [float(x) for x in row] for row in molgraph.unified_coulomb_matrix
        ]
    except AttributeError:
        data["unified_coulomb_matrix"] = []

    try:

        data["final_energies"] = [
            data["final_energy"] for key, data in molgraph.fragments_data.items()
        ]

    except AttributeError:

        data["final_energies"] = ["no fragments found."]

    timings = {}

    try:
        timings["total_timing"] = molgraph.total_time
    except AttributeError:
        timings["total_timing"] = 0.0
    try:
        timings["bond_order_assignment_timing"] = molgraph.bo_time
    except AttributeError:
        timings["bond_order_timing"] = 0.0

    try:
        timings["prioritisation_timing"] = molgraph.prio_time
    except AttributeError:
        timings["prioritisation_timing"] = 0.0

    try:
        timings["cycle_determination_timing"] = molgraph.cycle_time
    except AttributeError:
        timings["cycle_determination_timing"] = 0.0

    try:
        timings["unification_timing"] = molgraph.ff_time
    except AttributeError:
        timings["unification_timing"] = 0.0

    try:
        timings["topology_diag_timing"] = molgraph.topology_time
    except AttributeError:
        timings["topology_diag_timing"] = 0.0

    try:
        timings["topography_diag_timing"] = molgraph.topography_time
    except AttributeError:
        timings["topography_diag_timing"] = 0.0

    try:
        timings["chirality_diag_timing"] = molgraph.chirality_time
    except AttributeError:
        timings["chirality_diag_timing"] = 0.0

    data["timings"] = timings

    data["elements"] = {
        int(i + 1): element
        for i, element in enumerate(molgraph.return_node_data(attribute="elements"))
    }

    data["atomic_numbers"] = {
        int(i + 1): int(Z)
        for i, Z in enumerate(molgraph.return_node_data(attribute="atomic_numbers"))
    }

    data["degrees"] = {
        int(i + 1): int(degree) for i, degree in enumerate(molgraph.get_degree())
    }

    try:
        data["priorities"] = {
            int(i + 1): int(round(prio))
            for i, prio in enumerate(molgraph.return_node_data(attribute="priorities"))
        }

    except TypeError:

        data["priorities"] = []

    try:

        data["fragment_priorities"] = {
            int(i + 1): int(prio) for i, prio in enumerate(molgraph.fragment_priorities)
        }

    except AttributeError:

        data["fragment_priorities"] = []

    data["vsepr_classes"] = {
        i + 1: str(vsepr)
        for i, vsepr in enumerate(molgraph.return_node_data(attribute="vsepr_classes"))
    }

    data["single_bonds"] = sorted(
        [
            sorted([int(edge[0] + 1), int(edge[1] + 1)])
            for edge in molgraph.return_edges(attributes=["bo"], values=[1])
        ],
        key=lambda x: (x[0], x[1]),
    )

    data["double_bonds"] = sorted(
        [
            sorted([int(edge[0] + 1), int(edge[1] + 1)])
            for edge in molgraph.return_edges(attributes=["bo"], values=[2])
        ],
        key=lambda x: (x[0], x[1]),
    )

    data["triple_bonds"] = sorted(
        [
            sorted([int(edge[0] + 1), int(edge[1] + 1)])
            for edge in molgraph.return_edges(attributes=["bo"], values=[3])
        ],
        key=lambda x: (x[0], x[1]),
    )

    all_cycle_ids = molgraph.return_edge_data(
        attribute="cycle_id", nodes_visible=True, include_all=False
    )

    edges = molgraph.return_visible_edges(nodes_visible=True)

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

    edges_in_cycles = [
        [
            edge
            for i, edge in enumerate(edges)
            if (all_cycle_ids[i] is not None) and (cycle in all_cycle_ids[i])
        ]
        for cycle in cycles
    ]

    nodes_in_cycle = [
        set([node for edge in edges_in_cycle for node in edge])
        for edges_in_cycle in edges_in_cycles
    ]

    data["cycles"] = sorted(
        [sorted([int(node + 1) for node in cycle]) for cycle in nodes_in_cycle],
        key=lambda x: [x[i] for i in range(len(x))],
    )

    try:
        data["fragment_data"] = {
            int(fragment_id + 1): _debug_fragment(molgraph, fragment_id)
            for fragment_id in molgraph.fragments_data.keys()
        }

    except AttributeError:

        data["fragment_data"] = {"1": "No fragments found."}

    return data


def _debug_fragment(molgraph: MolGraph, fragment_id: int) -> dict:
    """
    Returns the debug data of a fragment.

    Args:
        molgraph (MolGraph): MolGraph object.
        fragment_id (int): fragment id for which the debug data should be returned.

    Returns:
        dict: debug data of the fragment.
    """

    molgraph.set_nodes_visible(include_all=True)

    fragment_ids = molgraph.return_node_data(attribute="fragment_id")

    visible_nodes = molgraph.return_visible_nodes()

    fragment = [node for node in visible_nodes if fragment_ids[node] == fragment_id]

    fragment_data = molgraph.fragments_data[fragment_id]

    adjacent_nodes = molgraph.return_adjacent_nodes(core_nodes=fragment)

    molgraph.set_nodes_visible(visible_nodes=fragment)

    all_angle_constraints = molgraph.return_node_data(attribute="angle_constraints")

    all_edge_dihedral_constraints = molgraph.return_edge_data(
        attribute="dihedral_constraints", nodes_visible=True
    )

    all_node_dihedral_constraints = molgraph.return_node_data(
        attribute="dihedral_constraints"
    )

    visible_edges = molgraph.return_visible_edges(nodes_visible=True)

    molgraph.set_nodes_visible(visible_nodes=fragment + adjacent_nodes)

    visible_nodes = molgraph.return_visible_nodes()

    n_atoms = molgraph.return_n_atoms()

    coulomb_matrix = fragment_data["coulomb_matrix"]

    atomic_numbers = molgraph.return_node_data(attribute="atomic_numbers")

    degrees = molgraph.get_degree()

    fragment_barcode = [
        int(x)
        for x in _calculate_barcode(n_atoms, coulomb_matrix, atomic_numbers, degrees)[0]
    ]

    index_conversion = {node: index for index, node in enumerate(visible_nodes)}

    angles, ideal_angles = _get_angle_constraints(
        fragment, all_angle_constraints, index_conversion
    )

    angle_data = [
        {"nodes": [int(node + 1) for node in angle], "value": value}
        for angle, value in zip(angles, ideal_angles)
    ]

    dihedrals, ideal_dihedrals = _get_dihedral_constraints(
        visible_edges,
        fragment,
        all_edge_dihedral_constraints,
        all_node_dihedral_constraints,
        index_conversion,
    )

    dihedral_data = [
        {"nodes": [int(node + 1) for node in dihedral], "value": value}
        for dihedral, value in zip(dihedrals, ideal_dihedrals)
    ]

    bonds = [
        [index_conversion[bond[0]], index_conversion[bond[1]]]
        for bond in molgraph.return_visible_edges(nodes_visible=True)
    ]

    ideal_bonds = molgraph.return_edge_data(attribute="dij", nodes_visible=True)

    bond_data = [
        {"nodes": [int(node + 1) for node in bond], "value": value}
        for bond, value in zip(bonds, ideal_bonds)
    ]

    fragment_data["atoms_in_fragment"] = [
        int(node + 1) for node in fragment_data["visible_nodes"]
    ]

    fragment_data["fragment_barcode"] = fragment_barcode

    fragment_data["final_geometry"] = [
        [x for x in row] for row in fragment_data["final_geometry"]
    ]

    fragment_data["elements"] = [x for x in fragment_data["elements"]]

    fragment_data["bonds"] = bond_data

    fragment_data["angles"] = angle_data

    fragment_data["dihedrals"] = dihedral_data

    fragment_data.pop("nodes")

    fragment_data.pop("coulomb_matrix")

    fragment_data.pop("visible_nodes")

    return fragment_data


def _print_debug_data(data: dict, debug_path: str) -> None:
    """
    Writes debug.json file containing all information about the molecule.

    Args:
        data (dict): debug data.
        debug_path (str): path to debug file.

    Returns:
        None
    """

    with open(os.path.join(debug_path, "output.json"), "w") as f:

        json.dump(data, f, indent=4)
