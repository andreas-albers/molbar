import os
import json
import shutil
from typing import Union
import numpy as np
import time
from molbar.molgraph.molgraph import MolGraph
from molbar.helper.debug import _get_debug_path, _print_debug_data, _mkdir_debug_path
from molbar.topography.unification import unify_fragments
from molbar.barcodes.barcodes import (
    concat_barcode,
    concat_topology_barcodes,
    get_molecular_formula,
    get_topology_barcode,
    get_heavy_atom_topology_barcode,
    get_topography_barcode,
    get_absconf_barcode,
)
from molbar.indices.determine_chirality import determine_absolute_configuration
from molbar.helper.debug import _prepare_debug_data


def determine_graph(molgraph: MolGraph, mode="mb"):
    """
    Determines the graph of the molecule.
    That includes the following steps:
    - Define edges
    - Calculate priorities
    - Detect cycles
    - Assign bond orders

    Args:
        molgraph (MolGraph): The MolGraph object.
        mode (str, optional): The mode of the calculation. Defaults to "mb".
        If "topo" is chosen, only the topology without cycles and bond orders of the molecular barcode is calculated.
    """
    _execute_method_with_timings(molgraph.define_edges, include_all=True)
    get_topology_barcode(molgraph)
    _, molgraph.prio_time = _execute_method_with_timings(
        molgraph.calculate_priorities, type="topology", include_all=True
    )
    if mode == "mb":
        molgraph.set_nodes_visible(include_all=True)
        molgraph.set_edges_visible(include_all=True)
        _, molgraph.cycle_time = _execute_method_with_timings(molgraph.detect_cycles)
        molgraph.set_nodes_visible(include_all=True)
        molgraph.set_edges_visible(include_all=True)
        _, molgraph.bo_time = _execute_method_with_timings(molgraph.assign_bond_orders)


def determine_3D_properties(molgraph: MolGraph):
    """
    Determines the 3D properties of the molecule.
    That includes the following steps:
    - Calculate repulsion pairs
    - Classify nodes geometry

    Args:
        molgraph (MolGraph): The MolGraph object.
    """
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph._get_repulsion_pairs()
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph.classify_nodes_geometry(is_adjacent_visible=True)


def determine_constraints(molgraph: MolGraph):
    """ "
    Determines the constraints of the molecule.
    That includes the following steps:
    - Constrain nodes angles
    - Define dihedral constraints based on double bonds
    - Define dihedral constraints based on allene bonds
    - Define dihedral constraints based on cycles

    Args:
        molgraph (MolGraph): The MolGraph object.

    """
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph.constrain_nodes_angles(is_adjacent_visible=True)
    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    # Define dihedral constraints based on double bonds.
    molgraph._add_db_dihedral_constraints()
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph._add_allene_dihedral_constraints()
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    # Define dihedral constraints based on cycles.
    molgraph._add_cycle_dihedral_constraints()


def get_unified_fragments(molgraph: MolGraph, write_trj=False, debug_path=None):
    """
    Unifies the fragments of the molecule.
    That includes the following steps:
    - Fragmentate the molecule
    - Unify the fragments

    Args:
        molgraph (MolGraph): The MolGraph object.
        write_trj (_type_): Whether to write a trajectory of the unification process.
    """
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph.rigid_fragmentation(include_all=True)
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    molgraph, molgraph.ff_time = _execute_method_with_timings(
        unify_fragments, molgraph, write_trj=write_trj, debug_path=debug_path
    )


def get_3D_barcode(molgraph: MolGraph):
    """
    Calculates the barcode that describes the 3D shape of the molecule.
    That includes the following steps:
    - Get the topography barcode
    - Determine the absolute configuration
    - Get the absolute configuration barcode

    Args:
        molgraph (MolGraph): The MolGraph object.

    Returns:
        str: The 3D barcode of the molecule.
    """
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    get_topography_barcode(molgraph)
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    determine_absolute_configuration(molgraph)
    molgraph.set_nodes_visible(include_all=True)
    molgraph.set_edges_visible(include_all=True)
    get_absconf_barcode(molgraph)


def _execute_method_with_timings(func, *args, **kwargs):
    start = time.time()
    result = func(*args, **kwargs)
    end = time.time()
    return result, np.round(end - start, 6)


def _get_molbar_from_molGraph(
    file=None,
    coordinates=None,
    elements=None,
    total_charge=0,
    edges=None,
    bond_orders=None,
    formal_charges=None,
    mode="mb",
    single_constraint=None,
    debug=False,
    from_file=False,
    write_trj=False,
) -> Union[list, str]:
    """
    Sets up the MolGraph object for the calculation of the molecular barcode for a single file with a single constraint.
    Args:
        file (str): The path to the input file to be processed.
        single_constraint (dict): A dictionary containing the constraint for the calculation. See documentation for more information.
        debug (bool): Whether to print debugging information.
        multiple (bool, optional): If filename should be returned. Needed for several files processed at once. Defaults to False.
        write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.

    Returns:
        Union[list, str]: Either list with filename and molecular barcode or only molecular barcode.
    """
    # Initialize MolGraph object.
    molgraph = MolGraph(total_charge=total_charge)
    molgraph.file = file
    molgraph._commandline = write_trj
    if not single_constraint:
        single_constraint = {}
    # Set the constraint for the calculation.
    molgraph.constraints = single_constraint
    if from_file:
        # Read the input file.
        molgraph.from_file(filepath=file)
    else:
        molgraph.from_coordinates(coordinates=coordinates, elements=elements, edges=edges, bond_orders=bond_orders, formal_charges=formal_charges)
    if write_trj and file:
        debug_path = _get_debug_path(file)
        if os.path.isdir(debug_path):
            shutil.rmtree(debug_path)
        _mkdir_debug_path(debug_path)
    else:
        debug_path = None
    start_total_time = time.time()
    determine_graph(molgraph, mode=mode)
    get_molecular_formula(molgraph)
    get_heavy_atom_topology_barcode(molgraph)
    if mode == "topo":
        molbar_str = concat_topology_barcodes(molgraph)
        if debug:
            debug_data = _prepare_debug_data(molbar_str, molgraph)
        else:
            debug_data = {}
        if write_trj:
            _print_debug_data(debug_data, debug_path)
        return molbar_str, debug_data
    determine_3D_properties(molgraph)
    determine_constraints(molgraph)
    get_unified_fragments(molgraph, write_trj=write_trj, debug_path=debug_path)
    get_3D_barcode(molgraph)
    molbar_str = concat_barcode(molgraph)
    end_total_time = time.time()
    molgraph.total_time = np.round(end_total_time - start_total_time, 6)
    debug_data = _prepare_debug_data(molbar_str, molgraph) if debug else {}
    if write_trj:
        _print_debug_data(debug_data, debug_path)
    return molbar_str, debug_data
