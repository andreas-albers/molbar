import sys
from typing import Union
import time
import numpy as np
from molbar.helper.input import _get_constraints, _transform_constraints
from molbar.helper.debug import _get_debug_path, _prepare_debug_data
from molbar.molgraph.molgraph import MolGraph
from molbar.barcodes.barcodes import get_topology_barcode
from molbar.topography.unification import unify_fragments


def idealize_structure_from_file(
    file: str, return_data=False, timing=False, input_constraint=None, write_trj=False
) -> Union[list, str]:
    """
    Idealizes a whole molecule with the MolBar force field and returns the final energy, coordinates and elements.

    Args:

          file (str): The path to the input file to be processed.
          return_data (bool): Whether to print MolBar data.
          timing (bool): Whether to print the duration of this calculation.
          input_constraint (str): The path to the input file containing the constraint for the calculation. See documentation for more information.
          write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.
      Returns:
          n_atoms (int): Number of atoms in the molecule.
          energy (float): Final energy of the molecule after idealization.
          coordinates (list): Final coordinates of the molecule after idealization.
          elements (list): Elements of the molecule.
          data (dict): Molbar data.
    """

    if input_constraint:

        input_constraint = _get_constraints(input_constraint)

        input_constraint = _transform_constraints(file, input_constraint)

    else:

        input_constraint = {}

    input_constraint["unique_repulsion"] = True

    start = time.time()
    result = _idealize_structure(
        file=file,
        single_constraint=input_constraint,
        debug=return_data,
        from_file=True,
        write_trj=write_trj,
    )
    end = time.time()

    if timing:

        print("Duration [s]: " + str(np.round(end - start, 3)), file=sys.stderr)

    return result


def idealize_structure_from_coordinates(
    coordinates: list,
    elements: list,
    edges=None,
    bond_orders=None,
    formal_charges = None,
    return_data=False,
    timing=False,
    input_constraint=None,
) -> Union[list, str]:
    """
    Idealizes a whole molecule with the MolBar force field and returns the final energy, coordinates and elements.

    Args:
        coordinates (list): Cartesian coordinates of the molecule.
        elements (list): Elements of the molecule.
        return_data (bool, optional): Whether to return more information. Defaults to False.
        timing (bool, optional): Whether to print timing information. Defaults to False.
        input_constraint (dict, optional): The constraint for the calculation. See documentation for more information. Defaults to None.
    Returns:
        n_atoms (int): Number of atoms in the molecule.
        energy (float): Final energy of the molecule after idealization.
        coordinates (list): Final coordinates of the molecule after idealization.
        elements (list): Elements of the molecule.
        data (dict):  MolBar data.
    """

    if not input_constraint:

        input_constraint = {}

    start = time.time()
    result = _idealize_structure(
        coordinates=coordinates,
        elements=elements,
        edges=edges,
        bond_orders=bond_orders,
        formal_charges = formal_charges,
        single_constraint=input_constraint,
        debug=return_data,
        write_trj=False,
    )
    end = time.time()

    if timing:

        print("Duration [s]: " + str(np.round(end - start, 3)), file=sys.stderr)

    return result


def _idealize_structure(
    file=None,
    coordinates=None,
    elements=None,
    edges=None,
    bond_orders=None,
    formal_charges = None,
    single_constraint=None,
    debug=False,
    from_file=False,
    write_trj=False,
) -> Union[list, str]:
    """
    Sets up the MolGraph object to idealize a whole molecule with the force field.
    Args:
        file (str): The path to the input file to be processed.
        single_constraint (dict): A dictionary containing the constraint for the calculation. See documentation for more information.
        debug (bool): Whether to print debugging information.
        multiple (bool, optional): If filename should be returned. Needed for several files processed at once. Defaults to False.

    Returns:
        Union[list, str]: Either list with filename and molecular barcode or only molecular barcode.
    """

    # Initialize MolGraph object.
    molgraph = MolGraph()

    molgraph._commandline = write_trj

    if single_constraint == False:

        single_constraint = {}

    # Set the constraint for the calculation.
    molgraph.constraints = single_constraint

    if from_file:
        # Read the input file.
        molgraph.from_file(filepath=file)

    else:

        molgraph.from_coordinates(coordinates=coordinates, elements=elements, edges=edges, bond_orders=bond_orders, formal_charges=formal_charges)

    molgraph.set_nodes_visible(include_all=True)

    if molgraph.constraints.get("set_edges") != False:
        molgraph.define_edges()
        get_topology_barcode(molgraph)
        start = time.time()
        molgraph.calculate_priorities(type="topology")
        end = time.time()
        molgraph.prio_time = np.round(end - start, 6)

    molgraph._add_bond_constraints()

    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)

    # Calculate bond orders if allowed by constraints.
    if molgraph.constraints.get("bond_order_assignment") != False:

        start = time.time()

        molgraph.assign_bond_orders()

        end = time.time()

        molgraph.bo_time = np.round(end - start, 6)

    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)

    # Detects cycles if allowed by constraints.
    if molgraph.constraints.get("cycle_detection") != False:

        start = time.time()

        molgraph.detect_cycles()

        end = time.time()

        molgraph.cycle_time = np.round(end - start, 6)

    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)

    if molgraph.constraints.get("set_angles") != False:

        # Calculate the VSEPR geometry of each nodes.
        molgraph.classify_nodes_geometry(is_adjacent_visible=True)

        # Define angle constraints and internal dihedral constraints based on VSEPR geometry.
        molgraph.constrain_nodes_angles(is_adjacent_visible=True)

    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)

    if molgraph.constraints.get("set_dihedrals") != False:

        # Define dihedral constraints based on double bonds.
        molgraph._add_db_dihedral_constraints()

        # Define dihedral constraints based on cycles.
        molgraph._add_cycle_dihedral_constraints()

    molgraph.set_nodes_visible(include_all=True)

    molgraph.set_edges_visible(include_all=True)

    if molgraph.constraints.get("set_repulsion") != False:

        molgraph._get_repulsion_pairs()

    # Set all nodes to visible.
    molgraph.set_nodes_visible(include_all=True)

    # All atoms are one fragment
    molgraph.add_node_data(
        attribute="fragment_id", new_data=[0] * molgraph.return_n_atoms()
    )

    molgraph.set_nodes_visible(include_all=True)

    molgraph.set_edges_visible(include_all=True)

    # Get debug path based on file path.
    debug_path = _get_debug_path(file)

    # Unification of fragemnts based on force field.
    start = time.time()
    molgraph = unify_fragments(molgraph, write_trj=write_trj, debug_path=debug_path)

    end = time.time()

    molgraph.ff_time = np.round(end - start, 6)

    if debug:

        debug_data = _prepare_debug_data("", molgraph)

    else:

        debug_data = {}

    n_atoms = molgraph.return_n_atoms()

    final_energy = molgraph.fragments_data[0]["final_energy"]

    final_geometry = molgraph.fragments_data[0]["final_geometry"]

    elements = molgraph.fragments_data[0]["elements"]

    return final_energy, final_geometry, elements, n_atoms, debug_data
