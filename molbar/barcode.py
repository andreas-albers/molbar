import sys
import copy
import numbers
from tqdm import tqdm
from typing import Union
from joblib import Parallel, delayed
from molbar.runcalc.runcalc import _run_molbar_calculation
import warnings

warnings.filterwarnings(
    "ignore", category=UserWarning, module="molbar.indices.absolute_conf"
)
warnings.filterwarnings(
    "ignore", category=UserWarning, module="molbar.molgraph.nodes.vsepr"
)


def get_molbars_from_coordinates(
    list_of_coordinates: list,
    list_of_elements: list,
    total_charges=0,
    list_of_edges=None,
    list_of_bond_orders=None,
    list_of_formal_charges = None,
    list_of_stereo_settings = None,
    return_data=False,
    threads=1,
    timing=False,
    input_constraints=None,
    progress=False,
    mode="mb",
) -> Union[list, Union[str, dict]]:
    """
    Sets up the MolGraph object for the calculation of the molecular barcode for a list of coordinates with a list of constraints.
    And runs the calculation in parallel if more than one thread is used.

    Args:

          list_of_coordinates (list): A list of molecular geometries provided by atomic Cartesian coordinates with shape (n_molecules, n_atoms, 3).
          list_of_elements (list): A list of element lists for each molecule in the list_of_coordinates with shape (n_molecules, n_atoms). Either the element symbols or atomic numbers can be used.
          return_data (bool): Whether to return MolBar data.
          threads (int): Number of threads to use for the calculation. If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.
          timing (bool):  Whether to print the duration of this calculation.
          input_constraints (list, optional): A list of constraints for the calculation. Each constraint in that list is a Python dict as shown in the documentation.
          progress (bool): Whether to show a progress bar.
          mode (str): Whether to calculate the molecular barcode ("mb") or the topology part of the molecular barcode ("topo").

      Returns:

          Union[list, Union[str, dict]]: Either MolBar or the MolBar and MolBar data.

    """

    if input_constraints is None:
        input_constraints = [None] * len(list_of_elements)
    if list_of_edges is None:
        list_of_edges = [None] * len(list_of_elements)
    if list_of_bond_orders is None:
        list_of_bond_orders = [None] * len(list_of_elements)
    if list_of_formal_charges is None:
        list_of_formal_charges = [None] * len(list_of_elements)
    if list_of_stereo_settings is None:
        list_of_stereo_settings = [None] * len(list_of_stereo_settings)
    if isinstance(total_charges, numbers.Integral):
        total_charges = [total_charges] * len(list_of_elements)
    assert len(list_of_elements) == len(
        list_of_coordinates
    ), "The number of elements and coordinates must be the same."
    assert len(input_constraints) == len(
        list_of_elements
    ), "The number of constraints must be the same as the number of molecules."
    assert len(total_charges) == len(
        list_of_elements
    ), "The number of total_charges must be the same as the number of molecules"

    if progress:

        results = Parallel(n_jobs=threads)(
            delayed(_run_molbar_calculation)(
                coordinates=coord,
                elements=elem,
                total_charge=total_charge,
                edges=edges,
                bond_orders=bond_orders,
                formal_charges=formal_charges,
                stereo_settings=stereo_settings,
                single_constraint=copy.deepcopy(single_constraint),
                debug=return_data,
                timing=timing,
                mode=mode,
            )
            for single_constraint, coord, elem, total_charge, edges, bond_orders, formal_charges, stereo_settings in tqdm(
                zip(
                    input_constraints,
                    list_of_coordinates,
                    list_of_elements,
                    total_charges,
                    list_of_edges,
                    list_of_bond_orders,
                    list_of_formal_charges,
                    list_of_stereo_settings,
                ),
                total=len(list_of_elements),
                file=sys.stdout,
            )
        )

    else:

        results = Parallel(n_jobs=threads)(
            delayed(_run_molbar_calculation)(
                coordinates=coord,
                elements=elem,
                total_charge=total_charge,
                edges=edges,
                bond_orders=bond_orders,
                formal_charges=formal_charges,
                single_constraint=copy.deepcopy(single_constraint),
                debug=return_data,
                timing=timing,
                mode=mode,
            )
            for single_constraint, coord, elem, total_charge, edges, bond_orders, formal_charges in zip(
                input_constraints,
                list_of_coordinates,
                list_of_elements,
                total_charges,
                list_of_edges,
                list_of_bond_orders,
                list_of_formal_charges,
            )
        )

    return results


def get_molbar_from_coordinates(
    coordinates: list,
    elements: list,
    total_charge=0,
    edges=None,
    bond_orders=None,
    formal_charges = None,
    stereo_settings = None,
    return_data=False,
    timing=False,
    input_constraint=None,
    mode="mb",
) -> Union[str, dict]:
    """
    Sets up the MolGraph object for the calculation of the molecular barcode for a single set of coordinates with a single constraint.

    Args:

          coordinates (list): Molecular geometry provided by atomic Cartesian coordinates with shape (n_atoms, 3).
          elements (list): A list of elements in that molecule. Either the element symbols or atomic numbers can be used.
          return_data (bool): Whether to return MolBar data.
          timing (bool): Whether to print the duration of this calculation.
          input_constraint (dict, optional): A dict of extra constraints for the calculation. See documentation for more information. USED ONLY IN EXCEPTIONAL CASES.
          mode (str): Whether to calculate the molecular barcode ("mb") or only the topology part of the molecular barcode ("topo").

      Returns:

          Union[str, dict]: Either MolBar or the MolBar and MolBar data.

    """

    return _run_molbar_calculation(
        coordinates=coordinates,
        elements=elements,
        total_charge=total_charge,
        edges=edges,
        bond_orders=bond_orders,
        formal_charges=formal_charges,
        stereo_settings=stereo_settings,
        single_constraint=input_constraint,
        debug=return_data,
        timing=timing,
        mode=mode,
    )


def get_molbars_from_files(
    files: list,
    total_charges=0,
    return_data=False,
    threads=1,
    timing=False,
    input_constraints=None,
    progress=False,
    mode="mb",
    write_trj=False,
) -> Union[list, Union[str, dict]]:
    """
    Sets up the MolGraph object for the calculation of the molecular barcode for a list of files with a list of constraints.
    And runs the calculation in parallel if more than one thread is used.

    Args:

          files (list): The list of paths to the files containing the molecule information (either .xyz/.sdf/.mol format).
          return_data (bool): Whether to return MolBar data.
          threads (int): Number of threads to use for the calculation. If you need to process multiple molecules at once, it is recommended to use this function and specify the number of threads that can be used to process multiple molecules simultaneously.
          timing (bool):  Whether to print the duration of this calculation.
          input_constraints (list, optional): A list of file paths to the input files for the calculation. Each constrained is specified by a file path to a .yml file, as shown in the documentation.
          progress (bool): Whether to show a progress bar.
          mode (str): Whether to calculate the molecular barcode ("mb") or the topology part of the molecular barcode ("topo").
          write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.

      Returns:

          Union[list, Union[str, dict]]: Either MolBar or the MolBar and MolBar data.

    """

    if not input_constraints:

        input_constraints = [None] * len(files)

    if isinstance(total_charges, numbers.Integral):
        total_charges = [total_charges] * len(files)
    assert len(input_constraints) == len(
        files
    ), "The number of constraints must be the same as the number of molecules."
    assert len(total_charges) == len(
        files
    ), "The number of total_charges must be the same as the number of molecules"

    if progress:
        results = Parallel(n_jobs=threads)(
            delayed(_run_molbar_calculation)(
                file=file,
                total_charge=charge,
                single_constraint=copy.deepcopy(single_constraint),
                debug=return_data,
                timing=timing,
                mode=mode,
                write_trj=write_trj,
            )
            for single_constraint, charge, file in tqdm(
                zip(input_constraints, total_charges, files),
                total=len(files),
                file=sys.stdout,
            )
        )
    else:
        results = Parallel(n_jobs=threads)(
            delayed(_run_molbar_calculation)(
                file=file,
                total_charge=charge,
                single_constraint=copy.deepcopy(single_constraint),
                debug=return_data,
                timing=timing,
                mode=mode,
                write_trj=write_trj,
            )
            for single_constraint, charge, file in zip(
                input_constraints, total_charges, files
            )
        )
    return results


def get_molbar_from_file(
    file: str,
    total_charge=0,
    return_data=False,
    timing=False,
    input_constraint=None,
    mode="mb",
    write_trj=False,
) -> Union[str, dict]:
    """
    Sets up the MolGraph object for the calculation of the molecular barcode for a single file with a single constraint.

    Args:
          file (str): The path to the file containing the molecule information (either .xyz/.sdf/.mol format).
          return_data (bool): Whether to return MolBar data.
          timing (bool): Whether to print the duration of this calculation.
          input_constraint (dict, optional): A dict of extra constraints for the calculation. See documentation for more information. USED ONLY IN EXCEPTIONAL CASES.
          mode (str): Whether to calculate the molecular barcode ("mb") or only the topology part of the molecular barcode ("topo").
          write_trj (bool, optional): Whether to write a trajectory of the unification process. Defaults to False.

      Returns:

          Union[str, dict]: Either MolBar or the MolBar and MolBar data.

    """

    return _run_molbar_calculation(
        file=file,
        total_charge=total_charge,
        single_constraint=input_constraint,
        debug=return_data,
        timing=timing,
        mode=mode,
        write_trj=write_trj,
    )
