import time
import sys
import numpy as np
from typing import Union
import numbers
from molbar.helper.input import _get_constraints, _transform_constraints
from molbar.molgraph.modify import _get_molbar_from_molGraph


element_list = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]


def _process_single_constraint(single_constraint=None, file=None):
    """
    Process a single constraint for the calculation.
    Can be a string, a dictionary or None.

    Args:
        single_constraint (Union[str, dict], optional): A single constraint for the calculation. Defaults to None.
        file (str, optional): Filepath to the input file. Defaults to None.

    Raises:
        ValueError: Invalid input constraint.

    Returns:
        Union[dict, bool]: The input constraint.
    """
    if not single_constraint:
        return {}
    elif isinstance(single_constraint, str):
        constraints = _get_constraints(single_constraint)
        return _transform_constraints(file, constraints)
    elif isinstance(single_constraint, dict):
        return single_constraint
    else:
        raise ValueError("Invalid input constraint.")


def _validate_coordinates(coordinates: Union[np.ndarray, list]):
    """
    Validate the input coordinates for the calculation.

    Args:
        coordinates (Union[np.ndarray, list]): The input coordinates for the calculation.

    Raises:
        ValueError: Invalid input for coordinates.
    """

    if isinstance(coordinates, np.ndarray):
        assert np.shape(coordinates)[1] == 3, "Coordinates must be a 3D array."
    elif isinstance(coordinates, list):
        assert all(
            len(coord) == 3 for coord in coordinates
        ), "Coordinates must be a 3D array."
    else:
        raise ValueError("Invalid input for coordinates.")


def _validate_elements(
    elements: Union[list, np.ndarray],
    coordinates: Union[list, np.ndarray],
    allowed_elements: list,
):
    """
    Validate the input elements for the calculation.

    Args:
        elements (Union[list, np.ndarray]): The input elements for the calculation.
        coordinates (Union[list, np.ndarray]): The input coordinates for the calculation.
        allowed_elements (list): A list of allowed elements.

    Raises:
        AssertionError: Number of elements must match number of coordinates.
        AssertionError: Invalid element in elements. Must be something like 'C', 'H', 'O', etc. or an atomic number.
    """
    assert len(elements) == len(
        coordinates
    ), "Number of elements must match number of coordinates."
    assert all(
        (element in allowed_elements) or isinstance(element, numbers.Integral)
        for element in elements
    ), "Invalid element in elements. Must be something like 'C', 'H', 'O', etc. or an atomic number."


def _process_input(file=None, coordinates=None, elements=None, single_constraint=None):
    """
    Process the input for the calculation.

    Args:
        file (str, optional): Filepath to the input file. Defaults to None.
        coordinates (Union[np.ndarray, list], optional): The input coordinates for the calculation. Defaults to None.
        elements (Union[list, np.ndarray], optional): The input elements for the calculation. Defaults to None.
        single_constraint (Union[str, dict], optional): A single constraint for the calculation. Defaults to None.

    Raises:
        ValueError: No input provided, either file or coordinates must be provided.
        ValueError: Both file and coordinates provided, only one input allowed.
        ValueError: File must be a filepath.
        ValueError: Elements must be provided with coordinates.

    Returns:
        Union[dict, bool]: The input constraint and whether the input is from a file.
    """
    input_constraint = _process_single_constraint(single_constraint, file)

    if file is not None:
        if coordinates is not None:
            raise ValueError(
                "Both file and coordinates provided, only one input allowed."
            )
        assert isinstance(file, str), "File must be a filepath."
        from_file = True
    elif coordinates is not None:
        if file is not None:
            raise ValueError(
                "Both file and coordinates provided, only one input allowed."
            )
        from_file = False
        allowed_elements = element_list
        _validate_coordinates(coordinates)
        if elements is not None:
            _validate_elements(elements, coordinates, allowed_elements)
        else:
            raise ValueError("Elements must be provided with coordinates.")
    else:
        raise ValueError(
            "No input provided, either file or coordinates must be provided."
        )

    return input_constraint, from_file


def _run_molbar_calculation(
    file=None,
    coordinates=None,
    elements=None,
    total_charge=0,
    edges=None,
    bond_orders=None,
    formal_charges=None,
    single_constraint=None,
    debug=False,
    mode="mb",
    timing=False,
    write_trj=False,
) -> Union[list, str]:
    """
    Calls the MolBar calculation for a single molecule with a single constraint and collects the results.
    Args:
        file (str, optional): path to coordinate file. Defaults to None.
        coordinates (_type_, optional): _description_. Defaults to None.
        elements (_type_, optional): _description_. Defaults to None.
        total_charge (int, optional): _description_. Defaults to 0.
        single_constraint (_type_, optional): _description_. Defaults to None.
        debug (bool, optional): _description_. Defaults to False.
        mode (str, optional): _description_. Defaults to "mb".
        timing (bool, optional): _description_. Defaults to False.
        write_trj (bool, optional): _description_. Defaults to False.
    Raises:
        ValueError: _description_
    Returns:
        Union[list, str]: _description_
    """
    input_constraint, from_file = _process_input(
        file, coordinates, elements, single_constraint
    )

    start = time.time()
    # try:
    barcode, debug_data = _get_molbar_from_molGraph(
        file,
        coordinates,
        elements,
        total_charge,
        edges=edges,
        bond_orders=bond_orders,
        formal_charges=formal_charges,
        single_constraint=input_constraint,
        debug=debug,
        mode=mode,
        from_file=from_file,
        write_trj=write_trj,
    )
    # except Exception as e:
    #     print("Error: " + str(e), file=sys.stderr)
    #     barcode = "Error in MolBar generation. Check STDERR for more information."
    #     debug_data = None
    end = time.time()
    if timing:
        print("Duration [s]: " + str(np.round(end - start, 3)), file=sys.stderr)
    if debug:
        return barcode, debug_data
    else:
        return barcode
