import numpy as np
import time
from decimal import Decimal, ROUND_HALF_UP
from molbar.molgraph.fraggraph import FragGraph
from molbar.molgraph.molgraph import MolGraph
from collections import Counter
from importlib.metadata import version


def get_topology_barcode(molgraph):
    molgraph.set_nodes_visible(include_all=True)
    start = time.time()
    molgraph.topology_spectrum, molgraph.topology_orbitals = _get_topology_barcode(
        molgraph
    )
    end = time.time()
    molgraph.topology_time = np.round(end - start, 6)


def get_heavy_atom_topology_barcode(molgraph):
    molgraph.set_nodes_visible(include_all=True)
    molgraph.heavy_atom_topology_spectrum = _get_heavy_atom_topology_barcode(molgraph)


def get_topography_barcode(molgraph):
    molgraph.set_nodes_visible(include_all=True)
    start = time.time()
    molgraph.topography_spectrum, molgraph.topography_orbitals = (
        _get_topography_barcode(molgraph)
    )
    end = time.time()
    molgraph.topography_time = np.round(end - start, 6)


def get_absconf_barcode(molgraph):
    molgraph.set_nodes_visible(include_all=True)
    start = time.time()
    molgraph.absconfiguration_spectrum = _get_absconf_barcode(molgraph)
    end = time.time()
    molgraph.chirality_time = np.round(end - start, 6)


def get_molecular_formula(molgraph):
    elements = molgraph.return_node_data(attribute="elements")
    molgraph.molecular_formula = _get_molecular_formular(elements)


def concat_barcode(molgraph: MolGraph) -> str:
    """
    Returns MolBar based on the provided MolGraph object.
    Args:
        molgraph (MolGraph): A MolGraph object.

    Returns:
        A string containing the barcode of the molecule.
    """
    mb_version = version("molbar")
    molecular_formula = molgraph.molecular_formula
    total_charge = int(molgraph.total_charge)
    topology_spectrum_str = _get_rounded_spectrum_as_str(molgraph.topology_spectrum)
    heavy_atom_topology_spectrum_str = _get_rounded_spectrum_as_str(
        molgraph.heavy_atom_topology_spectrum
    )
    topography_spectrum_str = _get_rounded_spectrum_as_str(molgraph.topography_spectrum)
    absconfiguration_spectrum_str = _get_rounded_spectrum_as_str(
        molgraph.absconfiguration_spectrum
    )
    return f"MolBar | {mb_version} | {molecular_formula} | {total_charge} | {topology_spectrum_str} | {heavy_atom_topology_spectrum_str} | {topography_spectrum_str} | {absconfiguration_spectrum_str} "


def _get_molecular_formular(elements: list) -> str:
    """
    Creates the molecular formula from a list of elements.

    Args:
        elements (list): A list of elements.

    Returns:
        str: molecular formula
    """
    element_counts = Counter(elements).items()
    sorted_elements = _sort_element_counts(element_counts)
    return "".join(
        f"{element}{count}" if count > 1 else element
        for element, count in sorted_elements
    )


def concat_topology_barcodes(molgraph):
    """
    Calculates the topology barcodes only of a molecule.

    Args:
        molgraph (MolGraph): A MolGraph object.

    Returns:
        A list containing the barcode of the molecule.
    """
    mb_version = version("molbar")
    molecular_formula = molgraph.molecular_formula
    topology_spectrum_str = _get_rounded_spectrum_as_str(molgraph.topology_spectrum)
    heavy_atom_topology_spectrum_str = _get_rounded_spectrum_as_str(
        molgraph.heavy_atom_topology_spectrum
    )
    return f"TopoBar | {mb_version} | {molecular_formula} | {topology_spectrum_str} | {heavy_atom_topology_spectrum_str}"


def _get_topology_barcode(molgraph: MolGraph) -> np.ndarray:
    """
    Calculates the topology barcode based on the MolGraph object.

    Args:
        molgraph (MolGraph): A MolGraph object.

    Returns:
        np.ndarray: A list containing the topology spectrum of the molecule.
    """

    molgraph.set_nodes_visible(include_all=True)

    n_atoms = molgraph.return_n_atoms()

    adjacency_matrix = molgraph.return_cn_matrix()

    atomic_numbers = molgraph.return_node_data(attribute="atomic_numbers")

    degrees = molgraph.get_degree()

    eps, C = _calculate_barcode(
        n_atoms, adjacency_matrix, atomic_numbers, degrees, return_eigenvec=True
    )

    return eps, C


def _get_heavy_atom_topology_barcode(molgraph: MolGraph) -> np.ndarray:
    """
    Calculates the heavy atom topology barcode based on the MolGraph object.

    Args:
        molgraph (MolGraph): A MolGraph object.

    Returns:
        np.ndarray: A list containing the heavy atom topology spectrum of the molecule.
    """

    elements = molgraph.return_node_data(attribute="elements")

    heavy_atom_indices = [i for i, element in enumerate(elements) if element != "H"]

    molgraph.set_nodes_visible(visible_nodes=heavy_atom_indices)

    n_atoms = molgraph.return_n_atoms(include_all=False)

    adjacency_matrix = molgraph.return_cn_matrix(include_all=False)

    atomic_numbers = molgraph.return_node_data(attribute="atomic_numbers")

    degrees = molgraph.get_degree(include_all=False)

    eps, _ = _calculate_barcode(n_atoms, adjacency_matrix, atomic_numbers, degrees)

    return eps


def _get_topography_barcode(molgraph: MolGraph) -> np.ndarray:
    """
    Calculates the topography barcode based on the MolGraph object.

    Args:
        molgraph (MolGraph): A MolGraph object.

    Returns:
        np.ndarray: A list containing the topography spectrum of the molecule.
    """

    molgraph.set_nodes_visible(include_all=True)

    n_atoms = molgraph.return_n_atoms()

    cm_matrix = molgraph.unified_coulomb_matrix

    atomic_numbers = molgraph.return_node_data(attribute="atomic_numbers")

    degrees = molgraph.get_degree()

    eps, C = _calculate_barcode(
        n_atoms, cm_matrix, atomic_numbers, degrees, return_eigenvec=True
    )

    return eps, C


def _get_absconf_barcode(molgraph: MolGraph) -> np.ndarray:
    """
    Calculates the absolute configuration barcode based on the MolGraph object.

    Args:

        molgraph (MolGraph): A MolGraph object.

    Returns:

        np.ndarray: A list containing the absolute configuration spectrum of the molecule.
    """
    molgraph.set_nodes_visible(include_all=True)

    fraggraph = FragGraph(molgraph)

    if fraggraph.n_fragments == 0:

        return np.array([])

    fragment_distance_matrix = fraggraph.get_fragment_distance_matrix()

    np.fill_diagonal(fragment_distance_matrix, 1)

    fragment_coulomb_matrix = 1 / fragment_distance_matrix

    n_fragments = fraggraph.return_n_fragments()

    abs_conf_indices = fraggraph.get_abs_conf_indices()

    fragment_priorities = fraggraph.return_fragment_priorities()

    molgraph.fragment_priorities = np.array(fragment_priorities, dtype=int)

    return _calculate_absconf_barcode(
        n_fragments, fragment_coulomb_matrix, abs_conf_indices, fragment_priorities
    )


def _calculate_barcode(
    n_atoms, matrix, nuclear_charges, degrees, return_eigenvec=False
) -> np.ndarray:
    """
    Calculates the adjacency spectrum of a molecule.

    Args:
        cn_matrix (np.ndarray): A matrix representing the coordination numbers of each pair of atoms.
        charges (list): A list of the charges for each atom in the molecule.
        cn (list): A list of the coordination numbers for each atom in the molecule.

    Returns:
        An array containing the adjacency spectrum of the molecule, rounded to one decimal place and
        multiplied by 10 and converted to integers.
    """

    # Calculate the matrix diagonal values
    matrix_diagonal = (degrees + 1) * nuclear_charges

    # Create the diagonal barcode matrix
    diagonal_barcode_matrix = np.diag(matrix_diagonal)

    # Create the off-diagonal barcode matrix
    off_diagonal_barcode_matrix = np.zeros((n_atoms, n_atoms))
    off_diagonal_barcode_matrix[:, :] = matrix_diagonal
    off_diagonal_barcode_matrix = 0.5 * (
        off_diagonal_barcode_matrix + off_diagonal_barcode_matrix.transpose()
    )

    # Create the barcode matrix
    barcode_matrix = off_diagonal_barcode_matrix * matrix + diagonal_barcode_matrix

    barcode_matrix = symmetrize(barcode_matrix)

    # Calculate the eigenvalues of the barcode matrix and round them to one decimal point
    if return_eigenvec:
        try:
            eps, C = np.linalg.eigh(barcode_matrix)
        except np.linalg.LinAlgError:
            eps = np.linalg.eigvalsh(barcode_matrix)
            C = compute_eigenvectors(barcode_matrix, eps)
        return eps * 10, C
    else:
        eps = np.linalg.eigvalsh(barcode_matrix)
        return eps * 10, None


def symmetrize(matrix):
    return (matrix + matrix.T) / 2


def compute_eigenvectors(A, eigenvalues):
    """
    Compute the eigenvectors of a matrix A given its eigenvalues.
    Used in case when the np.linalg.eigh function fails to compute the eigenvectors.
    Args:
        A (_type_): _description_
        eigenvalues (_type_): _description_

    Returns:
        _type_: _description_
    """

    n = A.shape[0]
    eigenvectors = []
    identity = np.eye(n)
    for eigenvalue in eigenvalues:
        # Solve (A - λI)v = 0
        M = A - eigenvalue * identity
        # Find a non-trivial solution to Mv = 0, using np.linalg.svd
        U, S, Vh = np.linalg.svd(M)
        null_space = Vh.T[:, -1]
        eigenvectors.append(null_space)
    return np.array(eigenvectors).T


def _calculate_absconf_barcode(n_fragments, matrix, abs_conf_indices, priorities):

    barcode_matrix = np.zeros((n_fragments, n_fragments))

    matrix_diagonal = abs_conf_indices * priorities

    atomic_matrix = np.zeros((n_fragments, n_fragments))

    atomic_matrix[:, :] = matrix_diagonal

    atomic_matrix = 0.5 * (atomic_matrix + atomic_matrix.transpose())

    np.fill_diagonal(atomic_matrix, matrix_diagonal)

    barcode_matrix = matrix * atomic_matrix

    spectrum = np.linalg.eigvalsh(barcode_matrix) * 10

    return spectrum


def _sort_element_counts(element_counts):
    # Sort the element_counts list by the first element, which is the element symbol
    element_counts_sorted = sorted(element_counts, key=lambda x: x[0])

    # Move the tuple with the "H" element symbol to the end of the list
    for i in range(len(element_counts_sorted)):
        if element_counts_sorted[i][0] == "H":
            element_counts_sorted.append(element_counts_sorted.pop(i))
            break

    return element_counts_sorted


def _get_rounded_spectrum_as_str(spectrum):

    decimal_spectrum = [
        Decimal(eps).quantize(Decimal("1.000"), rounding=ROUND_HALF_UP)
        for eps in spectrum
    ]

    rounded_spectrum = [
        eps.quantize(Decimal("1"), rounding=ROUND_HALF_UP) for eps in decimal_spectrum
    ]

    spectrum_str = " ".join([str(eps) for eps in rounded_spectrum])

    if "-0" in spectrum_str:

        spectrum_str = spectrum_str.replace("-0", "0")

    return spectrum_str
