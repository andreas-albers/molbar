import os
import random
from molbar.barcode import get_molbar_from_file, get_molbar_from_coordinates
from molbar.molgraph.molgraph import MolGraph


def test_c76():
    script_path = os.path.dirname(os.path.abspath(__file__))
    filepath1 = os.path.join(script_path, "../example/c76_m.xyz")
    molbar1 = get_molbar_from_file(filepath1)

    filepath2 = os.path.join(script_path, "../example/c76_p.xyz")
    molbar2 = get_molbar_from_file(filepath2)
    print(molbar1)
    print(molbar2)
    assert molbar1 != molbar2, "C76_m and C76_p should have different molecular barcodes."

def test_ru():
    script_path = os.path.dirname(os.path.abspath(__file__))
    filepath1 = os.path.join(script_path, "../example/ru_001.xyz")
    molbar1 = get_molbar_from_file(filepath1)
    filepath2 = os.path.join(script_path, "../example/ru_003.xyz")
    molbar2 = get_molbar_from_file(filepath2)
    filepath3 = os.path.join(script_path, "../example/ru_006.xyz")
    molbar3 = get_molbar_from_file(filepath3)
    assert molbar1 != molbar2, "Ru_001 and Ru_003 should have different molecular barcodes."
    assert molbar1 == molbar3, "Ru_001 and Ru_006 should have different molecular barcodes."
    assert molbar2 != molbar3, "Ru_003 and Ru_006 should have different molecular barcodes."


def test_permutation():
    script_path = os.path.dirname(os.path.abspath(__file__))
    filepath1 = os.path.join(script_path, "../example/ru_001.xyz")

    molGraph = MolGraph()
    molGraph.from_file(filepath1)
    coordinates = molGraph.return_node_data("coordinates")
    elements = molGraph.return_node_data("elements")
    n_atoms = molGraph.return_n_atoms()
    indices = [i for i in range(n_atoms)]
    molbars = []
    for i in range(10):
        random.shuffle(indices)
        permutated_coordinates = coordinates[indices]
        permutated_elements = elements[indices]
        molbar = get_molbar_from_coordinates(permutated_coordinates, permutated_elements)
        molbars.append(molbar)

    assert len(set(molbars)) == 1, "Permutations should have the same molecular barcode."