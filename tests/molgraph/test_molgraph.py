import pytest
import os
import numpy as np
from molbar.molgraph.molgraph import MolGraph  # Adjust the import based on your file structure

@pytest.fixture
def molgraph():
    return MolGraph()

def test_initialization(molgraph):
    assert molgraph.total_charge == 0
    assert molgraph.constraints == {}
    assert molgraph.unified_distance_matrix.size == 0
    assert molgraph.cycle_nodes == []
    assert molgraph.visible_repulsion_nodes == []
    assert molgraph.fragment_data == {}
    assert molgraph.is_2D is False

def test_from_coordinates(molgraph):
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    elements = np.array(["H", "H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    assert molgraph.return_n_atoms() == 3

def test_from_coordinates_and_atomic_numbers(molgraph):
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    elements = np.array([1, 1, 1])
    molgraph.from_coordinates(coordinates, elements)
    assert molgraph.return_node_data(attribute="elements")[0] == "H"

def test_from_file(molgraph):
    script_path = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_path, "../../example/binol_m.xyz")
    molgraph.from_file(filepath)
    assert molgraph.return_n_atoms() == 36

def test_return_n_atoms(molgraph):
    coordinates = np.array([[0, 0, 0], [1, 0, 0]])
    elements = np.array(["H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    assert molgraph.return_n_atoms() == 2

def test_return_cn_matrix(molgraph):
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    elements = np.array(["H", "H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    molgraph.define_edges()
    cn_matrix = molgraph.return_cn_matrix(include_all=True)
    assert cn_matrix.shape == (3, 3)

def test_get_degree(molgraph):
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    elements = np.array(["H", "H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    molgraph.define_edges()
    degrees = molgraph.get_degree(include_all=True)
    assert len(degrees) == 3

def test_connected_components(molgraph):
    coordinates = np.array([[0, 0, 0], [0.5, 0, 0]])
    elements = np.array(["H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    molgraph.define_edges()
    components = molgraph.connected_components()
    assert len(components) == 1

def test_rigid_fragmentation(molgraph):
    script_path = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_path, "../../example/binol_m.xyz")
    molgraph.from_file(filepath)
    molgraph.define_edges()
    molgraph.detect_cycles()
    molgraph.rigid_fragmentation()
    assert "fragment_id" in molgraph.nodes[0]

def test_set_visible(molgraph):
    coordinates = np.array([[0, 0, 0], [0.5, 0, 0]])
    elements = np.array(["H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    molgraph.set_all_visible()
    visible_nodes = [data["visible"] for node, data in molgraph.nodes(data=True)]
    assert all(visible_nodes)


if __name__ == "__main__":
    
    molgraph = MolGraph()
    test_from_coordinates_and_atomic_numbers(molgraph)