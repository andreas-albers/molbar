import pytest
import numpy as np
from scipy import special
from molbar.molgraph.edges.coordination import Coordination

class MockCoordination(Coordination):
    def return_node_data(self, attribute, include_all=False):
        if attribute == "rcov":
            return np.array([0.76, 0.66, 0.67])
        elif attribute == "atomic_numbers":
            return np.array([1, 6, 7])
        elif attribute == "elements":
            return np.array(["H", "C", "N"])
        elif attribute == "rnorm":
            return np.array([0.32, 0.77, 0.70])
        elif attribute == "cn_fak":
            return np.array([0.5, 1.5, 1.0])
        elif attribute == "per_fak":
            return np.array([[1.0, 2.0], [1.0, 2.0], [1.0, 2.0]])
        elif attribute == "en_pauling":
            return np.array([2.1, 2.5, 3.0])
        elif attribute == "normcn":
            return np.array([1.0, 1.0, 1.0])
        else:
            raise ValueError("Unknown attribute")

    def return_n_atoms(self, include_all=False):
        return 3

    def return_distance_matrix(self, include_all=False):
        return np.array([
            [0.0, 1.1, 2.2],
            [1.1, 0.0, 1.5],
            [2.2, 1.5, 0.0]
        ])

@pytest.fixture
def mock_coordination():
    return MockCoordination()

def test_define_rij0(mock_coordination):
    rij0 = mock_coordination._define_rij0()
    assert rij0.shape == (3, 3)
    assert rij0[0, 1] == pytest.approx(1.42)
    assert rij0[1, 2] == pytest.approx(1.33)

def test_return_d4(mock_coordination):
    rij = np.array([
        [0.0, 1.1, 2.2],
        [1.1, 0.0, 1.5],
        [2.2, 1.5, 0.0]
    ])
    rij0 = np.array([1.0, 1.0, 1.0
    ])
    en = np.array([2.1, 2.5, 3.0])
    covCN = mock_coordination._return_d4(rij, rij0, en)
    assert covCN.shape == (3,)
    assert covCN[0] > 0
    assert covCN[1] > 0

def test_get_edges_d3(mock_coordination):
    bonds, rij0 = mock_coordination._get_edges_d3()
    assert len(bonds) > 0
    assert isinstance(rij0, np.ndarray)

def test_get_edges_gfnff(mock_coordination):
    bonds, rij0 = mock_coordination._get_edges_gfnff()
    assert len(bonds) > 0
    assert isinstance(rij0, np.ndarray)

if __name__ == "__main__":
    test_return_d4(MockCoordination())