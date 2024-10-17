import pytest
import json
import networkx as nx
import numpy as np
from molbar.molgraph.nodes.nodes import Nodes

class MockGraph(nx.Graph, Nodes):
    def __init__(self, **attr):
        super().__init__(**attr)
        self._nodes = {}
        self.element_data = {
            "atomic_numbers": {"H": 1, "C": 6, "N": 7},
            "rcov": [0.31, 0.76, 0.71],
            "en_pauling2": [2.20, 2.55, 3.04],
            "cn_fak": [0.5, 1.0, 1.5],
            "per_fak": [0.4, 0.8, 1.2],
            "masses": [1.008, 12.01, 14.01],
            "is_metal": [False, False, False],
            "atomic_valences": [1, 4, 5],
            "normcn": [0.5, 1.0, 1.5],
            "rnorm": [0.3, 0.6, 0.9]
        }

@pytest.fixture
def mock_graph():
    g = MockGraph()
    g.add_nodes_from([
        (0, {"visible": True, "rcov": 0.31, "coordinates": [0, 0, 0]}),
        (1, {"visible": True, "rcov": 0.76, "coordinates": [1, 0, 0]}),
        (2, {"visible": False, "rcov": 0.71, "coordinates": [0, 1, 0]})
    ])
    return g
def test_return_nodes(mock_graph):
    nodes = mock_graph.return_nodes(attributes="rcov", values=0.31, include_all=True)
    assert nodes == [0]

def test_return_node_data(mock_graph):
    rcov = mock_graph.return_node_data(attribute="rcov")
    assert np.array_equal(rcov, [0.31, 0.76])

def test_add_node_data(mock_graph):
    mock_graph.add_node_data("test_attr", [1, 2])
    assert mock_graph.nodes(data=True)[0]["test_attr"] == 1
    assert mock_graph.nodes(data=True)[1]["test_attr"] == 2

def test_get_data_for_element(mock_graph):
    data = mock_graph._get_data_for_element("H")
    assert data["atomic_numbers"] == 1
    assert data["rcov"] == 0.31

def test_return_distance_matrix(mock_graph):
    dist_matrix = mock_graph.return_distance_matrix()
    assert dist_matrix.shape == (2, 2)

def test_set_nodes_visible(mock_graph):
    mock_graph.set_nodes_visible(attributes="rcov", values=0.31)
    assert mock_graph.nodes(data=True)[0]["visible"]

def test_set_nodes_unvisible(mock_graph):
    mock_graph.set_nodes_unvisible(attributes="rcov", values=0.31)
    assert not mock_graph.nodes(data=True)[0]["visible"]

def test_return_adjacent_nodes(mock_graph):
    mock_graph.add_edge(0, 1)
    adjacent_nodes = mock_graph.return_adjacent_nodes(core_nodes=0)
    assert 1 in adjacent_nodes

def test_return_visible_nodes(mock_graph):
    visible_nodes = mock_graph.return_visible_nodes()
    assert visible_nodes == [0, 1]

if __name__ == "__main__":

    g = MockGraph()
    g.add_nodes_from([
        (0, {"visible": True, "rcov": 0.31, "coordinates": [0, 0, 0]}),
        (1, {"visible": True, "rcov": 0.76, "coordinates": [1, 0, 0]}),
        (2, {"visible": False, "rcov": 0.71, "coordinates": [0, 1, 0]})
    ])

    nx.set_node_attributes(g, {0: {"visible": True, "rcov": 0.31, "coordinates": [0, 0, 0]}})
    
    test_set_nodes_visible(g)