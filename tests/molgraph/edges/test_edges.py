import pytest
import numpy as np
import networkx as nx
from molbar.molgraph.molgraph import MolGraph  # Adjust the import based on your file structure

@pytest.fixture
def molgraph():
    return MolGraph()

@pytest.fixture
def mock_graph():
    molgraph = MolGraph()
    coordinates = np.array([[0, 0, 0], [0.5, 0, 0], [1., 0, 0]])
    elements = np.array(["H", "H", "H"])
    molgraph.from_coordinates(coordinates, elements)
    return molgraph

def test_define_edges(mock_graph):
    mock_graph.define_edges()
    assert len(mock_graph.edges) == 2
    for edge in mock_graph.edges(data=True):
        assert edge[2]['bo'] == 1
        assert edge[2]['visible']
        assert 'dij' in edge[2]
        assert 'rigid' in edge[2]

def test_return_edges(mock_graph):
    mock_graph.define_edges()
    edges = mock_graph.return_edges(attributes='bo', values=1)
    assert len(edges) == 2

def test_return_edges_with_nodes_visible(mock_graph):
    mock_graph.define_edges()
    edges = mock_graph.return_edges(attributes='bo', values=1, nodes_visible=True)
    assert len(edges) == 2

def test_return_edges_with_include_all(mock_graph):
    mock_graph.define_edges()
    edges = mock_graph.return_edges(attributes='bo', values=1, include_all=True)
    assert len(edges) == 2

def test_return_edge_data(mock_graph):
    mock_graph.define_edges()
    edge_data = mock_graph.return_edge_data(attribute='bo')
    assert np.array_equal(edge_data, np.array([1, 1]))

def test_return_visible_edges(mock_graph):
    mock_graph.define_edges()
    edges = mock_graph.return_visible_edges()
    assert len(edges) == 2

def test_return_unvisible_edges(mock_graph):
    mock_graph.define_edges()
    for edge in mock_graph.edges(data=True):
        edge[2]['visible'] = False
    edges = mock_graph.return_unvisible_edges()
    assert len(edges) == 2

def test_set_edges_visible(mock_graph):
    mock_graph.define_edges()
    mock_graph.set_edges_visible(attributes='bo', values=1)
    for edge in mock_graph.edges(data=True):
        assert edge[2]['visible'] == True

def test_set_edges_unvisible(mock_graph):
    mock_graph.define_edges()
    mock_graph.set_edges_unvisible(attributes='bo', values=1)
    for edge in mock_graph.edges(data=True):
        assert edge[2]['visible'] == False

def test_get_distance_matrix(mock_graph):
    mock_graph.define_edges()
    distance_matrix = mock_graph.get_graph_distance_matrix()
    assert distance_matrix.shape == (3, 3)

def test_get_all_shortest_paths(mock_graph):
    mock_graph.define_edges()
    paths = mock_graph.get_all_shortest_paths()
    assert len(paths) == 3
    for source in paths:
        assert len(paths[source]) == 2

def test_get_shortest_path(mock_graph):
    mock_graph.define_edges()
    path = mock_graph.get_shortest_path(0, 2)
    assert path == [0, 1, 2]