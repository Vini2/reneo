import sys
import types

import pytest


def install_import_stubs(monkeypatch):
    bio = types.ModuleType("Bio")
    seqio = types.ModuleType("Bio.SeqIO")
    seqio.parse = lambda *args, **kwargs: ()
    bio.SeqIO = seqio

    agtools = types.ModuleType("agtools")
    core = types.ModuleType("agtools.core")
    unitig_graph = types.ModuleType("agtools.core.unitig_graph")
    unitig_graph.UnitigGraph = object

    monkeypatch.setitem(sys.modules, "Bio", bio)
    monkeypatch.setitem(sys.modules, "Bio.SeqIO", seqio)
    monkeypatch.setitem(sys.modules, "agtools", agtools)
    monkeypatch.setitem(sys.modules, "agtools.core", core)
    monkeypatch.setitem(sys.modules, "agtools.core.unitig_graph", unitig_graph)
    monkeypatch.setitem(sys.modules, "igraph", types.ModuleType("igraph"))


@pytest.fixture
def edge_graph_utils(monkeypatch):
    install_import_stubs(monkeypatch)
    import reneo_utils

    sys.modules.pop("reneo_utils.edge_graph_utils", None)
    if hasattr(reneo_utils, "edge_graph_utils"):
        delattr(reneo_utils, "edge_graph_utils")
    import reneo_utils.edge_graph_utils as module

    yield module
    sys.modules.pop("reneo_utils.edge_graph_utils", None)
    if hasattr(reneo_utils, "edge_graph_utils"):
        delattr(reneo_utils, "edge_graph_utils")


class DiGraph:
    def __init__(self, edges):
        self.edges = list(edges)
        self.nodes = sorted({node for edge in edges for node in edge})
        self.removed = []

    def in_degree(self, node):
        return sum(1 for _, target in self.edges if target == node)

    def out_degree(self):
        return lambda node: sum(1 for source, _ in self.edges if source == node)

    def remove_nodes_from(self, nodes):
        self.removed.extend(nodes)
        self.nodes = [node for node in self.nodes if node not in nodes]
        self.edges = [
            (source, target)
            for source, target in self.edges
            if source not in nodes and target not in nodes
        ]


def test_bidirectional_map_maintains_inverse(edge_graph_utils):
    mapping = edge_graph_utils.BidirectionalMap()
    mapping[1] = "edge_1"

    assert mapping[1] == "edge_1"
    assert mapping.inverse["edge_1"] == 1

    del mapping[1]
    assert "edge_1" not in mapping.inverse


def test_bidirectional_map_rejects_duplicate_values(edge_graph_utils):
    mapping = edge_graph_utils.BidirectionalMap()
    mapping[1] = "edge_1"

    with pytest.raises(edge_graph_utils.BidirectionalError):
        mapping[2] = "edge_1"


def test_get_name_based_links_translates_segment_ids(edge_graph_utils):
    graph = types.SimpleNamespace(
        oriented_links={0: {1: {("+", "-"), ("-", "+")}}},
        link_overlap={(0, "+", 1, "-"): 31},
    )
    names = {0: "edge_1", 1: "edge_2"}

    oriented_links, link_overlap = edge_graph_utils.get_name_based_links(graph, names)

    assert oriented_links["edge_1"]["edge_2"] == [("+", "-"), ("-", "+")]
    assert link_overlap[("edge_1+", "edge_2-")] == 31


def test_get_circular_uses_self_loop_lengths(edge_graph_utils):
    circular = edge_graph_utils.get_circular(
        ["edge_1", "edge_2"], {"edge_1": "ACGT", "edge_2": "AA"}
    )

    assert circular == {"edge_1": 4, "edge_2": 2}


def test_remove_dead_ends_returns_iteratively_pruned_nodes(edge_graph_utils):
    graph = DiGraph([("source", "middle"), ("middle", "sink"), ("cycle", "cycle")])

    assert edge_graph_utils.remove_dead_ends(graph) == {"source", "middle", "sink"}
