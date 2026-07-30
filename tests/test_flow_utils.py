import sys
import types

import pytest


class FakeNetworkX(types.ModuleType):
    @staticmethod
    def bfs_layers(graph, source):
        layers = [[source]]
        seen = {source}
        frontier = [source]
        while frontier:
            next_frontier = []
            for node in frontier:
                for successor in graph.successors(node):
                    if successor not in seen:
                        seen.add(successor)
                        next_frontier.append(successor)
            if next_frontier:
                layers.append(next_frontier)
            frontier = next_frontier
        return layers


@pytest.fixture
def flow_utils(monkeypatch):
    import reneo_utils

    fake_fd = types.ModuleType("reneo_utils.fd_inexact")
    fake_fd.SolveInstances = lambda graphs, max_paths, outputfile, recordfile, nthreads: {
        0: {"weight": 1, "path": [("s", "t")]}
    }
    monkeypatch.setitem(sys.modules, "networkx", FakeNetworkX("networkx"))
    monkeypatch.setitem(sys.modules, "reneo_utils.fd_inexact", fake_fd)
    sys.modules.pop("reneo_utils.flow_utils", None)
    if hasattr(reneo_utils, "flow_utils"):
        delattr(reneo_utils, "flow_utils")
    import reneo_utils.flow_utils as module

    yield module
    sys.modules.pop("reneo_utils.flow_utils", None)
    if hasattr(reneo_utils, "flow_utils"):
        delattr(reneo_utils, "flow_utils")


class DiGraph:
    def __init__(self, edges):
        self.edges = list(edges)
        self.nodes = sorted({node for edge in edges for node in edge})

    def predecessors(self, node):
        return [source for source, target in self.edges if target == node]

    def successors(self, node):
        return [target for source, target in self.edges if source == node]


def test_get_source_sink_linear_excludes_self_looped_unitigs(flow_utils):
    graph = DiGraph(
        [
            ("edge_1+", "edge_2+"),
            ("edge_3+", "edge_3-"),
            ("edge_4+", "edge_5+"),
        ]
    )

    sources, sinks = flow_utils.get_source_sink_linear(
        graph,
        graph_unitigs={
            "edge_1": "A" * 100,
            "edge_2": "A" * 100,
            "edge_3": "A" * 100,
            "edge_4": "A" * 1001,
            "edge_5": "A" * 1001,
        },
        self_looped_nodes={"edge_3"},
    )

    assert sources == ["edge_4+"]
    assert sinks == ["edge_5+"]


def test_get_source_sink_circular_finds_node_that_returns_to_start(flow_utils):
    graph = DiGraph([("edge_1+", "edge_2+"), ("edge_2+", "edge_1+")])

    candidates = flow_utils.get_source_sink_circular(
        graph,
        graph_unitigs={"edge_1": "A" * 1001, "edge_2": "A" * 1001},
        self_looped_nodes=set(),
    )

    assert candidates == ["edge_1+", "edge_2+"]


def test_solve_mfd_delegates_to_solve_instances(flow_utils):
    graph = {"list of edges": [("s", "t", "1", "2")], "subpaths": {}}

    assert flow_utils.solve_mfd(graph, max_paths=2, nthreads=1) == {
        0: {"weight": 1, "path": [("s", "t")]}
    }
