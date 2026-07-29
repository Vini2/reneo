import sys
import types

import pytest


class FakeDiGraph:
    def __init__(self):
        self._nodes = set()
        self._edges = {}

    def add_nodes_from(self, nodes):
        self._nodes.update(nodes)

    def add_edge(self, source, target, **attrs):
        self._nodes.update([source, target])
        self._edges[(source, target)] = attrs

    @property
    def nodes(self):
        return list(self._nodes)

    def out_edges(self, node):
        return [(source, target) for source, target in self._edges if source == node]

    def in_edges(self, node):
        return [(source, target) for source, target in self._edges if target == node]

    def out_degree(self, node):
        return len(self.out_edges(node))

    def in_degree(self, node):
        return len(self.in_edges(node))


@pytest.fixture
def fd_inexact(monkeypatch):
    import reneo_utils

    fake_networkx = types.ModuleType("networkx")
    fake_networkx.DiGraph = FakeDiGraph
    fake_flowpaths = types.ModuleType("flowpaths")
    fake_flowpaths.AbstractPathModelDAG = object
    fake_flowpaths.stDAG = lambda graph: graph

    monkeypatch.setitem(sys.modules, "networkx", fake_networkx)
    monkeypatch.setitem(sys.modules, "flowpaths", fake_flowpaths)
    sys.modules.pop("reneo_utils.fd_inexact", None)
    if hasattr(reneo_utils, "fd_inexact"):
        delattr(reneo_utils, "fd_inexact")
    import reneo_utils.fd_inexact as module

    yield module
    sys.modules.pop("reneo_utils.fd_inexact", None)
    if hasattr(reneo_utils, "fd_inexact"):
        delattr(reneo_utils, "fd_inexact")


def test_read_input_parses_graphs_and_subpaths(tmp_path, fd_inexact):
    graphfile = tmp_path / "graph.txt"
    graphfile.write_text(
        "\n".join(
            [
                "# graph",
                "3",
                "s a 1 2",
                "subpaths",
                "s a t ",
                "# filler",
                "# filler",
                "# filler",
                "# next",
            ]
        )
    )

    graphs = fd_inexact.read_input(str(graphfile), number_subpath=1)

    assert graphs[0]["Nodes"] == 3
    assert graphs[0]["list of edges"] == [("s", "a", "1", "2")]
    assert graphs[0]["subpaths"] == {0: ["s", "a", "t"]}


def test_get_subpath_constraints_keeps_only_existing_edges(fd_inexact):
    constraints = fd_inexact._get_subpath_constraints(
        {0: ["s", "a", "t"], 1: ["x", "y"]},
        {("s", "a"), ("a", "t")},
        lambda node: f"n_{node}",
    )

    assert constraints == [[("n_s", "n_a"), ("n_a", "n_t")]]


def test_path_to_edges_restores_original_labels(fd_inexact):
    edges = fd_inexact._path_to_edges(["1", "2", "3"], int)

    assert edges == [(1, 2), (2, 3)]


def test_fd_algorithm_returns_first_solved_decomposition(monkeypatch, fd_inexact):
    calls = []

    def fake_flow_multiple(data, k, nthreads):
        calls.append(k)
        data = data.copy()
        if k == 2:
            data["message"] = "solved"
            data["solution"] = [[("s", "t")]]
            data["weights"] = [7]
        else:
            data["message"] = "unsolved"
        return data

    monkeypatch.setattr(fd_inexact, "flowMultipleDecomposition", fake_flow_multiple)

    _, solution_paths = fd_inexact.FD_Algorithm({"message": {}}, max_paths=3, nthreads=1)

    assert calls == [1, 2]
    assert solution_paths == {0: {"weight": 7, "path": [("s", "t")]}}
