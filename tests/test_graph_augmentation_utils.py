import sys
import types
from collections import defaultdict

import pytest


class FakeNoCycle(Exception):
    pass


class FakeNetworkX(types.ModuleType):
    class exception:
        NetworkXNoCycle = FakeNoCycle

    class DiGraph:
        def __init__(self):
            self.edges = {}
            self.nodes = set()

        def add_edge(self, source, target, weight=None):
            self.edges[(source, target)] = {"weight": weight}
            self.nodes.update([source, target])

        def predecessors(self, node):
            return [source for source, target in self.edges if target == node]

        def successors(self, node):
            return [target for source, target in self.edges if source == node]

    @staticmethod
    def find_cycle(graph, orientation=None):
        raise FakeNoCycle


@pytest.fixture
def graph_augmentation_utils(monkeypatch):
    import reneo_utils

    fake_coverage = types.ModuleType("reneo_utils.coverage_utils")
    fake_coverage.get_opposite_orientation = lambda orientation: (
        "-" if orientation == "+" else "+"
    )
    fake_coverage.get_oriented_external_endpoint_read_support = (
        lambda *args, **kwargs: defaultdict(int)
    )
    fake_coverage.get_component_external_endpoint_read_support = (
        lambda *args, **kwargs: defaultdict(int)
    )
    fake_coverage.get_oriented_junction_pe_coverage_for_pairs = (
        lambda *args, **kwargs: defaultdict(int)
    )
    fake_coverage.get_oriented_spanning_read_coverage_for_pairs = (
        lambda *args, **kwargs: defaultdict(int)
    )
    fake_component = types.ModuleType("reneo_utils.component_utils")
    fake_component.get_components = lambda **kwargs: (kwargs["pruned_vs"], kwargs["comp_vogs"])
    fake_flow = types.ModuleType("reneo_utils.flow_utils")
    fake_flow.get_source_sink_linear = lambda graph, self_looped_nodes: (
        [
            node
            for node in graph.nodes
            if len(list(graph.predecessors(node))) == 0
            and len(list(graph.successors(node))) > 0
        ],
        [
            node
            for node in graph.nodes
            if len(list(graph.predecessors(node))) > 0
            and len(list(graph.successors(node))) == 0
        ],
    )

    monkeypatch.setitem(sys.modules, "networkx", FakeNetworkX("networkx"))
    monkeypatch.setitem(sys.modules, "reneo_utils.coverage_utils", fake_coverage)
    monkeypatch.setitem(sys.modules, "reneo_utils.component_utils", fake_component)
    monkeypatch.setitem(sys.modules, "reneo_utils.flow_utils", fake_flow)
    sys.modules.pop("reneo_utils.graph_augmentation_utils", None)
    if hasattr(reneo_utils, "graph_augmentation_utils"):
        delattr(reneo_utils, "graph_augmentation_utils")
    import reneo_utils.graph_augmentation_utils as module

    yield module
    sys.modules.pop("reneo_utils.graph_augmentation_utils", None)
    if hasattr(reneo_utils, "graph_augmentation_utils"):
        delattr(reneo_utils, "graph_augmentation_utils")


class AssemblyGraph:
    def __init__(self, components=None):
        self._components = components or []
        self.edges = set()

    def components(self):
        return self._components

    def are_connected(self, left, right):
        return tuple(sorted([left, right])) in self.edges

    def add_edge(self, left, right):
        self.edges.add(tuple(sorted([left, right])))


class Logger:
    def info(self, message):
        pass


def test_build_oriented_component_graph_uses_min_nonzero_coverage(graph_augmentation_utils):
    graph = graph_augmentation_utils.build_oriented_component_graph(
        [0, 1],
        unitig_names={0: "edge_1", 1: "edge_2"},
        unitig_names_rev={"edge_1": 0, "edge_2": 1},
        oriented_links={"edge_1": {"edge_2": [("+", "-")]}, "edge_2": {}},
        self_looped_nodes=set(),
        unitig_coverages={"edge_1": 0, "edge_2": 7},
        MAX_VAL=999,
    )

    assert graph.edges[("edge_1+", "edge_2-")] == {"weight": 7}


def test_candidate_pair_helpers_filter_by_threshold(graph_augmentation_utils):
    junctions = {
        ("terminal", "isolated"): 10,
        ("terminal", "weak"): 9,
        ("terminal", "other_terminal"): 12,
        ("x", "y", "z"): 20,
    }

    assert graph_augmentation_utils.filter_isolated_unitigs_by_junction_support(
        {"isolated", "weak"}, {"terminal"}, junction_pe_coverage=junctions
    ) == {"isolated"}
    assert graph_augmentation_utils.get_candidate_extension_pairs(
        {"isolated", "weak"},
        {"terminal", "other_terminal"},
        junction_pe_coverage=junctions,
    ) == {("isolated", "terminal"), ("other_terminal", "terminal")}


def test_is_viral_singleton_candidate_checks_all_requirements(graph_augmentation_utils):
    kwargs = {
        "smg_unitigs": set(),
        "unitig_vogs": {"edge_1": {"VOG1", "VOG2"}},
        "nvogs": 2,
        "edges_lengths": {"edge_1": 6000},
        "minlength": 5000,
    }

    assert graph_augmentation_utils.is_viral_singleton_candidate("edge_1", **kwargs)
    assert not graph_augmentation_utils.is_viral_singleton_candidate(
        "edge_2", **kwargs
    )


def test_add_inferred_oriented_link_adds_primary_reverse_and_graph_edge(graph_augmentation_utils):
    oriented_links = defaultdict(lambda: defaultdict(list))
    junction_pe_coverage = defaultdict(int)
    assembly_graph = AssemblyGraph()

    added = graph_augmentation_utils.add_inferred_oriented_link(
        "edge_1+",
        "edge_2-",
        15,
        oriented_links=oriented_links,
        link_overlap={},
        junction_pe_coverage=junction_pe_coverage,
        unitig_names_rev={"edge_1": 1, "edge_2": 2},
        assembly_graph=assembly_graph,
    )

    assert added is True
    assert oriented_links["edge_1"]["edge_2"] == [("+", "-")]
    assert oriented_links["edge_2"]["edge_1"] == [("+", "-")]
    assert junction_pe_coverage[("edge_1", "edge_2")] == 15
    assert junction_pe_coverage[("edge_2", "edge_1")] == 15
    assert assembly_graph.edges == {(1, 2)}


def test_add_complete_isolated_linear_unitigs_adds_no_support_candidates(
    monkeypatch, graph_augmentation_utils
):
    monkeypatch.setattr(
        graph_augmentation_utils,
        "get_oriented_external_endpoint_read_support",
        lambda *args, **kwargs: defaultdict(int),
    )
    pruned_vs = {}
    comp_vogs = {}

    count = graph_augmentation_utils.add_complete_isolated_linear_unitigs(
        assembly_graph=AssemblyGraph(components=[[0]]),
        unitig_names={0: "edge_1"},
        unitig_names_rev={"edge_1": 0},
        self_looped_nodes=set(),
        smg_unitigs=set(),
        unitig_vogs={"edge_1": {"VOG1"}},
        nvogs=1,
        edges_lengths={"edge_1": 6000},
        minlength=5000,
        inferred_case3_links={},
        pruned_vs=pruned_vs,
        comp_vogs=comp_vogs,
        bampath=".",
        output=".",
        nthreads=1,
        logger=Logger(),
    )

    assert count == 1
    assert pruned_vs == {0: [0]}
    assert comp_vogs == {0: {"VOG1"}}


def test_filter_unextended_case3_linear_components_by_endpoint_support(
    monkeypatch, graph_augmentation_utils
):
    endpoint_support = {"edge_1+": 3, "edge_6+": 0}
    monkeypatch.setattr(
        graph_augmentation_utils,
        "get_component_external_endpoint_read_support",
        lambda *args, **kwargs: endpoint_support,
    )
    pruned_vs = {5: [0, 1, 2], 6: [3, 4, 5]}
    comp_vogs = {5: {"VOG1"}, 6: {"VOG2"}}

    removed = graph_augmentation_utils.filter_unextended_case3_linear_components_by_endpoint_support(
        pruned_vs=pruned_vs,
        comp_vogs=comp_vogs,
        unitig_names={
            0: "edge_1",
            1: "edge_2",
            2: "edge_3",
            3: "edge_4",
            4: "edge_5",
            5: "edge_6",
        },
        unitig_names_rev={
            "edge_1": 0,
            "edge_2": 1,
            "edge_3": 2,
            "edge_4": 3,
            "edge_5": 4,
            "edge_6": 5,
        },
        oriented_links={
            "edge_1": {"edge_2": [("+", "+")]},
            "edge_2": {"edge_3": [("+", "+")]},
            "edge_3": {},
            "edge_4": {"edge_5": [("+", "+")]},
            "edge_5": {"edge_6": [("+", "+")]},
            "edge_6": {},
        },
        self_looped_nodes=set(),
        unitig_coverages={},
        MAX_VAL=999,
        compcount=200,
        inferred_case3_links={},
        bampath=".",
        output=".",
        nthreads=1,
        logger=Logger(),
    )

    assert removed == 1
    assert pruned_vs == {6: [3, 4, 5]}
    assert comp_vogs == {6: {"VOG2"}}
