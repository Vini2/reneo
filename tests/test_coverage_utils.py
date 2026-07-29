from collections import defaultdict
import importlib.util
import sys
import types
from types import SimpleNamespace

import pytest


@pytest.fixture
def coverage_utils(monkeypatch):
    import reneo_utils

    if importlib.util.find_spec("pysam") is None:
        monkeypatch.setitem(sys.modules, "pysam", types.ModuleType("pysam"))

    sys.modules.pop("reneo_utils.coverage_utils", None)
    if hasattr(reneo_utils, "coverage_utils"):
        delattr(reneo_utils, "coverage_utils")
    import reneo_utils.coverage_utils as module

    yield module
    sys.modules.pop("reneo_utils.coverage_utils", None)
    if hasattr(reneo_utils, "coverage_utils"):
        delattr(reneo_utils, "coverage_utils")


def fake_read(reference_name="edge_1", is_reverse=False, tags=None, **overrides):
    values = {
        "reference_name": reference_name,
        "is_reverse": is_reverse,
        "reference_start": 0,
        "query_alignment_start": 0,
        "query_alignment_end": 10,
    }
    values.update(overrides)
    tags = tags or {}

    read = SimpleNamespace(**values)
    read.has_tag = lambda tag: tag in tags
    read.get_tag = lambda tag: tags[tag]
    return read


def test_get_unitig_coverage_sums_sample_columns(tmp_path, coverage_utils):
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text("Contig sample_a sample_b\nedge_1 1.5 2.5\nedge_2 0 3\n")

    assert coverage_utils.get_unitig_coverage(str(coverage)) == {
        "edge_1": 4.0,
        "edge_2": 3.0,
    }


def test_orientation_helpers_use_reference_strand(coverage_utils):
    assert coverage_utils.get_read_orientation(fake_read(is_reverse=False)) == "+"
    assert coverage_utils.get_read_orientation(fake_read(is_reverse=True)) == "-"
    assert coverage_utils.get_opposite_orientation("+") == "-"
    assert coverage_utils.get_opposite_orientation("-") == "+"


def test_add_oriented_pair_support_counts_link_and_reverse_complement(coverage_utils):
    counts = defaultdict(int)

    coverage_utils.add_oriented_pair_support_from_fields(
        counts, "edge_1", "+", "edge_2", "-"
    )

    assert counts[("edge_1+", "edge_2+")] == 1
    assert counts[("edge_2-", "edge_1-")] == 1


def test_add_oriented_spanning_read_support_counts_link_and_reverse_complement(
    coverage_utils,
):
    counts = defaultdict(int)
    left = fake_read(reference_name="edge_1", is_reverse=False)
    right = fake_read(reference_name="edge_2", is_reverse=True)

    coverage_utils.add_oriented_spanning_read_support(counts, left, right)

    assert counts[("edge_1+", "edge_2-")] == 1
    assert counts[("edge_2+", "edge_1-")] == 1


def test_endpoint_pair_support_only_counts_target_contigs(coverage_utils):
    counts = defaultdict(int)

    coverage_utils.add_endpoint_pair_support_from_fields(
        counts, "edge_1", "+", "edge_2", "-", target_contigs={"edge_1"}
    )

    assert counts == {"edge_1+": 1, "edge_1-": 1}


def test_endpoint_spanning_read_support_counts_target_ends(coverage_utils):
    counts = defaultdict(int)
    left = fake_read(reference_name="edge_1", is_reverse=False)
    right = fake_read(reference_name="edge_2", is_reverse=True)

    coverage_utils.add_endpoint_spanning_read_support(
        counts, left, right, target_contigs={"edge_2"}
    )

    assert counts == {"edge_2-": 1, "edge_2+": 1}


def test_endpoint_sa_tag_support_ignores_same_contig_and_malformed_alignments(
    coverage_utils,
):
    counts = defaultdict(int)
    read = fake_read(
        reference_name="edge_1",
        is_reverse=False,
        tags={"SA": "edge_1,1,+,10M,60,0;edge_2,5,-,10M,60,0;bad;"},
    )

    coverage_utils.add_endpoint_sa_tag_support(
        counts, read, target_contigs={"edge_1", "edge_2"}
    )

    assert counts == {"edge_1+": 1, "edge_2-": 1, "edge_2+": 1, "edge_1-": 1}
