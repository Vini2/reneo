from types import SimpleNamespace

from reneo_utils.output_utils import (
    get_augmented_gfa_path,
    get_augmented_summary_path,
    init_files,
    write_augmented_gfa,
    write_augmented_summary,
    write_component_info,
    write_component_vog_info,
    write_path,
    write_path_fasta,
    write_res_genome_info,
    write_unitigs,
)


class Logger:
    def __init__(self):
        self.messages = []

    def info(self, message):
        self.messages.append(message)


def test_augmented_output_paths_strip_gfa_suffix(tmp_path):
    graph = tmp_path / "assembly.gfa"

    assert get_augmented_gfa_path(graph=str(graph), output=str(tmp_path)) == (
        f"{tmp_path}/assembly.augmented.gfa"
    )
    assert get_augmented_summary_path(graph=str(graph), output=str(tmp_path)) == (
        f"{tmp_path}/assembly.augmented.summary.tsv"
    )


def test_write_augmented_gfa_appends_sorted_inferred_links(tmp_path):
    graph = tmp_path / "graph.gfa"
    graph.write_text("S\tedge_1\tACGT")
    logger = Logger()

    output = write_augmented_gfa(
        graph=str(graph),
        output=str(tmp_path),
        logger=logger,
        inferred_case3_links={
            ("edge_2-", "edge_3+"): {"support": 12, "round": 2},
            ("edge_1+", "edge_2-"): {
                "support": 20,
                "pe_support": 15,
                "spanning_support": 5,
                "round": 1,
            },
        },
    )

    assert (tmp_path / "graph.augmented.gfa").read_text().splitlines() == [
        "S\tedge_1\tACGT",
        "L\tedge_1\t+\tedge_2\t-\t0M\tRC:i:20\tPE:i:15\tSR:i:5\tIR:i:1\tRN:Z:case3_extension",
        "L\tedge_2\t-\tedge_3\t+\t0M\tRC:i:12\tPE:i:0\tSR:i:0\tIR:i:2\tRN:Z:case3_extension",
    ]
    assert output == f"{tmp_path}/graph.augmented.gfa"
    assert logger.messages


def test_write_augmented_summary_includes_metadata_defaults(tmp_path):
    graph = tmp_path / "graph.gfa"
    graph.write_text("S\tedge_1\tA\n")

    write_augmented_summary(
        graph=str(graph),
        output=str(tmp_path),
        logger=Logger(),
        inferred_case3_links={
            ("edge_1+", "edge_2-"): {"support": 11, "round": 3},
        },
    )

    assert (tmp_path / "graph.augmented.summary.tsv").read_text().splitlines() == [
        "round\tleft_oriented_contig\tright_oriented_contig\tleft_contig\tright_contig\tlink_type\tleft_endpoint_type\tright_endpoint_type\tleft_component\tright_component\tsupport\tpe_support\tspanning_support",
        "3\tedge_1+\tedge_2-\tedge_1\tedge_2\tunknown\tunknown\tunknown\tNA\tNA\t11\t0\t0",
    ]


def test_write_sequence_outputs_wrap_fasta_lines(tmp_path):
    long_sequence = "A" * 61

    write_unitigs(
        nodes=[0],
        unitig_names={0: "edge_1"},
        graph_unitigs={"edge_1": long_sequence},
        filename="unitigs",
        output=str(tmp_path),
    )
    write_path(
        [SimpleNamespace(id="path_1", path=long_sequence, length=len(long_sequence))],
        str(tmp_path),
    )
    write_path_fasta(
        [SimpleNamespace(id="path_1", path=long_sequence, length=len(long_sequence))],
        str(tmp_path / "genomes"),
    )

    assert (tmp_path / "unitigs.fasta").read_text().splitlines() == [
        ">edge_1",
        "A" * 60,
        "A",
    ]
    assert (tmp_path / "resolved_paths.fasta").read_text().splitlines() == [
        ">path_1",
        "A" * 60,
        "A",
    ]
    assert (tmp_path / "genomes" / "path_1.fasta").read_text().splitlines() == [
        ">path_1",
        "A" * 60,
        "A",
    ]


def test_write_info_files(tmp_path):
    component = SimpleNamespace(
        id=1,
        n_nodes=2,
        n_paths=1,
        frac_unitigs=1.0,
        max_degree=2,
        min_degree=1,
        max_in_degree=1,
        max_out_degree=1,
        avg_degree=1.5,
        avg_in_degree=0.75,
        avg_out_degree=0.75,
        density=0.5,
        max_path_length=100,
        min_path_length=80,
        min_max_len_ratio=1.25,
        max_cov_path_length=100,
        min_cov_path_length=80,
        min_max_cov_len_ratio=1.25,
        max_cov=20,
        min_cov=10,
        min_max_cov_ratio=2,
    )
    genome_path = SimpleNamespace(
        id="path_1",
        bubble_case="case1",
        coverage=20,
        length=100,
        gc=0.5,
        node_order=["edge_1"],
    )

    assert write_component_info([component], str(tmp_path)) == "resolved_component_info.txt"
    assert write_res_genome_info([genome_path], str(tmp_path)) == "resolved_genome_info.txt"
    assert write_component_vog_info([1], {1: {"VOG1"}}, str(tmp_path)) == "component_vogs.txt"
    init_files(str(tmp_path))

    assert "Component\tNumber of nodes" in (
        tmp_path / "resolved_component_info.txt"
    ).read_text()
    assert "path_1\tcase1\t20\t100\t0.5\t['edge_1']" in (
        tmp_path / "resolved_genome_info.txt"
    ).read_text()
    assert (tmp_path / "component_vogs.txt").read_text().splitlines() == [
        "Component\tvog",
        "1\t{'VOG1'}",
    ]
    assert (tmp_path / "resolved_edges.fasta").exists()
    assert (tmp_path / "resolved_paths.fasta").exists()
