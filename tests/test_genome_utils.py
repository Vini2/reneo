from reneo_utils.genome_utils import GenomeComponent, GenomePath


def test_genome_path_stores_constructor_values():
    path = GenomePath(
        id="path_1",
        bubble_case="case1",
        node_order=["edge_1", "edge_2"],
        node_id_order=[1, 2],
        path="ACGT",
        coverage=12.5,
        length=4,
        gc=0.5,
    )

    assert path.id == "path_1"
    assert path.bubble_case == "case1"
    assert path.node_order == ["edge_1", "edge_2"]
    assert path.node_id_order == [1, 2]
    assert path.path == "ACGT"
    assert path.coverage == 12.5
    assert path.length == 4
    assert path.gc == 0.5


def test_genome_component_stores_summary_metrics():
    component = GenomeComponent(
        id=3,
        n_nodes=4,
        n_paths=2,
        max_degree=5,
        min_degree=1,
        max_in_degree=3,
        max_out_degree=4,
        avg_degree=2.5,
        avg_in_degree=1.25,
        avg_out_degree=1.25,
        density=0.75,
        max_path_length=1000,
        min_path_length=250,
        min_max_len_ratio=4.0,
        max_cov_path_length=900,
        min_cov_path_length=300,
        min_max_cov_len_ratio=3.0,
        max_cov=42.0,
        min_cov=6.0,
        min_max_cov_ratio=7.0,
        frac_unitigs=0.8,
    )

    assert component.id == 3
    assert component.n_nodes == 4
    assert component.n_paths == 2
    assert component.max_degree == 5
    assert component.min_degree == 1
    assert component.max_in_degree == 3
    assert component.max_out_degree == 4
    assert component.avg_degree == 2.5
    assert component.avg_in_degree == 1.25
    assert component.avg_out_degree == 1.25
    assert component.density == 0.75
    assert component.max_path_length == 1000
    assert component.min_path_length == 250
    assert component.min_max_len_ratio == 4.0
    assert component.max_cov_path_length == 900
    assert component.min_cov_path_length == 300
    assert component.min_max_cov_len_ratio == 3.0
    assert component.max_cov == 42.0
    assert component.min_cov == 6.0
    assert component.min_max_cov_ratio == 7.0
    assert component.frac_unitigs == 0.8
