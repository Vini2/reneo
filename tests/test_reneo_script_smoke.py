import ast
import csv
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace


def read_resolved_genome_paths(genome_info):
    with genome_info.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {
            row["Path"]: {
                "case": row["Case"],
                "length": int(row["Length"]),
                "node_order": ast.literal_eval(row["Node order"]),
            }
            for row in rows
        }


def read_fasta_ids(fasta):
    with fasta.open() as handle:
        return [line[1:].strip() for line in handle if line.startswith(">")]


def path_signature(paths):
    return {
        (
            path["case"],
            path["length"],
            tuple(path["node_order"]),
        )
        for path in paths.values()
    }


def test_reneo_script_smoke_with_mock_intermediates(tmp_path, monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "tests" / "data"
    script_dir = repo_root / "reneo" / "workflow" / "scripts"
    output_dir = tmp_path / "reneo_output"
    output_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "matplotlib"))

    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            graph=str(data_dir / "assembly_graph_after_simplification.gfa"),
            coverage=str(data_dir / "reneo.coverage.tsv"),
            pickle=str(data_dir / "PE_junctions.pkl"),
        ),
        params=SimpleNamespace(
            bampath=str(data_dir),
            genomes_folder=None,
            unitigs=None,
            hmmout=None,
            vogs=str(data_dir / "all.hmmVOG.tbl"),
            vogfunctions=str(data_dir / "vog.annotations.tsv"),
            minlength=5000,
            mincov=50,
            compcount=200,
            maxpaths=10,
            mgfrac=0.2,
            evalue=1e-10,
            hmmscore=50,
            nvogs=10,
            covtol=100,
            alpha=1.1,
            output=str(output_dir),
        ),
        threads=1,
        log=SimpleNamespace(stderr=str(tmp_path / "reneo_output.err")),
    )

    sys.path.insert(0, str(script_dir))
    try:
        runpy.run_path(
            str(script_dir / "reneo.py"),
            run_name="__main__",
            init_globals={"snakemake": snakemake},
        )
    finally:
        sys.path.remove(str(script_dir))

    expected_outputs = [
        output_dir / "resolved_paths.fasta",
        output_dir / "resolved_genome_info.txt",
        output_dir / "resolved_component_info.txt",
        output_dir / "component_vogs.txt",
        output_dir / "resolved_edges.fasta",
        output_dir / "unresolved_virus_like_edges.fasta",
        output_dir / "assembly_graph_after_simplification.augmented.gfa",
        output_dir / "assembly_graph_after_simplification.augmented.summary.tsv",
    ]
    for output in expected_outputs:
        assert output.exists(), f"Missing expected Reneo output: {output}"

    log_text = (tmp_path / "reneo_output.err").read_text()
    assert "Minimum coverage of paths to output: 50" in log_text
    assert "Coverage tolerance for extending subpaths: 100.0" in log_text
    assert "Minimum length of unitigs to consider: 5000" in log_text
    assert "Coverage multipler for flow interval modelling: 1.1" in log_text

    expected_paths = read_resolved_genome_paths(data_dir / "resolved_genome_info.txt")
    observed_paths = read_resolved_genome_paths(output_dir / "resolved_genome_info.txt")
    assert path_signature(observed_paths) == path_signature(expected_paths)
    assert read_fasta_ids(output_dir / "resolved_paths.fasta") == list(observed_paths)
