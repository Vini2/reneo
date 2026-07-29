from pathlib import Path

from reneo_utils.gene_utils import get_smg_unitigs, get_vog_unitigs


def test_get_smg_unitigs_filters_by_mapped_marker_fraction(tmp_path):
    hmmout = tmp_path / "markers.domtblout"
    hmmout.write_text(
        "\n".join(
            [
                "# ignored comment",
                "edge_1_10_20_+ unused unused markerA unused 100 unused unused unused unused unused unused unused unused unused 1 90",
                "edge_2_10_20_- unused unused markerB unused 100 unused unused unused unused unused unused unused unused unused 1 20",
            ]
        )
        + "\n"
    )

    assert get_smg_unitigs(str(hmmout), mg_frac=0.5) == {"edge_1"}


def test_get_vog_unitigs_filters_thresholds_and_hypothetical_annotations(tmp_path):
    annotations = tmp_path / "vog.annotations.tsv"
    annotations.write_text(
        "\n".join(
            [
                "# comment",
                "VOG0001\t-\t-\t-\tmajor capsid protein",
                "VOG0002\t-\t-\t-\thypothetical protein",
                "VOG0003\t-\t-\t-\tDNA polymerase",
            ]
        )
        + "\n"
    )
    vogs = tmp_path / "all.hmmVOG.tbl"
    vogs.write_text(
        "\n".join(
            [
                "# comment",
                "edge_1_gene1 x VOG0001 x 1e-20 80",
                "edge_1_gene2 x VOG0002 x 1e-20 80",
                "edge_2_gene1 x VOG0003 x 1e-2 80",
                "edge_3_gene1 x VOG0003 x 1e-20 20",
            ]
        )
        + "\n"
    )

    unitig_vogs, vog_dict = get_vog_unitigs(
        str(vogs), e_value=1e-10, hmm_score=50, vogfunctions=str(annotations)
    )

    assert unitig_vogs == {"edge_1": {"VOG0001"}}
    assert vog_dict["VOG0001"] == "major capsid protein"
