#!/usr/bin/env python3

import os
import subprocess

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.6.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


FASTA_LINE_LEN = 60


def get_augmented_gfa_path(**kwargs):
    graph_basename = os.path.basename(kwargs["graph"])

    if graph_basename.endswith(".gfa"):
        graph_basename = graph_basename[:-4]

    return f"{kwargs['output']}/{graph_basename}.augmented.gfa"


def get_augmented_summary_path(**kwargs):
    graph_basename = os.path.basename(kwargs["graph"])

    if graph_basename.endswith(".gfa"):
        graph_basename = graph_basename[:-4]

    return f"{kwargs['output']}/{graph_basename}.augmented.summary.tsv"


def write_augmented_gfa(**kwargs):
    """
    Write the original GFA plus inferred case 3 extension links.
    """

    augmented_gfa = get_augmented_gfa_path(**kwargs)
    inferred_links = kwargs.get("inferred_case3_links", {})

    with open(kwargs["graph"]) as source, open(augmented_gfa, "w") as target:
        for line in source:
            target.write(line)
            if not line.endswith("\n"):
                target.write("\n")

        for (left, right), metadata in sorted(
            inferred_links.items(), key=lambda item: (item[1]["round"], item[0])
        ):
            target.write(
                f"L\t{left[:-1]}\t{left[-1]}\t{right[:-1]}\t{right[-1]}\t0M"
            )
            target.write(f"\tRC:i:{metadata['support']}")
            target.write(f"\tPE:i:{metadata.get('pe_support', 0)}")
            target.write(f"\tSR:i:{metadata.get('spanning_support', 0)}")
            target.write(f"\tIR:i:{metadata['round']}")
            target.write("\tRN:Z:case3_extension\n")

    kwargs["logger"].info(
        f"Augmented assembly graph written to {augmented_gfa} with {len(inferred_links)} inferred case 3 links"
    )

    return augmented_gfa


def write_augmented_summary(**kwargs):
    """
    Write a TSV summary of inferred case 3 extension links.
    """

    summary_file = get_augmented_summary_path(**kwargs)
    inferred_links = kwargs.get("inferred_case3_links", {})

    with open(summary_file, "w") as target:
        target.write(
            "round\tleft_oriented_contig\tright_oriented_contig\t"
            "left_contig\tright_contig\tlink_type\tleft_endpoint_type\t"
            "right_endpoint_type\tleft_component\tright_component\t"
            "support\tpe_support\tspanning_support\n"
        )

        for (left, right), metadata in sorted(
            inferred_links.items(), key=lambda item: (item[1]["round"], item[0])
        ):
            target.write(
                f"{metadata['round']}\t{left}\t{right}\t{left[:-1]}\t{right[:-1]}\t"
            )
            target.write(
                f"{metadata.get('link_type', 'unknown')}\t"
                f"{metadata.get('left_endpoint_type', 'unknown')}\t"
                f"{metadata.get('right_endpoint_type', 'unknown')}\t"
                f"{metadata.get('left_component', 'NA')}\t"
                f"{metadata.get('right_component', 'NA')}\t"
            )
            target.write(
                f"{metadata['support']}\t"
                f"{metadata.get('pe_support', 0)}\t"
                f"{metadata.get('spanning_support', 0)}\n"
            )

    kwargs["logger"].info(
        f"Case 3 graph augmentation summary written to {summary_file}"
    )

    return summary_file


def write_unitigs(nodes, unitig_names, graph_unitigs, filename, output):
    """
    Write unitigs to FASTA file
    """

    with open(f"{output}/{filename}.fasta", "w+") as myfile:
        for node in nodes:
            unitig_name = unitig_names[node]
            edge_seq = str(graph_unitigs[unitig_name])
            myfile.write(f">{unitig_name}\n")

            chunks = [
                edge_seq[i : i + FASTA_LINE_LEN]
                for i in range(0, len(edge_seq), FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_component_info(all_components, output):
    """
    Write component information to file
    """

    with open(f"{output}/resolved_component_info.txt", "w") as myfile:
        myfile.write(f"Component\t")
        myfile.write(f"Number of nodes\t")
        myfile.write(f"Number of paths\t")
        myfile.write(f"Fraction of unitigs recovered\t")
        myfile.write(f"Maximum degree\t")
        myfile.write(f"Minimum degree\t")
        myfile.write(f"Maximum in degree\t")
        myfile.write(f"Maximum out degree\t")
        myfile.write(f"Average degree\t")
        myfile.write(f"Average in degree\t")
        myfile.write(f"Average out degree\t")
        myfile.write(f"Density\t")
        myfile.write(f"Maximum path length\t")
        myfile.write(f"Minimum path length\t")
        myfile.write(f"Length ratio (long/short)\t")
        myfile.write(f"Maximum coverage path length\t")
        myfile.write(f"Minimum coverage path length\t")
        myfile.write(f"Length ratio (highest cov/lowest cov)\t")
        myfile.write(f"Maximum coverage\t")
        myfile.write(f"Minimum coverage\t")
        myfile.write(f"Coverage ratio (highest/lowest)\n")

        if len(all_components) > 0:
            for component in all_components:
                myfile.write(f"{component.id}\t")
                myfile.write(f"{component.n_nodes}\t")
                myfile.write(f"{component.n_paths}\t")
                myfile.write(f"{component.frac_unitigs}\t")
                myfile.write(f"{component.max_degree}\t")
                myfile.write(f"{component.min_degree}\t")
                myfile.write(f"{component.max_in_degree}\t")
                myfile.write(f"{component.max_out_degree}\t")
                myfile.write(f"{component.avg_degree}\t")
                myfile.write(f"{component.avg_in_degree}\t")
                myfile.write(f"{component.avg_out_degree}\t")
                myfile.write(f"{component.density}\t")
                myfile.write(f"{component.max_path_length}\t")
                myfile.write(f"{component.min_path_length}\t")
                myfile.write(f"{component.min_max_len_ratio}\t")
                myfile.write(f"{component.max_cov_path_length}\t")
                myfile.write(f"{component.min_cov_path_length}\t")
                myfile.write(f"{component.min_max_cov_len_ratio}\t")
                myfile.write(f"{component.max_cov}\t")
                myfile.write(f"{component.min_cov}\t")
                myfile.write(f"{component.min_max_cov_ratio}\n")
        else:
            myfile.write(f"No complex components were resolved.")

    return "resolved_component_info.txt"


def write_res_genome_info(all_resolved_paths, output):
    """
    Write resolved genome information to file
    """

    with open(f"{output}/resolved_genome_info.txt", "w") as myfile:
        myfile.write(f"Path\tCase\tCoverage\tLength\tGC content\tNode order\n")
        for genomic_path in all_resolved_paths:
            myfile.write(
                f"{genomic_path.id}\t{genomic_path.bubble_case}\t{genomic_path.coverage}\t{genomic_path.length}\t{genomic_path.gc}\t{genomic_path.node_order}\n"
            )

    return "resolved_genome_info.txt"


def write_path(final_genomic_paths, output):
    """
    Write genomic paths to a single FASTA file
    """

    with open(f"{output}/resolved_paths.fasta", "a+") as myfile:
        for genomic_path in final_genomic_paths:
            myfile.write(f">{genomic_path.id}\n")

            chunks = [
                genomic_path.path[i : i + FASTA_LINE_LEN]
                for i in range(0, genomic_path.length, FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_path_fasta(final_genomic_paths, output_genomes_path):
    """
    Write genomic paths to individual FASTA files
    """

    if not os.path.isdir(f"{output_genomes_path}"):
        subprocess.run("mkdir -p " + output_genomes_path, shell=True)

    for genomic_path in final_genomic_paths:
        with open(f"{output_genomes_path}/{genomic_path.id}.fasta", "w+") as myfile:
            myfile.write(f">{genomic_path.id}\n")

            chunks = [
                genomic_path.path[i : i + FASTA_LINE_LEN]
                for i in range(0, genomic_path.length, FASTA_LINE_LEN)
            ]

            for chunk in chunks:
                myfile.write(f"{chunk}\n")


def write_component_vog_info(resolved_components, comp_vogs, output):
    """
    Write VOGs found in resolved components
    """

    with open(f"{output}/component_vogs.txt", "w") as myfile:
        myfile.write(f"Component\tvog\n")
        for comp in resolved_components:
            myfile.write(f"{comp}\t{comp_vogs[comp]}\n")

    return "component_vogs.txt"


def init_files(output):
    """
    Initialise files and folders
    """

    open(f"{output}/resolved_edges.fasta", "a").close()
    open(f"{output}/resolved_paths.fasta", "a").close()
    open(f"{output}/resolved_genome_info.txt", "a").close()
    open(f"{output}/resolved_component_info.txt", "a").close()
    open(f"{output}/component_vogs.txt", "a").close()

    # if not os.path.isdir(f"{output}/resolved_paths"):
    #     subprocess.run(f"mkdir -p {output}/resolved_paths", shell=True)
