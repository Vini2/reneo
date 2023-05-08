#!/usr/bin/env python3

import logging
import sys
import time

import networkx as nx
from igraph import *
from reneo_utils import component_utils, edge_graph_utils, flow_utils, gene_utils
from reneo_utils.coverage_utils import get_junction_pe_coverage, get_unitig_coverage
from reneo_utils.genome_utils import GenomeComponent, GenomePath
from reneo_utils.output_utils import (
    write_component_info,
    write_component_vog_info,
    write_path,
    write_path_fasta,
    write_res_genome_info,
    write_unitigs,
)
from tqdm import tqdm

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2022, Reneo Project"
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"

MAX_VAL = sys.maxsize
LEN_THRESHOLD = 0.95

# Reneo main code
# ----------------------------------------------------------------------


def main():
    # Get arguments
    # ----------------------------------------------------------------------
    graph = snakemake.params.graph
    coverage = snakemake.params.coverage
    bampath = snakemake.params.bampath
    hmmout = snakemake.params.hmmout
    vogs = snakemake.params.vogs
    minlength = int(snakemake.params.minlength)
    mincov = int(snakemake.params.mincov)
    compcount = int(snakemake.params.compcount)
    maxpaths = int(snakemake.params.maxpaths)
    mgfrac = float(snakemake.params.mgfrac)
    evalue = float(snakemake.params.evalue)
    hmmscore = float(snakemake.params.hmmscore)
    output = snakemake.params.output
    log = snakemake.params.log

    # Setup logger
    # ----------------------------------------------------------------------

    logger = logging.getLogger(__version__)
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    # Setup output path for log file
    if log is None:
        fileHandler = logging.FileHandler(f"{output}/reneo.log")
    else:
        fileHandler = logging.FileHandler(f"{log}")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to Reneo: Unraveling Viral Genomes from Metagenomes."
    )

    logger.info(f"Input arguments: ")
    logger.info(f"Assembly graph file: {graph}")
    logger.info(f"Unitig coverage file: {coverage}")
    logger.info(f"BAM files path: {bampath}")
    logger.info(f"Unitig .hmmout file: {hmmout}")
    logger.info(f"Unitig vog annotations file: {vogs}")
    logger.info(f"Minimum length of unitigs to consider: {minlength}")
    logger.info(f"Minimum coverage of paths to output: {mincov}")
    logger.info(f"Minimum unitig count to consider a component: {compcount}")
    logger.info(f"Maximum number of paths to resolve for a component: {maxpaths}")
    logger.info(f"Length threshold to consider single copy marker genes: {mgfrac}")
    logger.info(f"Maximum e-value for vog annotations: {evalue}")
    logger.info(f"Minimum hmm score for vog annotations: {hmmscore}")
    logger.info(f"Output folder: {output}")

    start_time = time.time()

    # Get assembly graph
    # ----------------------------------------------------------------------
    (
        assembly_graph,
        oriented_links,
        unitig_names,
        unitig_names_rev,
        graph_unitigs,
        self_looped_nodes,
        edges_lengths,
    ) = edge_graph_utils.build_assembly_graph(graph)

    logger.info(
        f"Total number of vertices in the assembly graph: {len(assembly_graph.vs)}"
    )
    logger.info(
        f"Total number of links in the assembly graph: {len(assembly_graph.es)}"
    )

    # Get circular unitigs
    # ----------------------------------------------------------------------
    circular = edge_graph_utils.get_circular(self_looped_nodes, graph_unitigs)

    # Get unitigs with bacterial single copy marker genes
    # ----------------------------------------------------------------------
    smg_unitigs = gene_utils.get_smg_unitigs(hmmout, mgfrac)

    # Get unitigs with VOGs
    # ----------------------------------------------------------------------
    unitig_vogs = gene_utils.get_vog_unitigs(
        vogs, evalue, hmmscore
    )

    # Get components with viral bubbles
    # ----------------------------------------------------------------------
    pruned_vs, comp_vogs = component_utils.get_components(
        assembly_graph,
        unitig_names,
        smg_unitigs,
        unitig_vogs,
        circular,
        edges_lengths,
        minlength,
    )
    logger.info(f"Total number of components found: {len(pruned_vs)}")

    # Get unitig and junction pe coverages
    # ----------------------------------------------------------------------

    unitig_coverages = get_unitig_coverage(coverage)
    junction_pe_coverage = get_junction_pe_coverage(bampath, output)

    # Resolve genomes
    # ----------------------------------------------------------------------

    resolved_edges = set()

    all_resolved_paths = []

    all_components = []

    cycle_components = set()
    resolved_components = set()
    circular_unitigs = set()
    resolved_cyclic = set()

    case1_found = set()
    case1_resolved = set()
    case2_found = set()
    case2_resolved = set()
    case3_found = set()
    case3_resolved = set()

    virus_like_edges = set()
    all_virus_like_edges = set()

    for my_count in tqdm(pruned_vs, desc="Resolving components"):
        component_time_start = time.time()

        my_genomic_paths = []

        original_candidate_nodes = pruned_vs[my_count]

        candidate_nodes = pruned_vs[my_count]

        pruned_graph = assembly_graph.subgraph(candidate_nodes)

        has_cycles = False

        # if 423 not in candidate_nodes:
        #     continue

        logger.debug(f"my_count: {my_count}")

        logger.debug(f"number of unitigs: {len(candidate_nodes)}")
        logger.debug(f"{candidate_nodes}")

        all_virus_like_edges = all_virus_like_edges.union(set(candidate_nodes))

        in_degree = []
        out_degree = []

        # Case 2 components
        if len(candidate_nodes) == 2:
            all_self_looped = True

            for node in candidate_nodes:
                if unitig_names[node] not in self_looped_nodes:
                    all_self_looped = False

            if all_self_looped:
                case2_found.add(my_count)

                cycle_components.add(my_count)

                virus_like_edges = virus_like_edges.union(set(candidate_nodes))

                unitig1 = ""
                unitig2 = ""

                for edge in pruned_graph.es:
                    source_vertex_id = edge.source
                    target_vertex_id = edge.target

                    if source_vertex_id != target_vertex_id:
                        unitig1 = candidate_nodes[source_vertex_id]
                        unitig2 = candidate_nodes[target_vertex_id]

                if unitig1 != "" and unitig2 != "":
                    unitig1_name = unitig_names[unitig1]
                    unitig2_name = unitig_names[unitig2]

                    unitig1_len = len(str(graph_unitigs[unitig1_name]))
                    unitig2_len = len(str(graph_unitigs[unitig2_name]))

                    unitig_to_consider = -1
                    unitig_name = ""

                    repeat_unitig = -1
                    repeat_unitig_name = ""

                    if unitig1_len > unitig2_len and unitig1_len > minlength:
                        unitig_to_consider = unitig1
                        unitig_name = unitig1_name
                        repeat_unitig = unitig2
                        repeat_unitig_name = unitig2_name
                    elif unitig2_len > unitig1_len and unitig2_len > minlength:
                        unitig_to_consider = unitig2
                        unitig_name = unitig2_name
                        repeat_unitig = unitig1
                        repeat_unitig_name = unitig1_name

                    if unitig_to_consider != -1:
                        logger.debug(
                            f"Case 2 bubble: {unitig1_name} is {unitig1_len} bp long and {unitig2_name} is {unitig2_len} bp long."
                        )
                        cycle_number = 1
                        resolved_edges.add(unitig_to_consider)
                        resolved_edges.add(repeat_unitig)
                        path_string = (
                            str(graph_unitigs[repeat_unitig_name])
                            + str(graph_unitigs[unitig_name])
                            + str(graph_unitigs[repeat_unitig_name].reverse_complement())
                        )
                        logger.debug(
                            f"Terminal repeat detected is {repeat_unitig_name}"
                        )

                        genome_path = GenomePath(
                            f"virus_comp_{my_count}_cycle_{cycle_number}",
                            "case2",
                            [
                                f"{repeat_unitig_name}+",
                                f"{unitig_name}+",
                                f"{repeat_unitig_name}-",
                            ],
                            [repeat_unitig, unitig_to_consider, repeat_unitig],
                            path_string,
                            int(unitig_coverages[unitig_name]),
                            len(path_string),
                            (path_string.count("G") + path_string.count("C"))
                            / len(path_string)
                            * 100,
                        )
                        my_genomic_paths.append(genome_path)
                        resolved_components.add(my_count)
                        resolved_cyclic.add(my_count)
                        case2_resolved.add(my_count)

        # Case 3 components
        elif len(candidate_nodes) > 2 and len(candidate_nodes) <= compcount:
            # Create initial directed graph with coverage values
            # ----------------------------------------------------------------------
            G_edge = nx.DiGraph()

            my_counter = 0

            node_indices = {}
            node_indices_rev = {}

            cycle_edges = {}

            clean_node_count = 0

            for vertex in pruned_graph.vs["id"]:
                unitig_name = unitig_names[vertex]

                if unitig_name not in self_looped_nodes:
                    clean_node_count += 1

                for node in oriented_links[unitig_name]:
                    consider_edge = False

                    if not (
                        unitig_name in self_looped_nodes and node in self_looped_nodes
                    ):
                        consider_edge = True

                    if consider_edge:
                        cov_1 = MAX_VAL
                        cov_2 = MAX_VAL

                        if unitig_name in unitig_coverages:
                            cov_1 = unitig_coverages[unitig_name]
                        if node in unitig_coverages:
                            cov_2 = unitig_coverages[node]

                        min_cov = min([cov_1, cov_2])

                        for edge in oriented_links[unitig_name][node]:
                            cycle_edges[(unitig_name + edge[0], node + edge[1])] = int(
                                min_cov
                            )

            logger.debug(f"clean_node_count: {clean_node_count}")

            for cedge in cycle_edges:
                G_edge.add_edge(cedge[0], cedge[1], weight=cycle_edges[cedge])

            two_comp = sorted(nx.weakly_connected_components(G_edge), key=len)
            logger.debug(f"No. of weakly connected components: {len(two_comp)}")

            if len(two_comp) >= 2:
                G_edge.remove_nodes_from(list(two_comp[0]))

            try:
                cycles_found = nx.find_cycle(G_edge, orientation="original")
                if len(cycles_found) > 0:
                    has_cycles = True
            except nx.exception.NetworkXNoCycle:
                logger.debug(f"No cycles found in component {my_count}")

            if has_cycles:
                logger.debug(
                    f"Potentially cycles can be detected in component {my_count}."
                )

                # Remove dead-ends (nodes with no incoming or no outgoing edges)
                # ----------------------------------------------------------------------
                dead_ends_to_remove = edge_graph_utils.remove_dead_ends(G_edge)

                if len(dead_ends_to_remove) > 0:
                    for node in dead_ends_to_remove:
                        node_id = unitig_names_rev[node[:-1]]
                        if node_id in candidate_nodes:
                            candidate_nodes.remove(node_id)

                    G_edge.remove_nodes_from(dead_ends_to_remove)

                    logger.debug(f"Dead-ends found and removed: {dead_ends_to_remove}")

                # Identify source/sink vertex
                # ----------------------------------------------------------------------

                source_sink_candidates = flow_utils.get_source_sink(
                    G_edge, graph_unitigs, minlength, self_looped_nodes
                )

                source_sink = 0

                logger.debug(f"Original candidate_nodes: {candidate_nodes}")
                logger.debug(
                    f"Identified candidate source_sinks from BFS: {source_sink_candidates}"
                )

                if len(source_sink_candidates) > 0:
                    # Identify the longest source/sink vertex
                    max_length = -1
                    max_length_st_vertex = -1

                    for vertex in source_sink_candidates:
                        if len(graph_unitigs[vertex[:-1]]) > max_length:
                            max_length = len(graph_unitigs[vertex[:-1]])
                            max_length_st_vertex = vertex

                    source_sink = unitig_names_rev[max_length_st_vertex[:-1]]
                    logger.debug(f"Identified source_sink from BFS: {source_sink}, {max_length_st_vertex}")

                    candidate_nodes.remove(source_sink)
                    candidate_nodes.insert(0, source_sink)
                    logger.debug(f"Ordered candidate_nodes: {candidate_nodes}")

                else:
                    logger.debug(f"No source/sink node detected")
                    continue

                # Create refined directed graph for flow network
                # ----------------------------------------------------------------------
                G = nx.DiGraph()

                for u, v, cov in G_edge.edges(data=True):
                    if u not in node_indices:
                        node_indices[u] = my_counter
                        node_indices_rev[my_counter] = u
                        my_counter += 1
                    if v not in node_indices:
                        node_indices[v] = my_counter
                        node_indices_rev[my_counter] = v
                        my_counter += 1

                    logger.debug(f"Edge: {u}, {v}, {cov['weight']}")

                    G.add_edge(node_indices[u], node_indices[v], weight=cov["weight"])

                # Get connections and degree information
                # ----------------------------------------------------------------------
                in_degree = []
                out_degree = []

                for node in list(G.nodes):
                    if node_indices_rev[node][:-1] not in self_looped_nodes:
                        clean_indegree = len(
                            [
                                x
                                for x in G.predecessors(node)
                                if node_indices_rev[x][:-1] not in self_looped_nodes
                            ]
                        )
                        in_degree.append(clean_indegree)

                        clean_outdegree = len(
                            [
                                x
                                for x in G.successors(node)
                                if node_indices_rev[x][:-1] not in self_looped_nodes
                            ]
                        )
                        out_degree.append(clean_outdegree)

                degrees = in_degree + out_degree

                if len(degrees) == 0:
                    logger.debug(
                        f"Skipping component as no clean connections were found"
                    )
                    continue

                # Create flow network
                # ----------------------------------------------------------------------
                network_edges = list()

                edge_list_indices = {}

                subpaths = {}
                subpath_count = 0

                visited_edges = []

                logger.debug(f"G_edge.nodes: {list(G_edge.nodes)}")
                logger.debug(f"G_edge.edges: {G_edge.edges(data=True)}")

                for u, v, cov in G_edge.edges(data=True):

                    u_name = unitig_names_rev[u[:-1]]
                    v_name = unitig_names_rev[v[:-1]]

                    original_edge = (u[:-1], v[:-1])
                    rev_original_edge = (v[:-1], u[:-1])

                    u_index = candidate_nodes.index(u_name)
                    v_index = candidate_nodes.index(v_name)

                    edge_list_indices[u_index] = u
                    edge_list_indices[v_index] = v

                    juction_cov = junction_pe_coverage[(u[:-1], v[:-1])]

                    if v_index == 0:
                        if (u_index, len(candidate_nodes)) not in visited_edges and (
                            len(candidate_nodes),
                            u_index,
                        ) not in visited_edges:
                            if juction_cov < cov["weight"]:
                                network_edges.append(
                                    (
                                        u_index,
                                        len(candidate_nodes),
                                        juction_cov,
                                        cov["weight"],
                                    )
                                )
                            else:
                                network_edges.append(
                                    (u_index, len(candidate_nodes), 0, cov["weight"])
                                )

                            visited_edges.append((u_index, len(candidate_nodes)))

                            if juction_cov != 0:
                                subpaths[subpath_count] = [
                                    u_index,
                                    len(candidate_nodes),
                                ]
                                subpath_count += 1
                    else:
                        if (u_index, v_index) not in visited_edges and (
                            v_index,
                            u_index,
                        ) not in visited_edges:
                            cov_upper_bound = cov["weight"]

                            if juction_cov <= cov_upper_bound:
                                network_edges.append(
                                    (u_index, v_index, juction_cov, cov_upper_bound)
                                )
                            else:
                                network_edges.append(
                                    (u_index, v_index, 0, cov_upper_bound)
                                )

                            visited_edges.append((u_index, v_index))

                            if juction_cov != 0:
                                subpaths[subpath_count] = [u_index, v_index]
                                subpath_count += 1

                logger.debug(f"edge_list_indices: {edge_list_indices}")

                # Create flow network and run MFD-ILP
                # ----------------------------------------------------------------------
                G_mfd = {
                    "Nodes": len(list(G_edge.nodes)),
                    "list of edges": network_edges,
                    "subpaths": subpaths,
                }
                logger.debug(f"G_mfd: {G_mfd}")
                solution_paths = flow_utils.solve_mfd(G_mfd, maxpaths, output)
                logger.debug(f"Number of paths found: {len(solution_paths)}")

                cycle_components.add(my_count)
                case3_found.add(my_count)

                # Iterate through solution paths
                # ----------------------------------------------------------------------
                if len(solution_paths) != 0:
                    virus_like_edges = virus_like_edges.union(
                        set(original_candidate_nodes)
                    )

                    cycle_number = 1

                    for solution_path in solution_paths:
                        coverage_val = solution_paths[solution_path]["weight"]

                        logger.debug(f"Path {cycle_number} coverage: {coverage_val}")

                        # Filter path by coverage
                        if coverage_val >= mincov:
                            # Create graph for path
                            G_path = nx.DiGraph()

                            # Fill graph with data
                            G_path.add_edges_from(solution_paths[solution_path]["path"])
                            logger.debug(
                                f"solution path: {solution_paths[solution_path]['path']}"
                            )

                            if 0 in list(G_path.nodes):
                                # Get all simple paths from node 0 to last node
                                candidate_paths = list(
                                    nx.all_simple_paths(G_path, 0, len(candidate_nodes))
                                )

                                if len(candidate_paths) > 0:
                                    logger.debug(
                                        f"candidate_paths: {candidate_paths[0]}"
                                    )

                                    # Get mapped unitigs in order from the flow network
                                    path_order = []
                                    for path_edge in candidate_paths[0]:
                                        if path_edge != len(candidate_nodes):
                                            path_order.append(
                                                edge_list_indices[path_edge]
                                            )

                                    logger.debug(f"path_order: {path_order}")

                                    # Get the order of unitigs in path
                                    path_string = ""
                                    total_length = 0

                                    for node in path_order:
                                        unitig_name = node[:-1]
                                        if node.endswith("+"):
                                            path_string += str(
                                                graph_unitigs[unitig_name]
                                            )
                                        else:
                                            path_string += str(
                                                graph_unitigs[
                                                    unitig_name
                                                ].reverse_complement()
                                            )
                                        total_length += len(
                                            str(graph_unitigs[unitig_name])
                                        )

                                    # Create GenomePath object with path details
                                    genome_path = GenomePath(
                                        f"virus_comp_{my_count}_cycle_{cycle_number}",
                                        "case3",
                                        [x for x in path_order],
                                        [unitig_names_rev[x[:-1]] for x in path_order],
                                        path_string,
                                        int(coverage_val),
                                        total_length,
                                        (
                                            path_string.count("G")
                                            + path_string.count("C")
                                        )
                                        / len(path_string)
                                        * 100,
                                    )
                                    my_genomic_paths.append(genome_path)
                                    logger.debug(f"total_length: {total_length}")

                                    cycle_number += 1

                    logger.debug(f"Number of paths selected: {cycle_number-1}")

                    if cycle_number > 1:
                        resolved_components.add(my_count)
                        resolved_cyclic.add(my_count)
                        case3_resolved.add(my_count)

                else:
                    logger.debug(f"No paths detected")
                    continue

            else:
                logger.debug(f"No cycles detected. Maybe a linear component?")

        # Case 1 components - circular unitigs
        elif len(candidate_nodes) == 1:
            case1_found.add(my_count)

            unitig_name = unitig_names[candidate_nodes[0]]

            resolved_edges.add(candidate_nodes[0])

            path_string = str(graph_unitigs[unitig_name])

            cycle_number = 1

            # Create GenomePath object with path details
            genome_path = GenomePath(
                f"virus_comp_{my_count}_cycle_{cycle_number}",
                "case1",
                [unitig_names[candidate_nodes[0]]],
                [candidate_nodes[0]],
                path_string,
                int(unitig_coverages[unitig_name]),
                len(graph_unitigs[unitig_name]),
                (path_string.count("G") + path_string.count("C"))
                / len(path_string)
                * 100,
            )
            my_genomic_paths.append(genome_path)
            resolved_components.add(my_count)
            circular_unitigs.add(my_count)
            case1_resolved.add(my_count)

            virus_like_edges = virus_like_edges.union(set(candidate_nodes))

        # Record final paths for the component
        # ----------------------------------------------------------------------

        # Order resolved paths in descending order of length
        my_genomic_paths.sort(key=lambda x: (x.length, x.coverage), reverse=True)

        final_genomic_paths = []
        comp_resolved_edges = set()
        visited_nodes = set()
        comp_resolved_paths = set()

        frac_unitigs = 1
        n_paths = 0

        if len(my_genomic_paths) > 1:
            # Get the degree of the component
            graph_degree = assembly_graph.degree(original_candidate_nodes)

            path_lengths = []
            path_coverages = []

            largest_length = my_genomic_paths[0].length

            # Filter genomic paths
            for genomic_path in my_genomic_paths:

                passed = False

                if genomic_path.length > largest_length * LEN_THRESHOLD:
                    passed = True

                path_node_order_string = ",".join(genomic_path.node_order)

                if path_node_order_string in comp_resolved_paths:
                    passed = False

                if passed:
                    logger.debug(
                        f"{genomic_path.id}\t{genomic_path.length}\t{genomic_path.coverage}"
                    )
                    logger.debug(f"{genomic_path.node_order}")
                    path_lengths.append(genomic_path.length)
                    path_coverages.append(genomic_path.coverage)
                    final_genomic_paths.append(genomic_path)
                    visited_nodes = visited_nodes.union(set(genomic_path.node_order))
                    comp_resolved_paths.add(path_node_order_string)
                    n_paths += 1

                    for path_node in genomic_path.node_id_order:
                        comp_resolved_edges.add(path_node)

            frac_unitigs = len(visited_nodes) / len(original_candidate_nodes)

            logger.debug(f"frac_unitigs: {frac_unitigs}")

            # Filter components
            if (
                len(final_genomic_paths) > 1
                and len(in_degree) > 0
                and len(out_degree) > 0
            ):
                # Create GenomeComponent object with component details
                genome_comp = GenomeComponent(
                    f"virus_comp_{my_count}",
                    len(original_candidate_nodes),
                    len(final_genomic_paths),
                    max(graph_degree),
                    min(graph_degree),
                    max(in_degree),
                    max(out_degree),
                    sum(graph_degree) / len(graph_degree),
                    sum(in_degree) / len(in_degree),
                    sum(out_degree) / len(out_degree),
                    pruned_graph.density(loops=False),
                    max(path_lengths),
                    min(path_lengths),
                    max(path_lengths) / min(path_lengths),
                    path_lengths[path_coverages.index(max(path_coverages))],
                    path_lengths[path_coverages.index(min(path_coverages))],
                    path_lengths[path_coverages.index(max(path_coverages))]
                    / path_lengths[path_coverages.index(min(path_coverages))],
                    max(path_coverages),
                    min(path_coverages),
                    max(path_coverages) / min(path_coverages),
                    frac_unitigs,
                )
                all_components.append(genome_comp)
                resolved_edges = resolved_edges.union(comp_resolved_edges)

            if len(final_genomic_paths) > 0:
                resolved_cyclic.add(my_count)
                resolved_components.add(my_count)
                all_resolved_paths += final_genomic_paths
                component_elapsed_time = time.time() - component_time_start
                logger.debug(
                    f"Elapsed time to resolve component {my_count} with {len(original_candidate_nodes)} nodes: {component_elapsed_time} seconds"
                )

        else:
            # Circular unitigs
            for genomic_path in my_genomic_paths:
                final_genomic_paths.append(genomic_path)
                all_resolved_paths.append(genomic_path)
                logger.debug(f"{genomic_path.id}\t{genomic_path.length}")
                resolved_components.add(my_count)

        # Write genome path to file
        # ----------------------------------------------------------------------
        write_path(final_genomic_paths, output)
        write_path_fasta(final_genomic_paths, f"{output}/resolved_viruses")

    # Log final summary information
    # ----------------------------------------------------------------------
    logger.info(f"Total number of cyclic components found: {len(cycle_components)}")
    logger.info(f"Total number of cyclic components resolved: {len(resolved_cyclic)}")
    logger.info(f"Circular unitigs identified: {len(circular_unitigs)}")
    logger.info(
        f"Total number of cyclic components found including circular unitigs: {len(cycle_components) + len(circular_unitigs)}"
    )
    logger.info(
        f"Total number of components resolved: {len(circular_unitigs)+len(resolved_cyclic)}"
    )
    logger.info(f"Case 1 (resolved/found): {len(case1_resolved)}/{len(case1_found)}")
    logger.info(f"Case 2 (resolved/found): {len(case2_resolved)}/{len(case2_found)}")
    logger.info(f"Case 3 (resolved/found): {len(case3_resolved)}/{len(case3_found)}")
    logger.info(f"Total number of genomes resolved: {len(all_resolved_paths)}")

    if len(all_resolved_paths) == 0:
        logger.info(f"No genomes were resolved.")
    else:
        logger.info(f"Resolved genomes can be found in {output}/resolved_paths.fasta")

    # Write edges to file
    # ----------------------------------------------------------------------

    write_unitigs(
        virus_like_edges, unitig_names, graph_unitigs, "virus_like_edges", output
    )
    write_unitigs(
        all_virus_like_edges, unitig_names, graph_unitigs, "all_virus_like_edges", output
    )
    write_unitigs(resolved_edges, unitig_names, graph_unitigs, "resolved_edges", output)

    # Record path information
    # ----------------------------------------------------------------------

    filename = write_res_genome_info(all_resolved_paths, output)
    if len(all_resolved_paths) > 0:
        logger.info(f"Resolved genome information can be found in {output}/{filename}")

    # Record component information
    # ----------------------------------------------------------------------

    filename = write_component_info(all_components, output)
    if len(all_components) > 0:
        logger.info(f"Resolved component information can be found in {output}/{filename}")

    filename = write_component_vog_info(resolved_components, comp_vogs, output)
    if len(resolved_components) > 0:
        logger.info(f"VOGs found in resolved components can be found in {output}/{filename}")

    # Get elapsed time
    # ----------------------------------------------------------------------

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Print elapsed time for the process
    logger.info(f"Elapsed time: {elapsed_time} seconds")

    # Exit program
    # ----------------------------------------------------------------------

    logger.info("Thank you for using Reneo!")


if __name__ == "__main__":
    main()
