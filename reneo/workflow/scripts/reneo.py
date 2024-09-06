#!/usr/bin/env python3

import logging
import pickle
import queue
import sys
import threading
import time

import networkx as nx
from igraph import *
from reneo_utils import component_utils, edge_graph_utils, flow_utils, gene_utils
from reneo_utils.coverage_utils import get_unitig_coverage
from reneo_utils.genome_utils import GenomeComponent, GenomePath
from reneo_utils.output_utils import (
    init_files,
    write_component_info,
    write_component_vog_info,
    write_path,
    write_path_fasta,
    write_res_genome_info,
    write_unitigs,
)

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.5.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


def setup_logging(**kwargs):

    logging.basicConfig(
        filename=kwargs["log"],
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.captureWarnings(True)
    logger = logging.getLogger(f"reneo {__version__}")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    return logger


def results_dict():
    results = {
        "resolved_edges": set(),
        "all_resolved_paths": [],
        "all_components": [],
        "cycle_components": set(),
        "linear_components": set(),
        "resolved_components": set(),
        "resolved_linear": set(),
        "single_unitigs": set(),
        "resolved_cyclic": set(),
        "case1_found": set(),
        "case1_resolved": set(),
        "case2_found": set(),
        "case2_resolved": set(),
        "case3_found": set(),
        "case3_resolved": set(),
        "virus_like_edges": set(),
        "all_virus_like_edges": set(),
        "genome_path_sets": set(),
    }
    return results


def merge_results(orig_res, new_res):
    for key in orig_res.keys():
        if isinstance(orig_res[key], list):
            orig_res[key] += new_res[key]
        else:
            orig_res[key] = orig_res[key].union(new_res[key])
    return orig_res


def worker_resolve_components(component_queue, results_queue, **kwargs):
    results = results_dict()

    while True:
        my_count = component_queue.get()
        if my_count is None:
            break

        component_time_start = time.time()
        my_genomic_paths = []
        original_candidate_nodes = kwargs["pruned_vs"][my_count]
        candidate_nodes = kwargs["pruned_vs"][my_count]
        pruned_graph = kwargs["assembly_graph"].subgraph(candidate_nodes)
        has_cycles = False
        results["all_virus_like_edges"] = results["all_virus_like_edges"].union(
            set(candidate_nodes)
        )
        comp_all_edges = set(set(candidate_nodes))
        comp_resolved_edges = set()

        kwargs["logger"].debug(f"my_count: {my_count}")
        kwargs["logger"].debug(f"number of unitigs: {len(candidate_nodes)}")
        kwargs["logger"].debug(f"{candidate_nodes}")

        in_degree = []
        out_degree = []

        case_name = ""

        # Case 2 components
        if len(candidate_nodes) == 2:
            all_self_looped = True
            one_circular = False

            if (
                kwargs["unitig_names"][candidate_nodes[0]]
                in kwargs["self_looped_nodes"]
                and kwargs["unitig_names"][candidate_nodes[1]]
                in kwargs["self_looped_nodes"]
            ):
                all_self_looped = True
            else:
                if (
                    kwargs["unitig_names"][candidate_nodes[0]]
                    in kwargs["self_looped_nodes"]
                ):
                    one_circular = True
                    all_self_looped = False
                if (
                    kwargs["unitig_names"][candidate_nodes[1]]
                    in kwargs["self_looped_nodes"]
                ):
                    one_circular = True
                    all_self_looped = False

            unitig1 = ""
            unitig2 = ""

            for edge in pruned_graph.es:
                source_vertex_id = edge.source
                target_vertex_id = edge.target

                if source_vertex_id != target_vertex_id:
                    unitig1 = candidate_nodes[source_vertex_id]
                    unitig2 = candidate_nodes[target_vertex_id]

            unitig1_name = kwargs["unitig_names"][unitig1]
            unitig2_name = kwargs["unitig_names"][unitig2]

            unitig1_len = len(str(kwargs["graph_unitigs"][unitig1_name]))
            unitig2_len = len(str(kwargs["graph_unitigs"][unitig2_name]))

            if unitig1 != "" and unitig2 != "":

                # Case 2 - both are circular
                if all_self_looped:

                    case_name = "case2_circular"

                    results["case2_found"].add(my_count)

                    results["cycle_components"].add(my_count)

                    results["virus_like_edges"] = results["virus_like_edges"].union(
                        set(candidate_nodes)
                    )
                    comp_resolved_edges = comp_resolved_edges.union(
                        set(candidate_nodes)
                    )

                    unitig_to_consider = -1
                    unitig_name = ""

                    repeat_unitig = -1
                    repeat_unitig_name = ""

                    if unitig1_len > unitig2_len and unitig1_len > kwargs["minlength"]:
                        unitig_to_consider = unitig1
                        unitig_name = unitig1_name
                        repeat_unitig = unitig2
                        repeat_unitig_name = unitig2_name
                    elif (
                        unitig2_len > unitig1_len and unitig2_len > kwargs["minlength"]
                    ):
                        unitig_to_consider = unitig2
                        unitig_name = unitig2_name
                        repeat_unitig = unitig1
                        repeat_unitig_name = unitig1_name

                    if unitig_to_consider != -1:
                        kwargs["logger"].debug(
                            f"Case 2 component: {unitig1_name} is {unitig1_len} bp long and {unitig2_name} is {unitig2_len} bp long."
                        )
                        cycle_number = 1
                        results["resolved_edges"].add(unitig_to_consider)
                        results["resolved_edges"].add(repeat_unitig)
                        path_string = (
                            str(kwargs["graph_unitigs"][repeat_unitig_name])
                            + str(
                                kwargs["graph_unitigs"][unitig_name][
                                    kwargs["link_overlap"][
                                        (repeat_unitig, unitig_to_consider)
                                    ] :
                                ]
                            )
                            + str(
                                kwargs["graph_unitigs"][repeat_unitig_name][
                                    kwargs["link_overlap"][
                                        (unitig_to_consider, repeat_unitig)
                                    ] :
                                ]
                            )
                        )
                        kwargs["logger"].debug(
                            f"Terminal repeat detected is {repeat_unitig_name}"
                        )

                        genome_path = GenomePath(
                            id=f"virus_comp_{my_count}_cycle_{cycle_number}",
                            bubble_case=case_name,
                            node_order=[
                                f"{repeat_unitig_name}+",
                                f"{unitig_name}+",
                                f"{repeat_unitig_name}-",
                            ],
                            node_id_order=[
                                repeat_unitig,
                                unitig_to_consider,
                                repeat_unitig,
                            ],
                            path=path_string,
                            coverage=int(kwargs["unitig_coverages"][unitig_name]),
                            length=len(path_string),
                            gc=(path_string.count("G") + path_string.count("C"))
                            / len(path_string)
                            * 100,
                        )
                        my_genomic_paths.append(genome_path)
                        results["resolved_components"].add(my_count)
                        results["resolved_cyclic"].add(my_count)
                        results["case2_resolved"].add(my_count)

                # Case 2 - only one is circular
                elif one_circular:

                    case_name = "case2_linear"

                    results["case2_found"].add(my_count)

                    results["cycle_components"].add(my_count)

                    results["virus_like_edges"] = results["virus_like_edges"].union(
                        set(candidate_nodes)
                    )
                    comp_resolved_edges = comp_resolved_edges.union(
                        set(candidate_nodes)
                    )

                    unitig_to_consider = -1
                    unitig_name = ""

                    repeat_unitig = -1
                    repeat_unitig_name = ""

                    if (
                        unitig1_len > unitig2_len
                        and unitig1_len > kwargs["minlength"]
                        and unitig2_name in kwargs["self_looped_nodes"]
                    ):
                        unitig_to_consider = unitig1
                        unitig_name = unitig1_name
                        repeat_unitig = unitig2
                        repeat_unitig_name = unitig2_name
                    elif (
                        unitig2_len > unitig1_len
                        and unitig2_len > kwargs["minlength"]
                        and unitig1_name in kwargs["self_looped_nodes"]
                    ):
                        unitig_to_consider = unitig2
                        unitig_name = unitig2_name
                        repeat_unitig = unitig1
                        repeat_unitig_name = unitig1_name

                    if unitig_to_consider != -1:
                        kwargs["logger"].debug(
                            f"Case 2 component: {unitig1_name} is {unitig1_len} bp long and {unitig2_name} is {unitig2_len} bp long."
                        )
                        cycle_number = 1
                        results["resolved_edges"].add(unitig_to_consider)
                        results["resolved_edges"].add(repeat_unitig)

                        # Get repeat count
                        repeat_count = max(
                            int(
                                kwargs["unitig_coverages"][repeat_unitig_name]
                                / kwargs["unitig_coverages"][unitig_name]
                            ),
                            1,
                        )
                        kwargs["logger"].debug(f"Repeat count: {repeat_count}")

                        path_string = (
                            str(
                                kwargs["graph_unitigs"][unitig_name][
                                    kwargs["link_overlap"][
                                        (repeat_unitig, unitig_to_consider)
                                    ] :
                                ]
                            )
                            + str(
                                kwargs["graph_unitigs"][repeat_unitig_name][
                                    kwargs["link_overlap"][
                                        (unitig_to_consider, repeat_unitig)
                                    ] :
                                ]
                            )
                            * repeat_count
                        )
                        kwargs["logger"].debug(
                            f"Terminal repeat detected is {repeat_unitig_name}"
                        )

                        # Format path node order
                        path_with_repeats = [f"{unitig_name}+"] + [
                            f"{repeat_unitig_name}+" for x in range(repeat_count)
                        ]

                        repeat_order = f"{repeat_unitig_name}:fwd," * repeat_count
                        path_with_repeats_human = (
                            f"{unitig_name}:fwd,{repeat_order[:-1]}"
                        )
                        node_id_order_with_repeats = [unitig_to_consider] + [
                            repeat_unitig for x in range(repeat_count)
                        ]

                        genome_path = GenomePath(
                            id=f"virus_comp_{my_count}_cycle_{cycle_number}",
                            bubble_case=case_name,
                            node_order=path_with_repeats,
                            node_id_order=node_id_order_with_repeats,
                            path=path_string,
                            coverage=int(kwargs["unitig_coverages"][unitig_name]),
                            length=len(path_string),
                            gc=(path_string.count("G") + path_string.count("C"))
                            / len(path_string)
                            * 100,
                        )
                        my_genomic_paths.append(genome_path)
                        results["resolved_components"].add(my_count)
                        results["resolved_linear"].add(my_count)
                        results["case2_resolved"].add(my_count)

        # Case 3 components
        elif len(candidate_nodes) > 2 and len(candidate_nodes) <= kwargs["compcount"]:

            case_name = "case3_circular"

            # Create initial directed graph with coverage values
            # ----------------------------------------------------------------------
            G_edge = nx.DiGraph()

            my_counter = 0

            node_indices = {}
            node_indices_rev = {}

            cycle_edges = {}

            clean_node_count = 0

            max_comp_cov = -1

            for vertex in pruned_graph.vs["id"]:
                unitig_name = kwargs["unitig_names"][vertex]

                # Find the maximum coverage within the component
                if (
                    unitig_name in kwargs["unitig_coverages"]
                    and kwargs["unitig_coverages"][unitig_name] > max_comp_cov
                ):
                    max_comp_cov = kwargs["unitig_coverages"][unitig_name]

                if unitig_name not in kwargs["self_looped_nodes"]:
                    clean_node_count += 1

                for node in kwargs["oriented_links"][unitig_name]:
                    consider_edge = False

                    if not (
                        unitig_name in kwargs["self_looped_nodes"]
                        and node in kwargs["self_looped_nodes"]
                    ):
                        consider_edge = True

                    if consider_edge:
                        cov_1 = kwargs["MAX_VAL"]
                        cov_2 = kwargs["MAX_VAL"]

                        if unitig_name in kwargs["unitig_coverages"]:
                            cov_1 = kwargs["unitig_coverages"][unitig_name]
                        if node in kwargs["unitig_coverages"]:
                            cov_2 = kwargs["unitig_coverages"][node]

                        min_cov = min([cov_1, cov_2]) if min([cov_1, cov_2]) != 0 else max([cov_1, cov_2])

                        for edge in kwargs["oriented_links"][unitig_name][node]:
                            cycle_edges[(unitig_name + edge[0], node + edge[1])] = int(
                                min_cov
                            )

            kwargs["logger"].debug(f"clean_node_count: {clean_node_count}")

            for cedge in cycle_edges:
                G_edge.add_edge(cedge[0], cedge[1], weight=cycle_edges[cedge])

            two_comp = sorted(nx.weakly_connected_components(G_edge), key=len)
            kwargs["logger"].debug(
                f"No. of weakly connected components: {len(two_comp)}"
            )

            if len(two_comp) >= 2:
                G_edge.remove_nodes_from(list(two_comp[0]))

            try:
                cycles_found = nx.find_cycle(G_edge, orientation="original")
                if len(cycles_found) > 0:
                    has_cycles = True
            except nx.exception.NetworkXNoCycle:
                kwargs["logger"].debug(f"No cycles found in component {my_count}")

            if has_cycles:
                kwargs["logger"].debug(
                    f"Potentially cycles can be detected in component {my_count}."
                )

                # Remove dead-ends (nodes with no incoming or no outgoing edges)
                # ----------------------------------------------------------------------
                dead_ends_to_remove = edge_graph_utils.remove_dead_ends(G_edge)

                if len(dead_ends_to_remove) > 0:
                    for node in dead_ends_to_remove:
                        node_id = kwargs["unitig_names_rev"][node[:-1]]
                        if node_id in candidate_nodes:
                            candidate_nodes.remove(node_id)

                    G_edge.remove_nodes_from(dead_ends_to_remove)

                    kwargs["logger"].debug(
                        f"Dead-ends found and removed: {dead_ends_to_remove}"
                    )

                # Identify source/sink vertex
                # ----------------------------------------------------------------------

                source_sink_candidates = flow_utils.get_source_sink_circular(
                    G_edge,
                    kwargs["graph_unitigs"],
                    kwargs["minlength"],
                    kwargs["self_looped_nodes"],
                )

                kwargs["logger"].debug(f"Original candidate_nodes: {candidate_nodes}")
                kwargs["logger"].debug(
                    f"Identified candidate source_sinks from BFS: {source_sink_candidates}"
                )

                if len(source_sink_candidates) > 0:
                    # Identify the longest source/sink vertex
                    max_length = -1
                    max_length_st_vertex = -1

                    for vertex in source_sink_candidates:
                        if len(kwargs["graph_unitigs"][vertex[:-1]]) > max_length:
                            max_length = len(kwargs["graph_unitigs"][vertex[:-1]])
                            max_length_st_vertex = vertex

                    source_sink = kwargs["unitig_names_rev"][max_length_st_vertex[:-1]]
                    kwargs["logger"].debug(
                        f"Identified source_sink from BFS: {source_sink}, {max_length_st_vertex}"
                    )

                    candidate_nodes.remove(source_sink)
                    candidate_nodes.insert(0, source_sink)
                    kwargs["logger"].debug(
                        f"Ordered candidate_nodes: {candidate_nodes}"
                    )

                else:
                    kwargs["logger"].debug(f"No source/sink node detected")
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

                    kwargs["logger"].debug(f"Edge: {u}, {v}, {cov['weight']}")

                    G.add_edge(node_indices[u], node_indices[v], weight=cov["weight"])

                # Get connections and degree information
                # ----------------------------------------------------------------------
                in_degree = []
                out_degree = []

                for node in list(G.nodes):
                    if node_indices_rev[node][:-1] not in kwargs["self_looped_nodes"]:
                        clean_indegree = len(
                            [
                                x
                                for x in G.predecessors(node)
                                if node_indices_rev[x][:-1]
                                not in kwargs["self_looped_nodes"]
                            ]
                        )
                        in_degree.append(clean_indegree)

                        clean_outdegree = len(
                            [
                                x
                                for x in G.successors(node)
                                if node_indices_rev[x][:-1]
                                not in kwargs["self_looped_nodes"]
                            ]
                        )
                        out_degree.append(clean_outdegree)

                degrees = in_degree + out_degree

                if len(degrees) == 0:
                    kwargs["logger"].debug(
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

                junctions_visited = []

                kwargs["logger"].debug(f"G_edge.nodes: {list(G_edge.nodes)}")
                kwargs["logger"].debug(f"G_edge.edges: {G_edge.edges(data=True)}")

                for u, v, cov in G_edge.edges(data=True):
                    u_name = kwargs["unitig_names_rev"][u[:-1]]
                    v_name = kwargs["unitig_names_rev"][v[:-1]]

                    u_index = candidate_nodes.index(u_name)
                    v_index = candidate_nodes.index(v_name)

                    edge_list_indices[u_index] = u
                    edge_list_indices[v_index] = v

                    juction_cov = kwargs["junction_pe_coverage"][(u[:-1], v[:-1])]

                    if v_index == 0:
                        final_vertex = len(candidate_nodes)
                    else:
                        final_vertex = v_index

                    if (u_index, final_vertex) not in visited_edges and (
                        final_vertex,
                        u_index,
                    ) not in visited_edges:
                        # Get coverage interval
                        cov_lower_bound = cov["weight"]
                        cov_upper_bound = int(max_comp_cov * kwargs["alpha"])

                        kwargs["logger"].debug(
                            f"({v}, {u}), {juction_cov}, {cov_lower_bound}, {cov_upper_bound}"
                        )

                        if juction_cov == 0:
                            network_edges.append(
                                (u_index, final_vertex, cov_lower_bound, cov_upper_bound)
                            )
                        else:
                            network_edges.append(
                                (
                                    u_index,
                                    final_vertex,
                                    cov_lower_bound,
                                    cov_upper_bound,
                                )
                            )

                        visited_edges.append((u_index, final_vertex))

                        # Add subpaths based on junction coverage
                        if juction_cov >= kwargs["mincov"]:
                            kwargs["logger"].debug(
                                f"Adding subpath {[u_name, v_name]}"
                            )
                            subpaths[subpath_count] = [u_index, final_vertex]
                            subpath_count += 1

                        # Extend subpaths using coverages of successors and predecessors
                        else:
                            # Extend subpaths of l=3 based on paired-end reads
                            # aligned to successors and predecessors

                            for u_pred in G_edge.predecessors(u):
                                if (
                                    kwargs["junction_pe_coverage"][
                                        (u_pred[:-1], v[:-1])
                                    ] > 0
                                    and [u_pred, u, v] not in junctions_visited
                                ):
                                    u_pred_name = kwargs["unitig_names_rev"][
                                        u_pred[:-1]
                                    ]
                                    u_pred_index = candidate_nodes.index(u_pred_name)
                                    if (
                                        final_vertex != 0
                                        and u_index != 0
                                        and u_pred_index != final_vertex
                                    ):
                                        subpaths[subpath_count] = [
                                            u_pred_index,
                                            u_index,
                                            final_vertex,
                                        ]
                                        kwargs["logger"].debug(
                                            f"Extending subpath {[u_pred, u, v]}"
                                        )
                                        subpath_count += 1
                                        junctions_visited.append([u_pred, u, v])

                            for v_succ in G_edge.successors(v):
                                if (
                                    kwargs["junction_pe_coverage"][
                                        (u[:-1], v_succ[:-1])
                                    ] > 0
                                    and [u, v, v_succ] not in junctions_visited
                                ):
                                    v_succ_name = kwargs["unitig_names_rev"][
                                        v_succ[:-1]
                                    ]
                                    v_succ_index = candidate_nodes.index(v_succ_name)
                                    if (
                                        v_succ_index != 0
                                        and u_index != 0
                                        and final_vertex != 0
                                        and final_vertex != len(candidate_nodes)
                                        and v_succ_index != u_index
                                    ):
                                        subpaths[subpath_count] = [
                                            u_index,
                                            final_vertex,
                                            v_succ_index,
                                        ]
                                        kwargs["logger"].debug(
                                            f"Extending subpath {[u, v, v_succ]}"
                                        )
                                        subpath_count += 1
                                        junctions_visited.append([u, v, v_succ])

                        # Extend subpath using coverages of predecessors
                        kwargs["logger"].debug(
                            f"Predecessors of {u}: {[x for x in G_edge.predecessors(u)]}"
                        )
                        for u_pred in G_edge.predecessors(u):
                            u_pred_name = kwargs["unitig_names_rev"][u_pred[:-1]]
                            u_pred_index = candidate_nodes.index(u_pred_name)
                            u_pred_cov = kwargs["unitig_coverages"][u_pred[:-1]]
                            u_cov = kwargs["unitig_coverages"][u[:-1]]
                            v_cov = kwargs["unitig_coverages"][v[:-1]]

                            if (
                                [u_pred, u, v] not in junctions_visited
                                and final_vertex != 0
                                and u_index != 0
                                and u_pred_index != final_vertex
                            ):
                                kwargs["logger"].debug(
                                    f"Coverage across predecessor - {u_pred}:{u_pred_cov}, {u}:{u_cov}, {v}:{v_cov}"
                                )

                                left_weight = min(u_pred_cov, u_cov) if u_cov > 0 else max(u_pred_cov, u_cov)
                                right_weight = min(u_cov, v_cov) if u_cov > 0 else max(u_cov, v_cov)

                                if (
                                    abs(left_weight - right_weight)
                                    < kwargs["covtol"]
                                ):
                                    subpaths[subpath_count] = [
                                        u_pred_index,
                                        u_index,
                                        final_vertex,
                                    ]
                                    kwargs["logger"].debug(
                                        f"Extending subpath based on predecessor coverage {[u_pred, u, v]}, {[u_pred_cov, u_cov, v_cov]}"
                                    )
                                    subpath_count += 1

                                    junctions_visited.append([u_pred, u, v])

                        # Extend subpath using coverages of successors
                        kwargs["logger"].debug(
                            f"Successors of {v}: {[x for x in G_edge.successors(v)]}"
                        )
                        for v_succ in G_edge.successors(v):
                            v_succ_name = kwargs["unitig_names_rev"][v_succ[:-1]]
                            v_succ_index = candidate_nodes.index(v_succ_name)
                            v_succ_cov = kwargs["unitig_coverages"][v_succ[:-1]]
                            v_cov = kwargs["unitig_coverages"][v[:-1]]
                            u_cov = kwargs["unitig_coverages"][u[:-1]]

                            if (
                                [u, v, v_succ] not in junctions_visited
                                and v_succ_index != 0
                                and u_index != 0
                                and final_vertex != 0
                                and final_vertex != len(candidate_nodes)
                                and v_succ_index != u_index
                            ):
                                kwargs["logger"].debug(
                                    f"Coverage across successor - {u}:{u_cov}, {v}:{v_cov}, {v_succ}:{v_succ_cov}"
                                )

                                left_weight = min(u_cov, v_cov) if v_cov > 0 else max(u_cov, v_cov)
                                right_weight = min(v_succ_cov, v_cov) if v_cov > 0 else max(v_succ_cov, v_cov)

                                if (
                                    abs(left_weight - right_weight)
                                    < kwargs["covtol"]
                                ):
                                    subpaths[subpath_count] = [
                                        u_index,
                                        final_vertex,
                                        v_succ_index,
                                    ]
                                    kwargs["logger"].debug(
                                        f"Extending subpath based on successor coverage {[u, v, v_succ]}, {[u_cov, v_cov, v_succ_cov]}"
                                    )
                                    subpath_count += 1

                                    junctions_visited.append([u, v, v_succ])

                kwargs["logger"].debug(f"edge_list_indices: {edge_list_indices}")
                kwargs["logger"].debug(f"subpaths: {subpaths}")

                # Create flow network and run MFD-ILP
                # ----------------------------------------------------------------------
                G_mfd = {
                    "Nodes": len(list(G_edge.nodes)),
                    "list of edges": network_edges,
                    "subpaths": subpaths,
                }
                kwargs["logger"].debug(f"G_mfd: {G_mfd}")
                solution_paths = flow_utils.solve_mfd(G_mfd, kwargs["maxpaths"], 1)
                kwargs["logger"].debug(f"Number of paths found: {len(solution_paths)}")

                results["cycle_components"].add(my_count)
                results["case3_found"].add(my_count)

                # Iterate through solution paths
                # ----------------------------------------------------------------------
                if len(solution_paths) != 0:
                    results["virus_like_edges"] = results["virus_like_edges"].union(
                        set(original_candidate_nodes)
                    )

                    cycle_number = 1

                    for solution_path in solution_paths:
                        coverage_val = solution_paths[solution_path]["weight"]

                        # Filter path by coverage
                        if coverage_val >= kwargs["mincov"]:
                            kwargs["logger"].debug(
                                f"Path {cycle_number} coverage: {coverage_val}"
                            )

                            # Create graph for path
                            G_path = nx.DiGraph()

                            # Fill graph with data
                            G_path.add_edges_from(solution_paths[solution_path]["path"])
                            kwargs["logger"].debug(
                                f"solution path: {solution_paths[solution_path]['path']}"
                            )

                            if 0 in list(G_path.nodes):
                                # Get all simple paths from node 0 to last node
                                try:
                                    candidate_paths = list(
                                        nx.all_simple_paths(
                                            G_path, 0, len(candidate_nodes)
                                        )
                                    )

                                    if len(candidate_paths) > 0:
                                        kwargs["logger"].debug(
                                            f"candidate_paths: {candidate_paths[0]}"
                                        )

                                        # Get mapped unitigs in order from the flow network
                                        path_order = []
                                        for path_edge in candidate_paths[0]:
                                            if path_edge != len(candidate_nodes):
                                                path_order.append(
                                                    edge_list_indices[path_edge]
                                                )

                                        kwargs["logger"].debug(
                                            f"path_order: {path_order}"
                                        )

                                        # Get the order of unitigs in path
                                        path_string = ""
                                        total_length = 0

                                        previous_edge = 0

                                        for nodeid in range(len(path_order)):
                                            node = path_order[nodeid]
                                            unitig_name = node[:-1]

                                            if node.endswith("+"):
                                                unitig_seq = str(
                                                    kwargs["graph_unitigs"][unitig_name]
                                                )
                                            else:
                                                unitig_seq = str(
                                                    kwargs["graph_unitigs"][
                                                        unitig_name
                                                    ].reverse_complement()
                                                )

                                            # If first node in path
                                            if previous_edge == 0:
                                                path_string += unitig_seq
                                                total_length += len(unitig_seq)

                                            else:
                                                trimmed_seq = unitig_seq[
                                                    kwargs["link_overlap"][
                                                        (previous_edge, node)
                                                    ] :
                                                ]
                                                path_string += trimmed_seq
                                                total_length += len(trimmed_seq)

                                            previous_edge = node

                                        # Create GenomePath object with path details
                                        genome_path = GenomePath(
                                            id=f"virus_comp_{my_count}_cycle_{cycle_number}",
                                            bubble_case=case_name,
                                            node_order=[x for x in path_order],
                                            node_id_order=[
                                                kwargs["unitig_names_rev"][x[:-1]]
                                                for x in path_order
                                            ],
                                            path=path_string,
                                            coverage=int(coverage_val),
                                            length=total_length,
                                            gc=(
                                                path_string.count("G")
                                                + path_string.count("C")
                                            )
                                            / len(path_string)
                                            * 100,
                                        )
                                        my_genomic_paths.append(genome_path)
                                        kwargs["logger"].debug(
                                            f"total_length: {total_length}"
                                        )

                                        cycle_number += 1

                                except nx.exception.NodeNotFound:
                                    kwargs["logger"].debug(
                                        f"Could not resolve a continuous path."
                                    )

                    kwargs["logger"].debug(
                        f"Number of paths selected: {cycle_number - 1}"
                    )

                    if cycle_number > 1:
                        results["resolved_components"].add(my_count)
                        results["resolved_cyclic"].add(my_count)
                        results["case3_resolved"].add(my_count)

                else:
                    kwargs["logger"].debug(f"No paths detected")
                    continue

            else:
                kwargs["logger"].debug(
                    f"No cycles detected. Found a complex linear component."
                )

                case_name = "case3_linear"

                results["linear_components"].add(my_count)

                # Identify source/sink vertex
                # ----------------------------------------------------------------------

                source_candidates, sink_candidates = flow_utils.get_source_sink_linear(
                    G_edge, kwargs["self_looped_nodes"]
                )

                kwargs["logger"].debug(f"Original candidate_nodes: {candidate_nodes}")
                kwargs["logger"].debug(
                    f"Identified candidate sources: {source_candidates}"
                )
                kwargs["logger"].debug(f"Identified candidate sinks: {sink_candidates}")

                if len(source_candidates) > 0 and len(sink_candidates) > 0:

                    source_node_indices = [
                        kwargs["unitig_names_rev"][x[:-1]] for x in source_candidates
                    ]
                    sink_node_indices = [
                        kwargs["unitig_names_rev"][x[:-1]] for x in sink_candidates
                    ]

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

                        kwargs["logger"].debug(f"Edge: {u}, {v}, {cov['weight']}")

                        G.add_edge(
                            node_indices[u], node_indices[v], weight=cov["weight"]
                        )

                    # Get connections and degree information
                    # ----------------------------------------------------------------------
                    in_degree = []
                    out_degree = []

                    for node in list(G.nodes):
                        if (
                            node_indices_rev[node][:-1]
                            not in kwargs["self_looped_nodes"]
                        ):
                            clean_indegree = len(
                                [
                                    x
                                    for x in G.predecessors(node)
                                    if node_indices_rev[x][:-1]
                                    not in kwargs["self_looped_nodes"]
                                ]
                            )
                            in_degree.append(clean_indegree)

                            clean_outdegree = len(
                                [
                                    x
                                    for x in G.successors(node)
                                    if node_indices_rev[x][:-1]
                                    not in kwargs["self_looped_nodes"]
                                ]
                            )
                            out_degree.append(clean_outdegree)

                    degrees = in_degree + out_degree

                    if len(degrees) == 0:
                        kwargs["logger"].debug(
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

                    junctions_visited = []

                    kwargs["logger"].debug(f"G_edge.nodes: {list(G_edge.nodes)}")
                    kwargs["logger"].debug(f"G_edge.edges: {G_edge.edges(data=True)}")

                    for u, v, cov in G_edge.edges(data=True):
                        u_name = kwargs["unitig_names_rev"][u[:-1]]
                        v_name = kwargs["unitig_names_rev"][v[:-1]]

                        u_index = candidate_nodes.index(u_name) + 1
                        v_index = candidate_nodes.index(v_name) + 1

                        edge_list_indices[u_index] = u
                        edge_list_indices[v_index] = v

                        juction_cov = kwargs["junction_pe_coverage"][(u[:-1], v[:-1])]

                        if (u_index, v_index) not in visited_edges and (
                            v_index,
                            u_index,
                        ) not in visited_edges:
                            # Get coverage interval
                            cov_lower_bound = cov["weight"]
                            cov_upper_bound = int(max_comp_cov * kwargs["alpha"])

                            kwargs["logger"].debug(
                                f"({v}, {u}), ({u_index}, {v_index}), {juction_cov}, {cov_lower_bound}, {cov_upper_bound}"
                            )

                            if juction_cov == 0:
                                network_edges.append(
                                    (u_index, v_index, 0, cov_upper_bound)
                                )
                            else:
                                network_edges.append(
                                    (
                                        u_index,
                                        v_index,
                                        cov_lower_bound,
                                        cov_upper_bound,
                                    )
                                )

                            visited_edges.append((u_index, v_index))

                            # Add subpaths based on junction coverage
                            if juction_cov >= kwargs["mincov"]:
                                kwargs["logger"].debug(
                                    f"Adding subpath {[u_name, v_name]}"
                                )
                                subpaths[subpath_count] = [u_index, v_index]
                                subpath_count += 1

                            # Extend subpaths using coverages of successors and predecessors
                            else:
                                # Extend subpaths of l=3 based on paired-end reads
                                # aligned to successors and predecessors

                                for u_pred in G_edge.predecessors(u):
                                    if (
                                        kwargs["junction_pe_coverage"][
                                            (u_pred[:-1], v[:-1])
                                        ] > 0
                                        and [u_pred, u, v] not in junctions_visited
                                    ):
                                        u_pred_name = kwargs["unitig_names_rev"][
                                            u_pred[:-1]
                                        ]
                                        u_pred_index = (
                                            candidate_nodes.index(u_pred_name) + 1
                                        )
                                        if (
                                            (v_index - 1) not in source_node_indices
                                            and (u_index - 1) not in source_node_indices
                                            and u_pred_index != v_index
                                        ):
                                            subpaths[subpath_count] = [
                                                u_pred_index,
                                                u_index,
                                                v_index,
                                            ]
                                            kwargs["logger"].debug(
                                                f"Extending subpath {[u_pred, u, v]}"
                                            )
                                            subpath_count += 1
                                            junctions_visited.append([u_pred, u, v])

                                for v_succ in G_edge.successors(v):
                                    if (
                                        kwargs["junction_pe_coverage"][
                                            (u[:-1], v_succ[:-1])
                                        ] > 0
                                        and [u, v, v_succ] not in junctions_visited
                                    ):
                                        v_succ_name = kwargs["unitig_names_rev"][
                                            v_succ[:-1]
                                        ]
                                        v_succ_index = (
                                            candidate_nodes.index(v_succ_name) + 1
                                        )
                                        if (
                                            (v_succ_index - 1)
                                            not in source_node_indices
                                            and (u_index - 1) not in source_node_indices
                                            and (v_index - 1) not in source_node_indices
                                            and (v_index - 1) not in sink_node_indices
                                            and v_succ_index != u_index
                                        ):
                                            subpaths[subpath_count] = [
                                                u_index,
                                                v_index,
                                                v_succ_index,
                                            ]
                                            kwargs["logger"].debug(
                                                f"Extending subpath {[u, v, v_succ]}"
                                            )
                                            subpath_count += 1
                                            junctions_visited.append([u, v, v_succ])

                            # Extend subpath using coverages of predecessors
                            kwargs["logger"].debug(
                                f"Predecessors of {u}: {[x for x in G_edge.predecessors(u)]}"
                            )
                            for u_pred in G_edge.predecessors(u):
                                u_pred_name = kwargs["unitig_names_rev"][
                                    u_pred[:-1]
                                ]
                                u_pred_index = (
                                    candidate_nodes.index(u_pred_name) + 1
                                )
                                u_pred_cov = kwargs["unitig_coverages"][u_pred[:-1]]
                                u_cov = kwargs["unitig_coverages"][u[:-1]]
                                v_cov = kwargs["unitig_coverages"][v[:-1]]

                                if (
                                    [u_pred, u, v] not in junctions_visited
                                    and (v_index - 1) not in source_node_indices
                                    and (u_index - 1) not in source_node_indices
                                    and u_pred_index != v_index
                                ):
                                    kwargs["logger"].debug(
                                        f"Coverage across predecessor - {u_pred}:{u_pred_cov}, {u}:{u_cov}, {v}:{v_cov}"
                                    )

                                    left_weight = min(u_pred_cov, u_cov) if u_cov > 0 else max(u_pred_cov, u_cov)
                                    right_weight = min(u_cov, v_cov) if u_cov > 0 else max(u_cov, v_cov)

                                    if (
                                        abs(left_weight - right_weight)
                                        < kwargs["covtol"]
                                    ):
                                        subpaths[subpath_count] = [
                                            u_pred_index,
                                            u_index,
                                            v_index,
                                        ]
                                        kwargs["logger"].debug(
                                            f"Extending subpath based on predecessor coverage {[u_pred, u, v]}, {[u_pred_cov, u_cov, v_cov]}"
                                        )
                                        subpath_count += 1

                                        junctions_visited.append([u_pred, u, v])

                            # Extend subpath using coverages of successors
                            kwargs["logger"].debug(
                                f"Successors of {v}: {[x for x in G_edge.successors(v)]}"
                            )
                            for v_succ in G_edge.successors(v):
                                v_succ_name = kwargs["unitig_names_rev"][
                                    v_succ[:-1]
                                ]
                                v_succ_index = (
                                    candidate_nodes.index(v_succ_name) + 1
                                )
                                v_succ_cov = kwargs["unitig_coverages"][v_succ[:-1]]
                                v_cov = kwargs["unitig_coverages"][v[:-1]]
                                u_cov = kwargs["unitig_coverages"][v[:-1]]

                                if (
                                    [u, v, v_succ] not in junctions_visited
                                    and (v_succ_index - 1) not in source_node_indices
                                    and (u_index - 1) not in source_node_indices
                                    and (v_index - 1) not in source_node_indices
                                    and (v_index - 1) not in sink_node_indices
                                    and v_succ_index != u_index
                                ):
                                    kwargs["logger"].debug(
                                        f"Coverage across successor - {u}:{u_cov}, {v}:{v_cov}, {v_succ}:{v_succ_cov}"
                                    )

                                    left_weight = min(u_cov, v_cov) if v_cov > 0 else max(u_cov, v_cov)
                                    right_weight = min(v_succ_cov, v_cov) if v_cov > 0 else max(v_succ_cov, v_cov)

                                    if (
                                        abs(left_weight - right_weight)
                                        < kwargs["covtol"]
                                    ):
                                        subpaths[subpath_count] = [
                                            u_index,
                                            v_index,
                                            v_succ_index,
                                        ]
                                        kwargs["logger"].debug(
                                            f"Extending subpath based on successor coverage {[u, v, v_succ]}, {[u_cov, v_cov, v_succ_cov]}"
                                        )
                                        subpath_count += 1

                                        junctions_visited.append([u, v, v_succ])

                    # Add common start to source links
                    for source_v in source_candidates:
                        source_node_index = (
                            candidate_nodes.index(
                                kwargs["unitig_names_rev"][source_v[:-1]]
                            )
                            + 1
                        )
                        source_node_cov = kwargs["unitig_coverages"][source_v[:-1]]
                        cov_upper_bound = int(max_comp_cov * kwargs["alpha"])

                        network_edges.append(
                            (
                                0,
                                source_node_index,
                                source_node_cov,
                                cov_upper_bound,
                            )
                        )

                        subpaths[subpath_count] = [0, source_node_index]
                        subpath_count += 1

                    # Add common sink to end links
                    for sink_v in sink_candidates:
                        sink_node_index = (
                            candidate_nodes.index(
                                kwargs["unitig_names_rev"][sink_v[:-1]]
                            )
                            + 1
                        )
                        sink_node_cov = kwargs["unitig_coverages"][sink_v[:-1]]
                        cov_upper_bound = int(max_comp_cov * kwargs["alpha"])

                        network_edges.append(
                            (
                                sink_node_index,
                                len(candidate_nodes) + 1,
                                sink_node_cov,
                                cov_upper_bound,
                            )
                        )

                        subpaths[subpath_count] = [
                            sink_node_index,
                            len(candidate_nodes) + 1,
                        ]
                        subpath_count += 1

                    kwargs["logger"].debug(f"edge_list_indices: {edge_list_indices}")
                    kwargs["logger"].debug(f"subpaths: {subpaths}")

                    # Create flow network and run MFD-ILP
                    # ----------------------------------------------------------------------
                    G_mfd = {
                        "Nodes": len(list(G_edge.nodes)),
                        "list of edges": network_edges,
                        "subpaths": subpaths,
                    }
                    kwargs["logger"].debug(f"G_mfd: {G_mfd}")
                    solution_paths = flow_utils.solve_mfd(G_mfd, kwargs["maxpaths"], 1)
                    kwargs["logger"].debug(
                        f"Number of paths found: {len(solution_paths)}"
                    )

                    results["case3_found"].add(my_count)

                    # Iterate through solution paths
                    # ----------------------------------------------------------------------
                    if len(solution_paths) != 0:
                        results["virus_like_edges"] = results["virus_like_edges"].union(
                            set(original_candidate_nodes)
                        )

                        cycle_number = 1

                        for solution_path in solution_paths:
                            coverage_val = solution_paths[solution_path]["weight"]

                            # Filter path by coverage
                            if coverage_val >= kwargs["mincov"]:
                                kwargs["logger"].debug(
                                    f"Path {cycle_number} coverage: {coverage_val}"
                                )

                                # Create graph for path
                                G_path = nx.DiGraph()

                                # Fill graph with data
                                G_path.add_edges_from(
                                    solution_paths[solution_path]["path"]
                                )
                                kwargs["logger"].debug(
                                    f"solution path: {solution_paths[solution_path]['path']}"
                                )

                                if 0 in list(G_path.nodes):
                                    # Get all simple paths from node 0 to last node
                                    try:
                                        candidate_paths = list(
                                            nx.all_simple_paths(
                                                G_path, 0, len(candidate_nodes) + 1
                                            )
                                        )

                                        if len(candidate_paths) > 0:
                                            kwargs["logger"].debug(
                                                f"candidate_paths: {candidate_paths[0]}"
                                            )

                                            # Get mapped unitigs in order from the flow network
                                            path_order = []
                                            for path_edge in candidate_paths[0]:
                                                if not (
                                                    path_edge == 0
                                                    or path_edge
                                                    == len(candidate_nodes) + 1
                                                ):
                                                    path_order.append(
                                                        edge_list_indices[path_edge]
                                                    )

                                            kwargs["logger"].debug(
                                                f"path_order: {path_order}"
                                            )

                                            # Get the order of unitigs in path
                                            path_string = ""
                                            total_length = 0

                                            previous_edge = 0

                                            for nodeid in range(len(path_order)):
                                                node = path_order[nodeid]
                                                unitig_name = node[:-1]

                                                if node.endswith("+"):
                                                    unitig_seq = str(
                                                        kwargs["graph_unitigs"][
                                                            unitig_name
                                                        ]
                                                    )
                                                else:
                                                    unitig_seq = str(
                                                        kwargs["graph_unitigs"][
                                                            unitig_name
                                                        ].reverse_complement()
                                                    )

                                                # If first node in path
                                                if previous_edge == 0:
                                                    path_string += unitig_seq
                                                    total_length += len(unitig_seq)

                                                else:
                                                    trimmed_seq = unitig_seq[
                                                        kwargs["link_overlap"][
                                                            (previous_edge, node)
                                                        ] :
                                                    ]
                                                    path_string += trimmed_seq
                                                    total_length += len(trimmed_seq)

                                                previous_edge = node

                                            # Create GenomePath object with path details
                                            genome_path = GenomePath(
                                                id=f"virus_comp_{my_count}_cycle_{cycle_number}",
                                                bubble_case=case_name,
                                                node_order=[x for x in path_order],
                                                node_id_order=[
                                                    kwargs["unitig_names_rev"][x[:-1]]
                                                    for x in path_order
                                                ],
                                                path=path_string,
                                                coverage=int(coverage_val),
                                                length=total_length,
                                                gc=(
                                                    path_string.count("G")
                                                    + path_string.count("C")
                                                )
                                                / len(path_string)
                                                * 100,
                                            )
                                            my_genomic_paths.append(genome_path)
                                            kwargs["logger"].debug(
                                                f"total_length: {total_length}"
                                            )

                                            cycle_number += 1

                                    except nx.exception.NodeNotFound:
                                        kwargs["logger"].debug(
                                            f"Could not resolve a continuous path."
                                        )

                        kwargs["logger"].debug(
                            f"Number of paths selected: {cycle_number - 1}"
                        )

                        if cycle_number > 1:
                            results["resolved_components"].add(my_count)
                            results["resolved_linear"].add(my_count)
                            results["case3_resolved"].add(my_count)

                    else:
                        kwargs["logger"].debug(f"No paths detected")
                        continue

        # Case 1 components - single unitigs
        elif len(candidate_nodes) == 1:
            unitig_name = kwargs["unitig_names"][candidate_nodes[0]]

            results["case1_found"].add(my_count)

            if unitig_name in kwargs["self_looped_nodes"]:
                case_name = "case1_circular"
            else:
                case_name = "case1_linear"

            results["resolved_edges"].add(candidate_nodes[0])
            comp_resolved_edges.add(candidate_nodes[0])

            path_string = str(kwargs["graph_unitigs"][unitig_name])

            cycle_number = 1

            # Create GenomePath object with path details
            genome_path = GenomePath(
                id=f"virus_comp_{my_count}_cycle_{cycle_number}",
                bubble_case=case_name,
                node_order=[kwargs["unitig_names"][candidate_nodes[0]]],
                node_id_order=[candidate_nodes[0]],
                path=path_string,
                coverage=int(kwargs["unitig_coverages"][unitig_name]),
                length=len(kwargs["graph_unitigs"][unitig_name]),
                gc=(path_string.count("G") + path_string.count("C"))
                / len(path_string)
                * 100,
            )
            my_genomic_paths.append(genome_path)
            results["resolved_components"].add(my_count)
            results["single_unitigs"].add(my_count)
            results["case1_resolved"].add(my_count)

            results["virus_like_edges"] = results["virus_like_edges"].union(
                set(candidate_nodes)
            )

        # Record final paths for the component
        # ----------------------------------------------------------------------

        # Order resolved paths in descending order of length
        my_genomic_paths.sort(key=lambda x: (x.length, x.coverage), reverse=True)

        final_genomic_paths = []
        visited_nodes = set()
        comp_resolved_paths = set()

        n_paths = 0

        if len(my_genomic_paths) > 0:
            # Get the degree of the component
            graph_degree = kwargs["assembly_graph"].degree(original_candidate_nodes)

            path_lengths = []
            path_coverages = []

            largest_length = my_genomic_paths[0].length

            # Filter genomic paths
            for genomic_path in my_genomic_paths:
                passed = False

                if genomic_path.length > largest_length * kwargs["LEN_THRESHOLD"]:
                    passed = True

                path_node_order_string = ",".join(genomic_path.node_order)

                if path_node_order_string in comp_resolved_paths:
                    passed = False

                if passed:
                    kwargs["logger"].debug(
                        f"{genomic_path.id}\t{genomic_path.length}\t{genomic_path.coverage}"
                    )
                    kwargs["logger"].debug(f"{genomic_path.node_order}")
                    path_lengths.append(genomic_path.length)
                    path_coverages.append(genomic_path.coverage)
                    final_genomic_paths.append(genomic_path)
                    visited_nodes = visited_nodes.union(set(genomic_path.node_order))
                    comp_resolved_paths.add(path_node_order_string)
                    n_paths += 1

                    for path_node in genomic_path.node_id_order:
                        comp_resolved_edges.add(path_node)

            frac_unitigs = len(visited_nodes) / len(original_candidate_nodes)

            results["resolved_edges"] = results["resolved_edges"].union(
                comp_resolved_edges
            )

            kwargs["logger"].debug(f"frac_unitigs: {frac_unitigs}")

            # Filter components
            if (
                len(final_genomic_paths) > 1
                and len(in_degree) > 0
                and len(out_degree) > 0
            ):
                coverage_frac = (
                    max(path_coverages) / min(path_coverages)
                    if min(path_coverages) > 0
                    else 1
                )

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
                    coverage_frac,
                    frac_unitigs,
                )
                results["all_components"].append(genome_comp)

            if len(final_genomic_paths) > 0:
                results["resolved_components"].add(my_count)
                results["all_resolved_paths"] += final_genomic_paths
                component_elapsed_time = time.time() - component_time_start
                kwargs["logger"].debug(
                    f"Elapsed time to resolve component {my_count} with {len(original_candidate_nodes)} nodes: {component_elapsed_time} seconds"
                )

        else:
            # single unitigs
            for genomic_path in my_genomic_paths:
                final_genomic_paths.append(genomic_path)
                results["all_resolved_paths"].append(genomic_path)
                kwargs["logger"].debug(f"{genomic_path.id}\t{genomic_path.length}")
                results["resolved_components"].add(my_count)

        # Add the paths for writing
        results["genome_path_sets"].add(tuple(final_genomic_paths))

    results_queue.put(results)


def main(**kwargs):
    # Setup logger
    # ----------------------------------------------------------------------

    kwargs["logger"] = setup_logging(**kwargs)

    kwargs["logger"].info(
        "Welcome to Reneo: Unraveling Viral Genomes from Metagenomes."
    )
    kwargs["logger"].info(f"Input arguments: ")
    kwargs["logger"].info(f"Assembly graph file: {kwargs['graph']}")
    kwargs["logger"].info(f"Unitig coverage file: {kwargs['coverage']}")
    kwargs["logger"].info(f"BAM files path: {kwargs['bampath']}")
    kwargs["logger"].info(f"Unitig .hmmout file: {kwargs['hmmout']}")
    kwargs["logger"].info(f"Unitig vog annotations file: {kwargs['vogs']}")
    kwargs["logger"].info(
        f"Minimum length of unitigs to consider: {kwargs['minlength']}"
    )
    kwargs["logger"].info(f"Minimum coverage of paths to output: {kwargs['mincov']}")
    kwargs["logger"].info(
        f"Maximum unitig count to consider a component: {kwargs['compcount']}"
    )
    kwargs["logger"].info(
        f"Maximum number of paths to resolve for a component: {kwargs['maxpaths']}"
    )
    kwargs["logger"].info(
        f"Length threshold to consider single copy marker genes: {kwargs['mgfrac']}"
    )
    kwargs["logger"].info(f"Maximum e-value for vog annotations: {kwargs['evalue']}")
    kwargs["logger"].info(
        f"Minimum hmm score for vog annotations: {kwargs['hmmscore']}"
    )
    kwargs["logger"].info(
        f"Minimum number of vogs to consider a component: {kwargs['nvogs']}"
    )
    kwargs["logger"].info(
        f"Coverage tolerance for extending subpaths: {kwargs['covtol']}"
    )
    kwargs["logger"].info(
        f"Coverage multipler for flow interval modelling: {kwargs['alpha']}"
    )
    kwargs["logger"].info(f"Number of threads to use: {kwargs['nthreads']}")
    kwargs["logger"].info(f"Output folder: {kwargs['output']}")

    start_time = time.time()

    # Init files
    # ----------------------------------------------------------------------
    init_files(kwargs["output"])

    # Get assembly graph
    # ----------------------------------------------------------------------
    (
        kwargs["assembly_graph"],
        kwargs["oriented_links"],
        kwargs["link_overlap"],
        kwargs["unitig_names"],
        kwargs["unitig_names_rev"],
        kwargs["graph_unitigs"],
        kwargs["self_looped_nodes"],
        kwargs["edges_lengths"],
    ) = edge_graph_utils.build_assembly_graph(kwargs["graph"])

    kwargs["logger"].info(
        f"Total number of vertices in the assembly graph: {len(kwargs['assembly_graph'].vs)}"
    )
    kwargs["logger"].info(
        f"Total number of links in the assembly graph: {len(kwargs['assembly_graph'].es)}"
    )

    # Get single unitigs
    # ----------------------------------------------------------------------
    kwargs["circular"] = edge_graph_utils.get_circular(
        kwargs["self_looped_nodes"], kwargs["graph_unitigs"]
    )

    # Get unitigs with bacterial single copy marker genes
    # ----------------------------------------------------------------------
    if kwargs["hmmout"]:
        kwargs["smg_unitigs"] = gene_utils.get_smg_unitigs(
            kwargs["hmmout"], kwargs["mgfrac"]
        )
    else:
        kwargs["smg_unitigs"] = set()

    # Get unitigs with PHROGs
    # ----------------------------------------------------------------------
    if kwargs["vogs"]:
        kwargs["unitig_vogs"] = gene_utils.get_vog_unitigs(
            kwargs["vogs"], kwargs["evalue"], kwargs["hmmscore"]
        )
    else:
        kwargs["unitig_vogs"] = kwargs["graph_unitigs"]

    # Get components with viral genes
    # ----------------------------------------------------------------------
    kwargs["pruned_vs"], kwargs["comp_vogs"] = component_utils.get_components(**kwargs)
    kwargs["logger"].info(
        f"Total number of components found: {len(kwargs['pruned_vs'])}"
    )

    # Get unitig and junction pe coverages
    # ----------------------------------------------------------------------
    kwargs["logger"].info("Getting unitig coverage")
    kwargs["unitig_coverages"] = get_unitig_coverage(kwargs["coverage"])

    kwargs["logger"].info("Getting junction pe coverage")
    with open(kwargs["pickle_file"], "rb") as handle:
        kwargs["junction_pe_coverage"] = pickle.load(handle)

    # Set up worker queues
    component_queue = queue.Queue()
    results_queue = queue.Queue()

    # Populate worker queue
    for my_count in kwargs["pruned_vs"]:
        component_queue.put(my_count)

    # Send finish signals to queue for each worker
    for _ in range(kwargs["nthreads"]):
        component_queue.put(None)

    # Set up multithreading
    worker_threads = []
    for _ in range(kwargs["nthreads"]):
        t = threading.Thread(
            target=worker_resolve_components,
            args=(
                component_queue,
                results_queue,
            ),
            kwargs=kwargs,
        )
        t.start()
        worker_threads.append(t)

    # Wait for workers to finish
    for t in worker_threads:
        t.join()

    # Combine the results from all workers
    results = results_dict()

    while not results_queue.empty():
        r = (
            results_queue.get()
        )  # Dequeue (get and remove) the element from the front of the queue
        results = merge_results(results, r)

    # Get unresolved edges
    results["unresolved_virus_like_edges"] = results["all_virus_like_edges"].difference(
        results["resolved_edges"]
    )

    # write all the final genomic paths
    for final_genomic_paths in results["genome_path_sets"]:
        write_path(final_genomic_paths, kwargs["output"])
        if kwargs["genomes_folder"]:
            write_path_fasta(
                final_genomic_paths, f"{kwargs['output']}/resolved_viruses"
            )

    # Log final summary information
    # ----------------------------------------------------------------------
    kwargs["logger"].info(
        f"Total number of cyclic components found: {len(results['cycle_components'])}"
    )
    kwargs["logger"].info(
        f"Total number of cyclic components resolved: {len(results['resolved_cyclic'])}"
    )
    kwargs["logger"].info(
        f"Single unitigs identified: {len(results['single_unitigs'])}"
    )
    kwargs["logger"].info(
        f"Total number of linear components found: {len(results['linear_components'])}"
    )
    kwargs["logger"].info(
        f"Total number of linear components resolved: {len(results['resolved_linear'])}"
    )
    kwargs["logger"].info(
        f"Total number of cyclic components found including single unitigs: {len(results['cycle_components']) + len(results['single_unitigs'])}"
    )
    kwargs["logger"].info(
        f"Total number of components resolved: {len(results['single_unitigs'])+len(results['resolved_cyclic'])+len(results['resolved_linear'])}"
    )
    kwargs["logger"].info(
        f"Case 1 (resolved/found): {len(results['case1_resolved'])}/{len(results['case1_found'])}"
    )
    kwargs["logger"].info(
        f"Case 2 (resolved/found): {len(results['case2_resolved'])}/{len(results['case2_found'])}"
    )
    kwargs["logger"].info(
        f"Case 3 (resolved/found): {len(results['case3_resolved'])}/{len(results['case3_found'])}"
    )
    kwargs["logger"].info(
        f"Total number of genomes resolved: {len(results['all_resolved_paths'])}"
    )

    if len(results["all_resolved_paths"]) == 0:
        kwargs["logger"].info(f"No genomes were resolved.")
    else:
        kwargs["logger"].info(
            f"Resolved genomes can be found in {kwargs['output']}/resolved_paths.fasta"
        )

    # Write edges to file
    # ----------------------------------------------------------------------
    if kwargs["unitigs"]:
        write_unitigs(
            results["virus_like_edges"],
            kwargs["unitig_names"],
            kwargs["graph_unitigs"],
            "virus_like_edges",
            kwargs["output"],
        )
        write_unitigs(
            results["all_virus_like_edges"],
            kwargs["unitig_names"],
            kwargs["graph_unitigs"],
            "all_virus_like_edges",
            kwargs["output"],
        )

    write_unitigs(
        results["unresolved_virus_like_edges"],
        kwargs["unitig_names"],
        kwargs["graph_unitigs"],
        "unresolved_virus_like_edges",
        kwargs["output"],
    )

    write_unitigs(
        results["resolved_edges"],
        kwargs["unitig_names"],
        kwargs["graph_unitigs"],
        "resolved_edges",
        kwargs["output"],
    )

    # Record path information
    # ----------------------------------------------------------------------

    filename = write_res_genome_info(results["all_resolved_paths"], kwargs["output"])
    if len(results["all_resolved_paths"]) > 0:
        kwargs["logger"].info(
            f"Resolved genome information can be found in {kwargs['output']}/{filename}"
        )

    # Record component information
    # ----------------------------------------------------------------------

    filename = write_component_info(results["all_components"], kwargs["output"])
    if len(results["all_components"]) > 0:
        kwargs["logger"].info(
            f"Resolved component information can be found in {kwargs['output']}/{filename}"
        )

    filename = write_component_vog_info(
        results["resolved_components"], kwargs["comp_vogs"], kwargs["output"]
    )
    if len(results["resolved_components"]) > 0:
        kwargs["logger"].info(
            f"PHROGs found in resolved components can be found in {kwargs['output']}/{filename}"
        )

    # Get elapsed time
    # ----------------------------------------------------------------------

    # Determine elapsed time
    elapsed_time = time.time() - start_time

    # Print elapsed time for the process
    kwargs["logger"].info(f"Elapsed time: {elapsed_time} seconds")

    # Exit program
    # ----------------------------------------------------------------------

    kwargs["logger"].info("Thank you for using Reneo!")


if __name__ == "__main__":
    main(
        graph=snakemake.input.graph,
        coverage=snakemake.input.coverage,
        pickle_file=snakemake.input.pickle,
        bampath=snakemake.params.bampath,
        genomes_folder=snakemake.params.genomes_folder,
        unitigs=snakemake.params.unitigs,
        hmmout=snakemake.params.hmmout,
        vogs=snakemake.params.vogs,
        minlength=int(snakemake.params.minlength),
        mincov=int(snakemake.params.mincov),
        compcount=int(snakemake.params.compcount),
        maxpaths=int(snakemake.params.maxpaths),
        mgfrac=float(snakemake.params.mgfrac),
        evalue=float(snakemake.params.evalue),
        hmmscore=float(snakemake.params.hmmscore),
        nvogs=int(snakemake.params.nvogs),
        covtol=float(snakemake.params.covtol),
        alpha=float(snakemake.params.alpha),
        output=snakemake.params.output,
        nthreads=snakemake.threads,
        log=snakemake.log.stderr,
        MAX_VAL=sys.maxsize,
        LEN_THRESHOLD=0.95,
    )
