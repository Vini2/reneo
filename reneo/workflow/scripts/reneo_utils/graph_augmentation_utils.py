#!/usr/bin/env python3

import networkx as nx
from reneo_utils import component_utils, flow_utils
from reneo_utils.coverage_utils import (
    get_component_external_endpoint_read_support,
    get_opposite_orientation,
    get_oriented_external_endpoint_read_support,
    get_oriented_junction_pe_coverage_for_pairs,
    get_oriented_spanning_read_coverage_for_pairs,
)

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2026, Reneo Project"
__license__ = "MIT"
__version__ = "0.6.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


JUNCTION_PE_THRESHOLD = 10


def build_oriented_component_graph(candidate_nodes, **kwargs):
    """
    Build the directed, orientation-aware graph used to classify a component.
    """

    G_edge = nx.DiGraph()
    cycle_edges = {}

    for vertex in candidate_nodes:
        unitig_name = kwargs["unitig_names"][vertex]

        for node in kwargs["oriented_links"][unitig_name]:
            if node not in kwargs["unitig_names_rev"]:
                continue

            if kwargs["unitig_names_rev"][node] not in candidate_nodes:
                continue

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

                min_cov = (
                    min([cov_1, cov_2])
                    if min([cov_1, cov_2]) != 0
                    else max([cov_1, cov_2])
                )

                for edge in kwargs["oriented_links"][unitig_name][node]:
                    cycle_edges[(unitig_name + edge[0], node + edge[1])] = int(min_cov)

    for cedge in cycle_edges:
        G_edge.add_edge(cedge[0], cedge[1], weight=cycle_edges[cedge])

    return G_edge


def get_case3_linear_terminals(**kwargs):
    """
    Find source and sink terminal contigs for incomplete case 3 linear components.
    """

    source_terminals = []
    sink_terminals = []

    for component_id, candidate_nodes in kwargs["pruned_vs"].items():
        if not (
            len(candidate_nodes) > 2 and len(candidate_nodes) <= kwargs["compcount"]
        ):
            continue

        G_edge = build_oriented_component_graph(candidate_nodes, **kwargs)

        if len(G_edge.nodes) == 0:
            continue

        try:
            nx.find_cycle(G_edge, orientation="original")
            continue
        except nx.exception.NetworkXNoCycle:
            pass

        source_candidates, sink_candidates = flow_utils.get_source_sink_linear(
            G_edge,
            kwargs["graph_unitigs"],
            kwargs["self_looped_nodes"],
        )

        if len(source_candidates) == 0 or len(sink_candidates) == 0:
            continue

        for source in source_candidates:
            source_terminals.append(
                {
                    "component_id": component_id,
                    "oriented_node": source,
                    "unitig_name": source[:-1],
                }
            )

        for sink in sink_candidates:
            sink_terminals.append(
                {
                    "component_id": component_id,
                    "oriented_node": sink,
                    "unitig_name": sink[:-1],
                }
            )

    return source_terminals, sink_terminals


def get_isolated_linear_unitigs(**kwargs):
    """
    Find single-contig linear components in the full assembly graph.
    """

    isolated = set()

    for component in kwargs["assembly_graph"].components():
        if len(component) != 1:
            continue

        unitig_name = kwargs["unitig_names"][component[0]]
        if unitig_name not in kwargs["self_looped_nodes"]:
            isolated.add(unitig_name)

    return isolated


def filter_isolated_unitigs_by_junction_support(
    isolated_unitigs, terminal_unitigs, **kwargs
):
    """
    Keep isolated unitigs with enough existing PE support to a case 3 terminal.
    """

    supported_isolated_unitigs = set()

    for contigs, count in kwargs["junction_pe_coverage"].items():
        if count < JUNCTION_PE_THRESHOLD or len(contigs) != 2:
            continue

        contig_1, contig_2 = contigs

        if contig_1 in terminal_unitigs and contig_2 in isolated_unitigs:
            supported_isolated_unitigs.add(contig_2)
        elif contig_2 in terminal_unitigs and contig_1 in isolated_unitigs:
            supported_isolated_unitigs.add(contig_1)

    return supported_isolated_unitigs


def get_candidate_extension_pairs(isolated_unitigs, terminal_unitigs, **kwargs):
    """
    Get unordered contig pairs worth checking for strand-aware support.
    """

    candidate_pairs = set()

    for contigs, count in kwargs["junction_pe_coverage"].items():
        if count < JUNCTION_PE_THRESHOLD or len(contigs) != 2:
            continue

        contig_1, contig_2 = contigs

        if contig_1 in terminal_unitigs and contig_2 in terminal_unitigs:
            candidate_pairs.add(tuple(sorted([contig_1, contig_2])))
        elif contig_1 in terminal_unitigs and contig_2 in isolated_unitigs:
            candidate_pairs.add(tuple(sorted([contig_1, contig_2])))
        elif contig_2 in terminal_unitigs and contig_1 in isolated_unitigs:
            candidate_pairs.add(tuple(sorted([contig_1, contig_2])))

    return candidate_pairs


def is_viral_singleton_candidate(unitig_name, **kwargs):
    if unitig_name in kwargs["smg_unitigs"]:
        return False

    if unitig_name not in kwargs["unitig_vogs"]:
        return False

    if len(kwargs["unitig_vogs"][unitig_name]) < kwargs["nvogs"]:
        return False

    if kwargs["edges_lengths"][unitig_name] <= kwargs["minlength"]:
        return False

    return True


def add_complete_isolated_linear_unitigs(**kwargs):
    """
    Add isolated linear viral unitigs with no external end-linking evidence.
    """

    isolated_unitigs = get_isolated_linear_unitigs(**kwargs)
    augmented_unitigs = set()

    for left, right in kwargs.get("inferred_case3_links", {}):
        augmented_unitigs.add(left[:-1])
        augmented_unitigs.add(right[:-1])

    candidate_unitigs = sorted(
        [
            unitig
            for unitig in isolated_unitigs
            if unitig not in augmented_unitigs
            and is_viral_singleton_candidate(unitig, **kwargs)
        ]
    )

    if len(candidate_unitigs) == 0:
        kwargs.setdefault("complete_isolated_linear_unitigs", {}).clear()
        kwargs["logger"].info(
            "Added 0 complete isolated linear unitigs with no external end-linking evidence as case1 candidates"
        )
        return 0

    endpoint_support = get_oriented_external_endpoint_read_support(
        kwargs["bampath"], kwargs["output"], candidate_unitigs, kwargs["nthreads"]
    )

    next_component_id = max(kwargs["pruned_vs"].keys(), default=-1) + 1
    complete_unitigs = kwargs.setdefault("complete_isolated_linear_unitigs", {})
    complete_unitigs.clear()

    for unitig in candidate_unitigs:
        plus_support = endpoint_support.get(f"{unitig}+", 0)
        minus_support = endpoint_support.get(f"{unitig}-", 0)

        if plus_support + minus_support != 0:
            continue

        component_id = next_component_id
        next_component_id += 1
        vertex_id = kwargs["unitig_names_rev"][unitig]
        kwargs["pruned_vs"][component_id] = [vertex_id]
        kwargs["comp_vogs"][component_id] = kwargs["unitig_vogs"].get(unitig, set())
        complete_unitigs[unitig] = {
            "component_id": component_id,
            "plus_external_support": plus_support,
            "minus_external_support": minus_support,
        }

    kwargs["complete_isolated_linear_unitigs"] = complete_unitigs
    kwargs["logger"].info(
        f"Added {len(complete_unitigs)} complete isolated linear unitigs with no external end-linking evidence as case1 candidates"
    )

    return len(complete_unitigs)


def get_augmented_unitigs(**kwargs):
    augmented_unitigs = set()

    for left, right in kwargs.get("inferred_case3_links", {}):
        augmented_unitigs.add(left[:-1])
        augmented_unitigs.add(right[:-1])

    return augmented_unitigs


def get_unextended_case3_linear_component_terminals(**kwargs):
    """
    Find source/sink terminals for case 3 linear components not touched by augmentation.
    """

    augmented_unitigs = get_augmented_unitigs(**kwargs)
    component_terminals = {}
    contig_component_ids = {}
    terminal_unitigs = set()

    for component_id, candidate_nodes in kwargs["pruned_vs"].items():
        if not (
            len(candidate_nodes) > 2 and len(candidate_nodes) <= kwargs["compcount"]
        ):
            continue

        component_unitigs = set(
            kwargs["unitig_names"][node] for node in candidate_nodes
        )
        if len(component_unitigs.intersection(augmented_unitigs)) > 0:
            continue

        G_edge = build_oriented_component_graph(candidate_nodes, **kwargs)

        if len(G_edge.nodes) == 0:
            continue

        try:
            nx.find_cycle(G_edge, orientation="original")
            continue
        except nx.exception.NetworkXNoCycle:
            pass

        source_candidates, sink_candidates = flow_utils.get_source_sink_linear(
            G_edge,
            kwargs["graph_unitigs"],
            kwargs["self_looped_nodes"],
        )

        terminals = source_candidates + sink_candidates
        if len(terminals) == 0:
            continue

        component_terminals[component_id] = {
            "sources": source_candidates,
            "sinks": sink_candidates,
            "terminals": terminals,
        }
        for unitig in component_unitigs:
            contig_component_ids[unitig] = component_id
        for terminal in terminals:
            terminal_unitigs.add(terminal[:-1])

    return component_terminals, terminal_unitigs, contig_component_ids


def get_unextended_case3_linear_closure_support(component_terminals, **kwargs):
    """
    Count strand-aware PE and spanning-read evidence from sinks back to sources.
    """

    target_pairs = set()
    pair_components = {}

    for component_id, terminal_roles in component_terminals.items():
        for sink in terminal_roles["sinks"]:
            for source in terminal_roles["sources"]:
                sink_unitig = sink[:-1]
                source_unitig = source[:-1]

                if sink_unitig == source_unitig:
                    continue

                target_pair = tuple(sorted([sink_unitig, source_unitig]))
                target_pairs.add(target_pair)
                pair_components[(sink, source)] = component_id

    if len(target_pairs) == 0:
        return {}

    oriented_junction_pe_coverage = get_oriented_junction_pe_coverage_for_pairs(
        kwargs["bampath"],
        kwargs["output"],
        target_pairs,
        kwargs["nthreads"],
    )
    oriented_spanning_read_coverage = get_oriented_spanning_read_coverage_for_pairs(
        kwargs["bampath"],
        kwargs["output"],
        target_pairs,
        kwargs["nthreads"],
    )

    closure_support = {}
    for link, component_id in pair_components.items():
        support = oriented_junction_pe_coverage.get(
            link, 0
        ) + oriented_spanning_read_coverage.get(link, 0)
        if support < JUNCTION_PE_THRESHOLD:
            continue

        if (
            component_id not in closure_support
            or support > closure_support[component_id]["support"]
        ):
            closure_support[component_id] = {
                "link": link,
                "support": support,
                "junction_pe_support": oriented_junction_pe_coverage.get(link, 0),
                "spanning_support": oriented_spanning_read_coverage.get(link, 0),
            }

    return closure_support


def filter_unextended_case3_linear_components_by_endpoint_support(**kwargs):
    """
    Remove unextended case 3 linear components unless all terminal ends lack
    external evidence, or sink-to-source read evidence suggests missed
    circularization.
    """

    (
        component_terminals,
        terminal_unitigs,
        contig_component_ids,
    ) = get_unextended_case3_linear_component_terminals(**kwargs)

    if len(component_terminals) == 0:
        kwargs["logger"].info(
            "Filtered 0 unextended case 3 linear components with external end-linking evidence"
        )
        return 0

    endpoint_support = get_component_external_endpoint_read_support(
        kwargs["bampath"],
        kwargs["output"],
        terminal_unitigs,
        contig_component_ids,
        kwargs["nthreads"],
    )
    closure_support = get_unextended_case3_linear_closure_support(
        component_terminals, **kwargs
    )

    removed_components = []
    complete_components = {}
    circularized_components = {}

    for component_id, terminal_roles in component_terminals.items():
        terminal_support = {
            terminal: endpoint_support.get(terminal, 0)
            for terminal in terminal_roles["terminals"]
        }

        if sum(terminal_support.values()) == 0:
            complete_components[component_id] = terminal_support
            continue

        if component_id in closure_support:
            circularized_components[component_id] = {
                "terminal_support": terminal_support,
                "closure_support": closure_support[component_id],
            }
            continue

        removed_components.append(component_id)

    for component_id in removed_components:
        kwargs["pruned_vs"].pop(component_id, None)
        kwargs["comp_vogs"].pop(component_id, None)

    kwargs["complete_unextended_case3_linear_components"] = complete_components
    kwargs["circularized_unextended_case3_linear_components"] = circularized_components
    kwargs["logger"].info(
        f"Kept {len(complete_components)} unextended case 3 linear components with zero external end-linking evidence"
    )
    kwargs["logger"].info(
        f"Kept {len(circularized_components)} unextended case 3 linear components with sink-to-source read evidence suggesting missed circularization"
    )
    kwargs["logger"].info(
        f"Filtered {len(removed_components)} unextended case 3 linear components with external end-linking evidence"
    )

    return len(removed_components)


def add_inferred_oriented_link(left, right, support, **kwargs):
    """
    Add an inferred, zero-overlap oriented link and its reverse complement.
    """

    left_name = left[:-1]
    right_name = right[:-1]
    left_orientation = left[-1]
    right_orientation = right[-1]

    if left_name == right_name:
        return False

    if right_name not in kwargs["oriented_links"][left_name]:
        kwargs["oriented_links"][left_name][right_name] = []

    added_primary_link = False
    if (
        left_orientation,
        right_orientation,
    ) not in kwargs["oriented_links"][
        left_name
    ][right_name]:
        kwargs["oriented_links"][left_name][right_name].append(
            (left_orientation, right_orientation)
        )
        added_primary_link = True

    kwargs["link_overlap"][(left, right)] = 0

    rc_left = f"{right_name}{get_opposite_orientation(right_orientation)}"
    rc_right = f"{left_name}{get_opposite_orientation(left_orientation)}"

    if left_name not in kwargs["oriented_links"][right_name]:
        kwargs["oriented_links"][right_name][left_name] = []

    rc_orientation = (rc_left[-1], rc_right[-1])
    if rc_orientation not in kwargs["oriented_links"][right_name][left_name]:
        kwargs["oriented_links"][right_name][left_name].append(rc_orientation)

    kwargs["link_overlap"][(rc_left, rc_right)] = 0

    kwargs["junction_pe_coverage"][(left_name, right_name)] = max(
        kwargs["junction_pe_coverage"][(left_name, right_name)], support
    )
    kwargs["junction_pe_coverage"][(right_name, left_name)] = max(
        kwargs["junction_pe_coverage"][(right_name, left_name)], support
    )

    left_id = kwargs["unitig_names_rev"][left_name]
    right_id = kwargs["unitig_names_rev"][right_name]

    added_graph_edge = False
    if not kwargs["assembly_graph"].are_connected(left_id, right_id):
        kwargs["assembly_graph"].add_edge(left_id, right_id)
        added_graph_edge = True

    return added_primary_link or added_graph_edge


def augment_case3_linear_components_once(round_id, **kwargs):
    """
    Run one case 3 linear terminal extension round.
    """

    source_terminals, sink_terminals = get_case3_linear_terminals(**kwargs)
    isolated_unitigs = get_isolated_linear_unitigs(**kwargs)
    terminal_unitigs = set([x["unitig_name"] for x in source_terminals])
    terminal_unitigs.update([x["unitig_name"] for x in sink_terminals])
    isolated_unitigs = filter_isolated_unitigs_by_junction_support(
        isolated_unitigs, terminal_unitigs, **kwargs
    )

    if len(source_terminals) == 0 and len(sink_terminals) == 0:
        kwargs["logger"].info(
            "No incomplete case 3 linear terminals found for graph augmentation"
        )
        return 0

    extension_contigs = set(isolated_unitigs)
    extension_contigs.update(terminal_unitigs)
    candidate_pairs = get_candidate_extension_pairs(
        isolated_unitigs, terminal_unitigs, **kwargs
    )

    kwargs["logger"].info(
        f"Checking case 3 linear graph extension round {round_id} using {len(extension_contigs)} candidate contigs and {len(candidate_pairs)} candidate PE-supported pairs"
    )

    if len(candidate_pairs) == 0:
        kwargs["logger"].info(
            f"Added 0 inferred case 3 graph extension links using junction_pe_threshold={JUNCTION_PE_THRESHOLD}"
        )
        return 0

    oriented_junction_pe_coverage = get_oriented_junction_pe_coverage_for_pairs(
        kwargs["bampath"], kwargs["output"], candidate_pairs, kwargs["nthreads"]
    )
    oriented_spanning_read_coverage = get_oriented_spanning_read_coverage_for_pairs(
        kwargs["bampath"], kwargs["output"], candidate_pairs, kwargs["nthreads"]
    )

    inferred_links = {}
    source_components = {
        source["oriented_node"]: source["component_id"] for source in source_terminals
    }
    sink_components = {
        sink["oriented_node"]: sink["component_id"] for sink in sink_terminals
    }

    observed_links = set(oriented_junction_pe_coverage.keys())
    observed_links.update(oriented_spanning_read_coverage.keys())

    for link in observed_links:
        support = oriented_junction_pe_coverage.get(
            link, 0
        ) + oriented_spanning_read_coverage.get(link, 0)

        if support < JUNCTION_PE_THRESHOLD or len(link) != 2:
            continue

        left, right = link
        left_name = left[:-1]
        right_name = right[:-1]
        link_type = None
        left_endpoint_type = None
        right_endpoint_type = None
        left_component = "NA"
        right_component = "NA"

        if (
            left in sink_components
            and right in source_components
            and sink_components[left] != source_components[right]
        ):
            link_type = "case3_to_case3"
            left_endpoint_type = "case3_sink"
            right_endpoint_type = "case3_source"
            left_component = sink_components[left]
            right_component = source_components[right]
        elif left in sink_components and right_name in isolated_unitigs:
            link_type = "case3_to_isolated"
            left_endpoint_type = "case3_sink"
            right_endpoint_type = "isolated"
            left_component = sink_components[left]
        elif left_name in isolated_unitigs and right in source_components:
            link_type = "isolated_to_case3"
            left_endpoint_type = "isolated"
            right_endpoint_type = "case3_source"
            right_component = source_components[right]

        if link_type is not None and support > inferred_links.get(link, {}).get(
            "support", 0
        ):
            inferred_links[link] = {
                "support": support,
                "link_type": link_type,
                "left_endpoint_type": left_endpoint_type,
                "right_endpoint_type": right_endpoint_type,
                "left_component": left_component,
                "right_component": right_component,
            }

    links_added = 0
    for (left, right), metadata in sorted(inferred_links.items()):
        support = metadata["support"]
        if add_inferred_oriented_link(left, right, support, **kwargs):
            kwargs.setdefault("inferred_case3_links", {})[(left, right)] = {
                "support": support,
                "pe_support": oriented_junction_pe_coverage.get((left, right), 0),
                "spanning_support": oriented_spanning_read_coverage.get(
                    (left, right), 0
                ),
                "round": round_id,
                "link_type": metadata["link_type"],
                "left_endpoint_type": metadata["left_endpoint_type"],
                "right_endpoint_type": metadata["right_endpoint_type"],
                "left_component": metadata["left_component"],
                "right_component": metadata["right_component"],
            }
            kwargs["logger"].info(
                f"Added inferred case 3 extension round {round_id}: {left} -> {right} with {support} strand-aware read links"
            )
            links_added += 1

    kwargs["logger"].info(
        f"Added {links_added} inferred case 3 graph extension links in round {round_id} using junction_pe_threshold={JUNCTION_PE_THRESHOLD}"
    )

    return links_added


def augment_case3_linear_components(**kwargs):
    """
    Iteratively connect case 3 linear terminals to other case 3 terminals
    or isolated contigs until no more reliable extensions are found.
    """

    total_links_added = 0
    round_id = 1

    while True:
        links_added = augment_case3_linear_components_once(round_id, **kwargs)

        if links_added == 0:
            break

        total_links_added += links_added
        kwargs["pruned_vs"], kwargs["comp_vogs"] = component_utils.get_components(
            **kwargs
        )
        kwargs["logger"].info(
            f"Total number of components found after graph augmentation round {round_id}: {len(kwargs['pruned_vs'])}"
        )
        round_id += 1

    kwargs["logger"].info(
        f"Added {total_links_added} inferred case 3 graph extension links across {round_id - 1} augmentation rounds using junction_pe_threshold={JUNCTION_PE_THRESHOLD}"
    )

    return total_links_added, kwargs["pruned_vs"], kwargs["comp_vogs"]
