#!/usr/bin/env python3

import copy
import logging
from collections import defaultdict

from Bio import SeqIO
from agtools.core.unitig_graph import UnitigGraph
from igraph import *

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2026, Reneo Project"
__license__ = "MIT"
__version__ = "0.6.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


# Create logger
logger = logging.getLogger(f"reneo {__version__}")


class BidirectionalError(Exception):
    """Must set a unique value in a BijectiveMap."""

    def __init__(self, value):
        self.value = value
        msg = 'The value "{}" is already in the mapping.'
        super().__init__(msg.format(value))


class BidirectionalMap(dict):
    """Invertible map."""

    def __init__(self, inverse=None):
        if inverse is None:
            inverse = self.__class__(inverse=self)
        self.inverse = inverse

    def __setitem__(self, key, value):
        if value in self.inverse:
            raise BidirectionalError(value)

        self.inverse._set_item(value, key)
        self._set_item(key, value)

    def __delitem__(self, key):
        self.inverse._del_item(self[key])
        self._del_item(key)

    def _del_item(self, key):
        super().__delitem__(key)

    def _set_item(self, key, value):
        super().__setitem__(key, value)


def get_unitig_lengths(edge_file):
    """
    Get length of the unitigs
    """

    unitig_lengths = {}

    for index, record in enumerate(SeqIO.parse(edge_file, "fasta")):
        unitig_lengths[record.id] = len(record.seq)

    return unitig_lengths


def get_name_based_links(unitig_graph, contig_names):
    """
    Return oriented links and overlaps keyed by segment names.
    """

    oriented_links = defaultdict(lambda: defaultdict(list))
    link_overlap = defaultdict(int)

    for from_segment_id, linked_segments in unitig_graph.oriented_links.items():
        from_segment_name = contig_names[from_segment_id]
        for to_segment_id, orientations in linked_segments.items():
            to_segment_name = contig_names[to_segment_id]
            for from_orientation, to_orientation in sorted(orientations):
                oriented_links[from_segment_name][to_segment_name].append(
                    (from_orientation, to_orientation)
                )

    for (
        from_segment_id,
        from_orientation,
        to_segment_id,
        to_orientation,
    ), overlap in unitig_graph.link_overlap.items():
        from_segment_name = contig_names[from_segment_id]
        to_segment_name = contig_names[to_segment_id]
        link_overlap[
            (
                f"{from_segment_name}{from_orientation}",
                f"{to_segment_name}{to_orientation}",
            )
        ] = overlap

    return oriented_links, link_overlap


def build_assembly_graph(assembly_graph_file):
    """
    Build the assembly graph
    """

    unitig_graph = UnitigGraph.from_gfa(assembly_graph_file)

    contig_names = BidirectionalMap()
    for segment_id, segment_name in enumerate(unitig_graph.segment_names):
        contig_names[segment_id] = segment_name

    # Get reverse mapping of contig identifiers
    contig_names_rev = contig_names.inverse

    assembly_graph = unitig_graph.graph.copy()
    for i in range(unitig_graph.vcount):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["name"] = contig_names[i]
        assembly_graph.vs[i]["label"] = contig_names[i] + "\nID:" + str(i)

    oriented_links, link_overlap = get_name_based_links(unitig_graph, contig_names)
    graph_contigs = {
        segment_name: unitig_graph.get_segment_sequence(segment_name)
        for segment_name in unitig_graph.segment_names
    }
    self_looped_nodes = [
        contig_names[segment_id] for segment_id in unitig_graph.self_loops
    ]
    edges_lengths = dict(unitig_graph.segment_lengths)

    return (
        assembly_graph,
        oriented_links,
        link_overlap,
        contig_names,
        contig_names_rev,
        graph_contigs,
        self_looped_nodes,
        edges_lengths,
    )


def get_circular(self_looped_nodes, graph_unitigs):
    """
    Get circular unitigs
    """

    circular = {}

    for unitig in self_looped_nodes:
        circular[unitig] = len(str(graph_unitigs[unitig]))

    # with open(paths, "r") as myfile:

    #     for line in myfile.readlines():
    #         if not line.startswith("#"):
    #             strings = line.strip().split()

    #             if strings[3] == "Y":
    #                 contig_name = strings[0].replace("contig", "edge")
    #                 contig_length = int(strings[1])
    #                 circular[contig_name] = contig_length

    return circular


def remove_dead_ends(G_edge):
    """
    Remove dead-ends from the component
    """

    new_G = copy.deepcopy(G_edge)

    has_dead_ends = True

    dead_ends_to_remove = []

    while has_dead_ends:
        to_remove = []

        for node in list(new_G.nodes):
            if not (new_G.in_degree(node) > 0 and new_G.out_degree()(node)) > 0:
                to_remove.append(node)

        if len(to_remove) > 0:
            new_G.remove_nodes_from(to_remove)
            logger.debug(f"Removing dead-ends: {to_remove}")
        else:
            has_dead_ends = False

        dead_ends_to_remove += to_remove

    return set(dead_ends_to_remove)
