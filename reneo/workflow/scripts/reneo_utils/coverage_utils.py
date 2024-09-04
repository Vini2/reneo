#!/usr/bin/env python3

import glob
import os
import pickle
import queue
import threading
from collections import defaultdict

import pysam

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.5.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


def get_unitig_coverage(coverage):
    """
    Get coverage values of unitigs
    """

    unitig_coverages = {}

    with open(coverage, "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("Contig"):
                strings = line.strip().split()

                unitig_name = strings[0]

                coverage_sum = sum([float(x) for x in strings[1:]])

                unitig_coverages[unitig_name] = coverage_sum

    return unitig_coverages


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])

    for read in bam.fetch(region=region_string):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

    return read_dict


def find_links_in_bam(bam_queue, result_queue):
    link_counts = defaultdict(int)
    while True:
        bam_file = bam_queue.get()
        if bam_file is None:
            break
        bam = pysam.AlignmentFile(bam_file, "rb")
        read_pairs = read_pair_generator(bam)
        for read1, read2 in read_pairs:
            if read1.reference_name != read2.reference_name:
                link_counts[(read1.reference_name, read2.reference_name)] += 1
    result_queue.put(link_counts)


def get_junction_pe_coverage(bam_path, output, nthreads):
    """
    Get number of paired end reads supporting a junction
    """

    bam_queue = queue.Queue()
    result_queue = queue.Queue()

    if os.path.isfile(f"{output}/junction_pe_coverage.pickle"):
        with open(f"{output}/junction_pe_coverage.pickle", "rb") as handle:
            link_counts = pickle.load(handle)

    else:
        bam_files = glob.glob(bam_path + "/*.bam")

        # populate queue
        for bam_file in bam_files:
            bam_queue.put(bam_file)

        # send finish signal for workers and spawn workers
        threads = []
        for _ in range(nthreads):
            bam_queue.put(None)
            thread = threading.Thread(
                target=find_links_in_bam, args=(bam_queue, result_queue)
            )
            threads.append(thread)
            thread.start()

        # join workers
        for thread in threads:
            thread.join()

        # combine results
        link_counts = defaultdict(int)
        while not result_queue.empty():
            links = result_queue.get()
            for ctgs, count in links.items():
                link_counts[ctgs] += count

        with open(f"{output}/junction_pe_coverage.pickle", "wb") as handle:
            pickle.dump(link_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return link_counts


def get_graph_spanning_reads(gaf_path, output):
    """
    Get number of reads spanning across a junction
    """

    junction_reads = defaultdict(int)

    if os.path.isfile(f"{output}/graph_spanning_reads.pickle"):
        with open(f"{output}/graph_spanning_reads.pickle", "rb") as handle:
            junction_reads = pickle.load(handle)

    else:
        gaf_files = glob.glob(gaf_path + "/*.gaf")

        for gaf_file in gaf_files:
            with open(gaf_file, "r") as myfile:
                for line in myfile.readlines():
                    strings = line.strip().split("\t")

                    if strings[5].count(">") == 2:
                        edges = strings[5].split(">")[1:]
                        junction_reads[(edges[0], edges[1])] += 1

                    elif strings[5].count("<") == 2:
                        edges = strings[5].split("<")[1:]
                        junction_reads[(edges[1], edges[0])] += 1

        with open(f"{output}/graph_spanning_reads.pickle", "wb") as handle:
            pickle.dump(junction_reads, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return junction_reads
