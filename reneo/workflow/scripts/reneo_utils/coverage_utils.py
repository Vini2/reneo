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


def get_read_orientation(read):
    """
    Return the reference strand used by a read alignment.
    """

    return "-" if read.is_reverse else "+"


def get_opposite_orientation(orientation):
    """
    Return the opposite orientation sign.
    """

    return "-" if orientation == "+" else "+"


def add_oriented_pair_support_from_fields(
    link_counts, read_name, read_orientation, mate_name, mate_orientation
):
    """
    Count oriented adjacency from a read record and its mate fields.
    """

    left = f"{read_name}{read_orientation}"
    right = f"{mate_name}{get_opposite_orientation(mate_orientation)}"
    link_counts[(left, right)] += 1

    rc_left = f"{mate_name}{mate_orientation}"
    rc_right = f"{read_name}{get_opposite_orientation(read_orientation)}"
    link_counts[(rc_left, rc_right)] += 1


def add_oriented_spanning_read_support(link_counts, left_read, right_read):
    """
    Count the oriented adjacency implied by two alignments of the same read.
    """

    left_orientation = get_read_orientation(left_read)
    right_orientation = get_read_orientation(right_read)

    left = f"{left_read.reference_name}{left_orientation}"
    right = f"{right_read.reference_name}{right_orientation}"
    link_counts[(left, right)] += 1

    rc_left = f"{right_read.reference_name}{get_opposite_orientation(right_orientation)}"
    rc_right = f"{left_read.reference_name}{get_opposite_orientation(left_orientation)}"
    link_counts[(rc_left, rc_right)] += 1


def add_endpoint_support(endpoint_counts, left, right, target_contigs):
    left_name = left[:-1]
    right_name = right[:-1]

    if left_name in target_contigs:
        endpoint_counts[left] += 1
    if right_name in target_contigs:
        endpoint_counts[right] += 1


def add_endpoint_pair_support_from_fields(
    endpoint_counts,
    read_name,
    read_orientation,
    mate_name,
    mate_orientation,
    target_contigs,
):
    left = f"{read_name}{read_orientation}"
    right = f"{mate_name}{get_opposite_orientation(mate_orientation)}"
    add_endpoint_support(endpoint_counts, left, right, target_contigs)

    rc_left = f"{mate_name}{mate_orientation}"
    rc_right = f"{read_name}{get_opposite_orientation(read_orientation)}"
    add_endpoint_support(endpoint_counts, rc_left, rc_right, target_contigs)


def add_endpoint_spanning_read_support(
    endpoint_counts, left_read, right_read, target_contigs
):
    left_orientation = get_read_orientation(left_read)
    right_orientation = get_read_orientation(right_read)

    left = f"{left_read.reference_name}{left_orientation}"
    right = f"{right_read.reference_name}{right_orientation}"
    add_endpoint_support(endpoint_counts, left, right, target_contigs)

    rc_left = f"{right_read.reference_name}{get_opposite_orientation(right_orientation)}"
    rc_right = f"{left_read.reference_name}{get_opposite_orientation(left_orientation)}"
    add_endpoint_support(endpoint_counts, rc_left, rc_right, target_contigs)


def add_endpoint_sa_tag_support(endpoint_counts, read, target_contigs):
    if not read.has_tag("SA"):
        return

    for alignment in read.get_tag("SA").split(";"):
        if alignment == "":
            continue

        fields = alignment.split(",")
        if len(fields) < 3:
            continue

        sa_reference_name = fields[0]
        sa_orientation = fields[2]

        if sa_reference_name == read.reference_name:
            continue

        left = f"{read.reference_name}{get_read_orientation(read)}"
        right = f"{sa_reference_name}{sa_orientation}"
        add_endpoint_support(endpoint_counts, left, right, target_contigs)

        rc_left = f"{sa_reference_name}{get_opposite_orientation(sa_orientation)}"
        rc_right = (
            f"{read.reference_name}{get_opposite_orientation(get_read_orientation(read))}"
        )
        add_endpoint_support(endpoint_counts, rc_left, rc_right, target_contigs)


def find_oriented_external_endpoint_read_support_in_bam(
    bam_queue, result_queue, target_contigs
):
    endpoint_counts = defaultdict(int)
    target_contigs = set(target_contigs)

    while True:
        bam_file = bam_queue.get()
        if bam_file is None:
            break

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            references = set(bam.references)
            reads = defaultdict(list)

            for contig in sorted(target_contigs.intersection(references)):
                for read in bam.fetch(contig):
                    if read.is_unmapped or read.is_secondary:
                        continue

                    add_endpoint_sa_tag_support(endpoint_counts, read, target_contigs)

                    if (
                        not read.is_supplementary
                        and not read.mate_is_unmapped
                        and read.reference_name != read.next_reference_name
                    ):
                        mate_orientation = "-" if read.mate_is_reverse else "+"
                        add_endpoint_pair_support_from_fields(
                            endpoint_counts,
                            read.reference_name,
                            get_read_orientation(read),
                            read.next_reference_name,
                            mate_orientation,
                            target_contigs,
                        )

                    read_end = 1 if read.is_read1 else 2 if read.is_read2 else 0
                    reads[(read.query_name, read_end)].append(read)

            for alignments in reads.values():
                if len(alignments) < 2:
                    continue

                alignments = sorted(
                    alignments,
                    key=lambda read: (
                        read.query_alignment_start,
                        read.query_alignment_end,
                        read.reference_name,
                        read.reference_start,
                    ),
                )

                seen_pairs = set()
                for i in range(len(alignments) - 1):
                    left = alignments[i]
                    right = alignments[i + 1]

                    if left.reference_name == right.reference_name:
                        continue

                    read_pair = (
                        left.reference_name,
                        left.reference_start,
                        right.reference_name,
                        right.reference_start,
                    )
                    if read_pair in seen_pairs:
                        continue

                    seen_pairs.add(read_pair)
                    add_endpoint_spanning_read_support(
                        endpoint_counts, left, right, target_contigs
                    )

    result_queue.put(endpoint_counts)


def get_oriented_external_endpoint_read_support(bam_path, output, contigs, nthreads):
    """
    Count strand-aware external PE and split-read support touching selected contig ends.
    """

    contigs = sorted(set(contigs))
    cache_file = f"{output}/oriented_external_endpoint_read_support.pickle"

    if os.path.isfile(cache_file):
        with open(cache_file, "rb") as handle:
            cached = pickle.load(handle)

        if (
            isinstance(cached, dict)
            and cached.get("contigs") == contigs
            and "counts" in cached
        ):
            return cached["counts"]

    bam_queue = queue.Queue()
    result_queue = queue.Queue()
    bam_files = glob.glob(bam_path + "/*.bam")

    for bam_file in bam_files:
        bam_queue.put(bam_file)

    threads = []
    for _ in range(nthreads):
        bam_queue.put(None)
        thread = threading.Thread(
            target=find_oriented_external_endpoint_read_support_in_bam,
            args=(bam_queue, result_queue, contigs),
        )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    endpoint_counts = defaultdict(int)
    while not result_queue.empty():
        counts = result_queue.get()
        for endpoint, count in counts.items():
            endpoint_counts[endpoint] += count

    with open(cache_file, "wb") as handle:
        pickle.dump(
            {"contigs": contigs, "counts": endpoint_counts},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    return endpoint_counts


def find_oriented_links_for_pairs_in_bam(bam_queue, result_queue, target_pairs):
    link_counts = defaultdict(int)
    target_pairs = set(target_pairs)
    target_contigs = set()
    for pair in target_pairs:
        target_contigs.update(pair)

    while True:
        bam_file = bam_queue.get()
        if bam_file is None:
            break

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            references = set(bam.references)
            for contig in sorted(target_contigs.intersection(references)):
                for read in bam.fetch(contig):
                    if (
                        read.is_secondary
                        or read.is_supplementary
                        or read.is_unmapped
                        or read.mate_is_unmapped
                        or not read.is_read1
                        or read.reference_name == read.next_reference_name
                        or read.next_reference_name not in target_contigs
                    ):
                        continue

                    pair = tuple(
                        sorted([read.reference_name, read.next_reference_name])
                    )
                    if pair in target_pairs:
                        mate_orientation = "-" if read.mate_is_reverse else "+"
                        add_oriented_pair_support_from_fields(
                            link_counts,
                            read.reference_name,
                            get_read_orientation(read),
                            read.next_reference_name,
                            mate_orientation,
                        )

    result_queue.put(link_counts)


def find_oriented_spanning_reads_for_pairs_in_bam(
    bam_queue, result_queue, target_pairs
):
    link_counts = defaultdict(int)
    target_pairs = set(target_pairs)
    target_contigs = set()
    for pair in target_pairs:
        target_contigs.update(pair)

    while True:
        bam_file = bam_queue.get()
        if bam_file is None:
            break

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            reads = defaultdict(list)
            references = set(bam.references)

            for contig in sorted(target_contigs.intersection(references)):
                for read in bam.fetch(contig):
                    if read.is_unmapped or read.is_secondary:
                        continue

                    read_end = 1 if read.is_read1 else 2 if read.is_read2 else 0
                    reads[(read.query_name, read_end)].append(read)

            for alignments in reads.values():
                if len(alignments) < 2:
                    continue

                alignments = sorted(
                    alignments,
                    key=lambda read: (
                        read.query_alignment_start,
                        read.query_alignment_end,
                        read.reference_name,
                        read.reference_start,
                    ),
                )

                seen_pairs = set()
                for i in range(len(alignments) - 1):
                    left = alignments[i]
                    right = alignments[i + 1]

                    if left.reference_name == right.reference_name:
                        continue

                    pair = tuple(sorted([left.reference_name, right.reference_name]))
                    if pair not in target_pairs:
                        continue

                    read_pair = (
                        left.reference_name,
                        left.reference_start,
                        right.reference_name,
                        right.reference_start,
                    )
                    if read_pair in seen_pairs:
                        continue

                    seen_pairs.add(read_pair)
                    add_oriented_spanning_read_support(link_counts, left, right)

    result_queue.put(link_counts)


def get_oriented_junction_pe_coverage_for_pairs(
    bam_path, output, target_pairs, nthreads
):
    """
    Get strand-aware PE support for selected unordered contig pairs.
    """

    target_pairs = sorted(set([tuple(sorted(x)) for x in target_pairs]))
    cache_file = f"{output}/oriented_junction_pe_coverage.pickle"

    if os.path.isfile(cache_file):
        with open(cache_file, "rb") as handle:
            cached = pickle.load(handle)

        if (
            isinstance(cached, dict)
            and cached.get("target_pairs") == target_pairs
            and "counts" in cached
        ):
            return cached["counts"]

    bam_queue = queue.Queue()
    result_queue = queue.Queue()
    bam_files = glob.glob(bam_path + "/*.bam")

    for bam_file in bam_files:
        bam_queue.put(bam_file)

    threads = []
    for _ in range(nthreads):
        bam_queue.put(None)
        thread = threading.Thread(
            target=find_oriented_links_for_pairs_in_bam,
            args=(bam_queue, result_queue, target_pairs),
        )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    link_counts = defaultdict(int)
    while not result_queue.empty():
        links = result_queue.get()
        for ctgs, count in links.items():
            link_counts[ctgs] += count

    with open(cache_file, "wb") as handle:
        pickle.dump(
            {"target_pairs": target_pairs, "counts": link_counts},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    return link_counts


def get_oriented_spanning_read_coverage_for_pairs(
    bam_path, output, target_pairs, nthreads
):
    """
    Get strand-aware split-read support for selected unordered contig pairs.
    """

    target_pairs = sorted(set([tuple(sorted(x)) for x in target_pairs]))
    cache_file = f"{output}/oriented_spanning_read_coverage.pickle"

    if os.path.isfile(cache_file):
        with open(cache_file, "rb") as handle:
            cached = pickle.load(handle)

        if (
            isinstance(cached, dict)
            and cached.get("target_pairs") == target_pairs
            and "counts" in cached
        ):
            return cached["counts"]

    bam_queue = queue.Queue()
    result_queue = queue.Queue()
    bam_files = glob.glob(bam_path + "/*.bam")

    for bam_file in bam_files:
        bam_queue.put(bam_file)

    threads = []
    for _ in range(nthreads):
        bam_queue.put(None)
        thread = threading.Thread(
            target=find_oriented_spanning_reads_for_pairs_in_bam,
            args=(bam_queue, result_queue, target_pairs),
        )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    link_counts = defaultdict(int)
    while not result_queue.empty():
        links = result_queue.get()
        for ctgs, count in links.items():
            link_counts[ctgs] += count

    with open(cache_file, "wb") as handle:
        pickle.dump(
            {"target_pairs": target_pairs, "counts": link_counts},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    return link_counts


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
