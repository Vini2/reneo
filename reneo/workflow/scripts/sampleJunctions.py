import logging
import pickle
import sys
from collections import defaultdict

import pysam
from reneo_utils.coverage_utils import read_pair_generator

__author__ = "Michael Roach"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.5.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


def setup_logging(**kwargs):
    """Setup logging"""
    logging.basicConfig(
        filename=kwargs["log"],
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.captureWarnings(True)
    logger = logging.getLogger(__version__)
    return logger


def find_links_in_bam(**kwargs):
    """Return a dictionary of contig pairs and their PE junction counts"""
    link_counts = defaultdict(int)

    bam_file = pysam.AlignmentFile(kwargs["bam"], "rb")
    read_pairs = read_pair_generator(bam_file)

    for read1, read2 in read_pairs:
        if read1.reference_name != read2.reference_name:
            link_counts[(read1.reference_name, read2.reference_name)] += 1

    return link_counts


def pickle_and_dump(**kwargs):
    """Write a pickle file of the dictionary of contig pairs and their PE junction counts"""
    with open(kwargs["out"], "wb") as handle:
        pickle.dump(kwargs["links"], handle, protocol=pickle.HIGHEST_PROTOCOL)


def main(**kwargs):
    logger = setup_logging(**kwargs)

    logger.info(f"Finding junctions for {kwargs['smpl']}")
    kwargs["links"] = find_links_in_bam(**kwargs)

    logger.info("Writing junctions to pickle file")
    pickle_and_dump(**kwargs)

    logger.info("Done!")
    sys.exit(0)


if __name__ == "__main__":
    main(
        bam=snakemake.input.bam,
        out=snakemake.output.pkl,
        smpl=snakemake.wildcards.sample,
        log=snakemake.log.stderr,
    )
