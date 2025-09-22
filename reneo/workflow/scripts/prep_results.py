#!/usr/bin/python3

"""prep_results.py: Format the initial binning result from an existing binning tool.

Format the initial binning result from an existing binning tool in the .csv format
with seq ID and bin ID.

"""


import csv
import logging
import os
import subprocess
import sys

from Bio import SeqIO


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"


def main():

    # Get arguments
    # -----------------------

    bins_folder = snakemake.params.bins
    output_path = snakemake.params.output

    # Setup logger
    # ----------------------------------------------------------------------
    logging.basicConfig(
        filename=snakemake.log.stderr,
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.captureWarnings(True)

    logger = logging.getLogger("prep_results")

    # Get list of files in the folder path of binning result.
    files = os.listdir(bins_folder)


    # Check if folder path of binning result is empty.
    # ---------------------------------------------------
    if len(files) == 0:
        logger.error(
            "Folder containing the initial binning result is empty. Please enter a valid path to the folder containing the initial binning result."
        )
        logger.info("Exiting prep_results.py... Bye...!")
        sys.exit(1)


    # Check if binning result folder contains fasta files.
    # ---------------------------------------------------
    isFasta = False
    for myfile in files:
        if myfile.lower().endswith((".fasta", ".fa", ".fna")):
            isFasta = True

    if not isFasta:
        logger.error(
            "Make sure the folder containing the initial binning result contains fasta files (.fasta, .fa or .fna)."
        )
        logger.info("Exiting prep_results.py... Bye...!")
        sys.exit(1)

    # Format binning results.
    # ---------------------------------------------------

    logger.info("Formatting initial binning results")

    seq_bins = []

    for bin_file in files:
        if bin_file.lower().endswith((".fasta", ".fa", ".fna")):
            fasta_file = bins_folder + bin_file
            for record in SeqIO.parse(fasta_file, "fasta"):
                line = [record.id, str(bin_file)]
                seq_bins.append(line)


    # Write binning results to output file.
    # ---------------------------------------------------

    logger.info("Writing initial binning results to output file")

    with open(output_path, mode="w") as bins_file:
        seq_writer = csv.writer(
            bins_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )

        seq_writer.writerow(["unitig_name", "bin_name"])

        for row in seq_bins:
            seq_writer.writerow(row)

    logger.info(f"Formatted initial binning results can be found at {bins_file.name}")


    # Exit program
    # --------------

    logger.info("Thank you for using prep_results!")

if __name__ == "__main__":
    main()