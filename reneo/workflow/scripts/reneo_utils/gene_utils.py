#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2023, Reneo Project"
__license__ = "MIT"
__version__ = "0.5.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


def get_smg_unitigs(hmmout, mg_frac):
    """
    Get unitigs containing bacterial single-copy marker genes
    """

    # Commands
    # run_FragGeneScan.pl -genome=edges.fasta -out=edges.fasta.frag -complete=0 -train=complete -thread=8 1>edges.fasta.frag.out 2>edges.fasta.frag.err
    # hmmsearch --domtblout edges.fasta.hmmout --cut_tc --cpu 8 /home/mall0133/software/MetaCoAG/metacoag_utils/auxiliary/marker.hmm edges.fasta.frag.faa 1>edges.fasta.hmmout.out 2> edges.fasta.hmmout.err

    smg_unitigs = set()

    unitig_smgs = {}

    with open(hmmout, "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#") and line.startswith("edge_"):
                strings = line.strip().split()

                unitig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                name_strings = unitig.split("_")
                name_strings = name_strings[: len(name_strings) - 3]

                # unitig name
                unitig_name = "_".join(name_strings)

                if mapped_marker_length > marker_gene_length * mg_frac:
                    smg_unitigs.add(unitig_name)

                    if unitig_name not in unitig_smgs:
                        unitig_smgs[unitig_name] = set()
                        unitig_smgs[unitig_name].add(marker_gene)
                    else:
                        unitig_smgs[unitig_name].add(marker_gene)

    return smg_unitigs


def get_vog_unitigs(vogs, e_value, hmm_score):
    """
    Get unitigs containing VOGs
    """

    # Get unitigs containing vogs
    unitig_vogs = {}

    with open(vogs, "r") as myfile:
        for line in myfile.readlines():

            if not line.startswith("#"):

                strings = line.strip().split()

                name = "_".join(strings[0].split("_")[:-1])
                vog_id = strings[2]
                evalue = float(strings[4])
                hmmScore = float(strings[5])

                if evalue < e_value and hmmScore > hmm_score:
                    if name not in unitig_vogs:
                        unitig_vogs[name] = set([vog_id])
                    else:
                        unitig_vogs[name].add(vog_id)

    return unitig_vogs
