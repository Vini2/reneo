"""
Reneo: Unraveling Viral Genomes from Metagenomes.

2023, Vijini Mallawaarachchi

This is the main Snakefile to run reneo.
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")


"""PREFLIGHT CHECKS"""
include: os.path.join("rules", "00_database_preflight.smk")
include: os.path.join("rules", "02_reneo_preflight.smk")


"""TARGETS"""
include: os.path.join("rules", "02_reneo_targets.smk")


"""Target rules"""
target_rules = []

def targetRule(fn):
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn

localrules: all, preprocess, reneo, koverage_tsv, print_stages


"""Run stages"""
@targetRule
rule all:
    input:
        preprocessTargets,
        reneoTargets


@targetRule
rule preprocess:
    input:
        preprocessTargets


@targetRule
rule reneo:
    input:
        reneoTargets


@targetRule
rule print_stages:
    run:
        print("\nIndividual Reneo stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)


"""RULES"""
# Step 2: Obtain unitig sequences from assembly graph
include: os.path.join("rules", "gfa2fasta.smk")


# Step 3: Calculate coverage of unitig sequences
include: os.path.join("rules", "coverage.smk")


# Step 4: Scan unitig sequences for single-copy marker genes and PHROGs
include: os.path.join("rules", "genes.smk")


# Step 5: Run Reneo
include: os.path.join("rules", "reneo.smk")
