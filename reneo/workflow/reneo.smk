"""
Reneo: Unraveling Viral Genomes from Metagenomes.

2023, Vijini Mallawaarachchi

This is the main Snakefile to run reneo.
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")
config.update(config["reneo"])


"""PREFLIGHT CHECKS"""
include: os.path.join("rules", "databases.smk")
include: os.path.join("rules", "preflight.smk")


"""TARGETS"""
include: os.path.join("rules", "targets.smk")


"""Target rules"""
target_rules = []

def targetRule(fn):
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn

if config["profile"]:
    localrules: koverage, koverage_genomes


"""Run stages"""
@targetRule
rule all:
    input:
        # preprocessTargets,
        reneoTargets,
        postprocessTargets


@targetRule
rule preprocess:
    input:
        preprocessTargets


@targetRule
rule reneo:
    input:
        reneoTargets


@targetRule
rule postprocess:
    input:
        postprocessTargets


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


# Step 6: Postprocess genomes
include: os.path.join("rules", "postprocess.smk")