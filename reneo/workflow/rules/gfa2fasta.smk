"""
Run gfa2fasta to obtain the sequences corresponding to unitigs in the assembly graphs in FASTA format.
The assembly graph file with .GFA extension should be provided as inputs.
"""

rule run_gfa2fasta:
    input:
        GRAPH_FILE
    output:
        EDGES_FILE
    params:
        graph = GRAPH_FILE,
        output = OUTDIR,
        log = os.path.join(LOGSDIR, "gfa2fasta.log")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem_mb = config["resources"]["ram"]["mem"],
        mem = str(config["resources"]["ram"]["mem"]) + "MB",
        time = config["resources"]["ram"]["time"]
    log:
        os.path.join(LOGSDIR, "gfa2fasta.log")
    conda: 
        os.path.join("..", "envs", "reneo.yaml")
    script:
        os.path.join('..', 'scripts', 'gfa2fasta.py')