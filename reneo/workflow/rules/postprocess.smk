rule combine_genomes_and_unresolved_edges:
    """Combine resolved genomes and unresolved edges"""
    input:
        genomes = RESOLVED_GENOMES,
        unresolved_edges = os.path.join(RESDIR, "unresolved_virus_like_edges.fasta")
    output:
        os.path.join(RESDIR, "genomes_and_unresolved_edges.fasta")
    shell:
        """
        cat {input.genomes} {input.unresolved_edges} > {output}
        """


rule koverage_genomes:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(OUTDIR,"reneo.samples.tsv"),
        sequences = os.path.join(RESDIR, "genomes_and_unresolved_edges.fasta")
    params:
        out_dir = OUTDIR,
        profile = lambda wildcards: "--profile " + config["profile"] if config["profile"] else "",
    output:
        os.path.join(RESDIR, "sample_coverage.tsv")
    threads:
        lambda w: 1 if config["profile"] else config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time = config["resources"]["big"]["time"]
    conda:
        os.path.join("..", "envs", "koverage.yaml")
    shell:
        """
        koverage run \
            --reads {input.tsv} \
            --ref {input.sequences} \
            --threads {threads} \
            --output {params.out_dir} \
            {params.profile}
        """
