"""
Use raw_coverage to map to calculate coverage of unitigs.
Use combine_cov to combine the coverage values of multiple samples into one file.
"""

rule koverage_tsv:
    """Generate TSV of samples and reads for Koverage"""
    output:
        os.path.join(OUTDIR,"reneo.samples.tsv")
    params:
        SAMPLE_READS
    localrule: True
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0], output[0])


rule koverage:
    """Get coverage statistics with Koverage + CoverM"""
    input:
        tsv = os.path.join(OUTDIR,"reneo.samples.tsv"),
        edges = EDGES_FILE
    params:
        out_dir = OUTDIR,
        profile = lambda wildcards: "--profile " + config["profile"] if config["profile"] else "",
    output:
        expand(os.path.join(OUTDIR, "temp", "{sample}.{ext}"),
               sample=SAMPLE_NAMES,
               ext=["bam","bam.bai"]),
        COVERAGE_FILE
    threads:
        lambda w: 1 if config["profile"] else config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time= config["resources"]["big"]["time"]
    conda:
        os.path.join("..", "envs", "koverage.yaml")
    shell:
        """
        koverage run reneo_coverage \
            --reads {input.tsv} \
            --ref {input.edges} \
            --threads {threads} \
            --output {params.out_dir} \
            {params.profile}
        """


rule find_junctions_and_pickle:
    """Parse the BAM file and create a python dictionary pickle of PE junctions for contig pairs"""
    input:
        bam = os.path.join(OUTDIR, "temp", "{sample}.bam"),
        bai = os.path.join(OUTDIR, "temp", "{sample}.bam.bai"),
    output:
        pkl = temp(os.path.join(OUTDIR, "temp", "{sample}.pkl"))
    log:
        err = os.path.join(LOGSDIR, "find_junctions_and_pickle.{sample}.log")
    benchmark:
        os.path.join(BENCH, "find_junctions_and_pickle.{sample}.log")
    conda:
        os.path.join("..", "envs", "reneo.yaml")
    script:
        os.path.join("..", "scripts", "sampleJunctions.py")


rule combine_junction_pickles:
    input:
        pkls = expand(os.path.join(OUTDIR, "temp", "{sample}.pkl"), sample=SAMPLE_NAMES)
    output:
        pkl = protected(PICKLE_FILE)
    log:
        err=os.path.join(LOGSDIR,"combine_junction_pickles.log")
    benchmark:
        os.path.join(BENCH,"combine_junction_pickles.log")
    conda:
        os.path.join("..","envs","reneo.yaml")
    script:
        os.path.join("..","scripts","combineJunctionPickles.py")
