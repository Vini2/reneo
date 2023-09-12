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
        os.path.join(OUTDIR, "results", "sample_coverm_coverage.tsv")
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time= config["resources"]["big"]["time"]
    conda:
        os.path.join("..", "envs", "koverage.yaml")
    shell:
        """
        koverage run coverm \
            --reads {input.tsv} \
            --ref {input.edges} \
            --threads {threads} \
            --output {params.out_dir} \
            {params.profile}
        """


rule run_combine_cov:
    """Sample\tContig\tCount\tRPKM\tTPM\tMean\tCovered_bases\tVariance\n"""
    input:
        os.path.join(OUTDIR, "results", "sample_coverm_coverage.tsv")
    output:
        os.path.join(OUTDIR, "coverage.tsv")
    shell:
        """
        sed -i '1d' {input}
        awk -F '\t' '{{ sum[$2] += $6 }} END {{ for (key in sum) print key, sum[key] }}' {input} > {output}
        """
