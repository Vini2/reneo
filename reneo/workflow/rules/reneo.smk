rule run_reneo:
    input:
        graph = GRAPH_FILE,
        coverage = COVERAGE_FILE,
        bams = expand(os.path.join(BAM_PATH, "{sample}.{ext}"), sample=SAMPLE_NAMES, ext=["bam","bam.bai"])
        # VOG_ANNOT,
        # SMG_FILE,
        # preprocessTargets
    output:
        genomes_fasta = os.path.join(OUTDIR, "resolved_paths.fasta"),
        genomes_folder = directory(os.path.join(OUTDIR, "resolved_viruses")),
        genome_info = os.path.join(OUTDIR, "resolved_genome_info.txt"),
        unitigs = os.path.join(OUTDIR, "resolved_edges.fasta"),
        component_info = os.path.join(OUTDIR, "resolved_component_info.txt"),
        vog_comp_info = os.path.join(OUTDIR, "component_vogs.txt"),
        unresolved_edges = os.path.join(OUTDIR, "unresolved_virus_like_edges.fasta"),
    params:
        # graph = GRAPH_FILE,
        hmmout = None,
        vogs = None,
        # coverage = COVERAGE_FILE,
        bampath = BAM_PATH,
        minlength = ML,
        mincov = MC,
        compcount = CC,
        maxpaths = MP,
        mgfrac = MGF,
        evalue = EV,
        hmmscore = HS,
        covtol = CT,
        alpha = AL,
        output = OUTDIR,
        # nthreads = config["resources"]["big"]["cpu"],
        log = os.path.join(LOGSDIR, "reneo_output.log")
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time = config["resources"]["big"]["time"]
    log:
        os.path.join(LOGSDIR, "reneo_output.log")
    conda:
        os.path.join("..", "envs", "reneo.yaml")
    script:
        os.path.join("..", "scripts", "reneo.py")
