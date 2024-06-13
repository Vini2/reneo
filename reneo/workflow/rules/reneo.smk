rule run_reneo:
    input:
        graph = GRAPH_FILE,
        coverage = COVERAGE_FILE,
        pickle = PICKLE_FILE,
        bams = expand(os.path.join(BAM_PATH, "{sample}.{ext}"), sample=SAMPLE_NAMES, ext=["bam","bam.bai"]),
        other = preprocessTargets
    output:
        genomes_fasta = os.path.join(RESDIR, "resolved_paths.fasta"),
        genome_info = os.path.join(RESDIR, "resolved_genome_info.txt"),
        component_info = os.path.join(RESDIR, "resolved_component_info.txt"),
        vog_comp_info = os.path.join(RESDIR, "component_vogs.txt"),
        unresolved_edges= os.path.join(RESDIR,"unresolved_virus_like_edges.fasta"),
    params:
        genomes_folder = lambda w: directory(os.path.join(RESDIR,"resolved_viruses")) if config["split_paths"] else None,
        unitigs= lambda w: os.path.join(RESDIR,"resolved_edges.fasta") if config["unitigs"] else None,
        vogs=lambda w: VOG_ANNOT if config["hmmsearch"] else None,
        hmmout=lambda w: SMG_FILE if config["hmmsearch"] else None,
        bampath = BAM_PATH,
        minlength = config['minlength'],
        mincov = config['mincov'],
        compcount = config['compcount'],
        maxpaths = config['maxpaths'],
        mgfrac = config['mgfrac'],
        evalue = config['evalue'],
        hmmscore = config['hmmscore'],
        nvogs = config['nvogs'],
        covtol = config['covtol'],
        alpha = config['alpha'],
        output = RESDIR,
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time = config["resources"]["big"]["time"]
    log:
        stderr = os.path.join(LOGSDIR, "reneo_output.err"),
        stdout = os.path.join(LOGSDIR, "reneo_output.out")
    conda:
        os.path.join("..", "envs", "reneo.yaml")
    script:
        os.path.join("..", "scripts", "reneo.py")
