"""
Use FragGeneScan and HMMER to scan for bacterial single-copy marker genes in unitigs.
Use Prodigal and hmmsearch to scan for VOGs in unitigs.
"""

rule scan_smg:
    input:
        genome = EDGES_FILE,
        hmm = os.path.join(DBPATH, "marker.hmm"),
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time = config["resources"]["big"]["time"]
    output:
        hmmout = os.path.join(TMPDIR, "edges.fasta.hmmout")
    params:
        frag = EDGES_FILE + ".frag",
        frag_faa = EDGES_FILE + ".frag.faa",
    log:
        frag_out=os.path.join(LOGSDIR, "smg_scan_frag_out.log"),
        frag_err=os.path.join(LOGSDIR, "smg_scan_frag_err.log"),
        hmm_out=os.path.join(LOGSDIR, "smg_scan_hmm_out.log"),
        hmm_err=os.path.join(LOGSDIR, "smg_scan_hmm_err.log")
    conda: 
        os.path.join("..", "envs", "genes.yaml")
    shell:
        """
            run_FragGeneScan.pl -genome={input.genome} -out={params.frag} -complete=0 -train=complete -thread={threads} 1>{log.frag_out} 2>{log.frag_err}
            hmmsearch --domtblout {output.hmmout} --cut_tc --cpu {threads} {input.hmm} {params.frag_faa} 1>{log.hmm_out} 2> {log.hmm_err}
        """


rule scan_vogs:
    input:
        genomes = EDGES_FILE,
        db = os.path.join(DBPATH,"AllVOG.hmm")
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem_mb = config["resources"]["big"]["mem"],
        mem = str(config["resources"]["big"]["mem"]) + "MB",
        time = config["resources"]["big"]["time"]
    output:
        os.path.join(TMPDIR, "all.hmmVOG.tbl")
    params:
        genes = os.path.join(TMPDIR, "genes.fna"),
        proteins = os.path.join(TMPDIR, "proteins.faa"),
        hmm_out = os.path.join(TMPDIR, "vog_hmm_output.txt"),
    log:
        os.path.join(LOGSDIR, "vogs_scan.log")
    conda: 
        os.path.join("..", "envs", "genes.yaml")
    shell:
        """
        prodigal -i {input.genomes} -d {params.genes} -a {params.proteins} -p meta -g 11 &> {log}
        hmmsearch --cpu {threads} -E 1.0e-05 -o {params.hmm_out} --tblout {output} {input.db} {params.proteins}
        """