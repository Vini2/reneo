
preprocessTargets = []
reneoTargets = []


"""PREPROCESSING TARGETS"""
EDGES_FILE = os.path.join(OUTDIR, "edges.fasta")
preprocessTargets.append(EDGES_FILE)

BAM_PATH = os.path.join(OUTDIR, 'bam_files/')
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam"), sample=SAMPLES))
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam.bai"), sample=SAMPLES))

COVERAGE_PATH = os.path.join(OUTDIR, 'coverage_rpkm/')
preprocessTargets.append(expand(os.path.join(COVERAGE_PATH, "{sample}_rpkm.tsv"), sample=SAMPLES))
preprocessTargets.append(os.path.join(OUTDIR, "coverage.tsv"))
preprocessTargets.append(os.path.join(OUTDIR, "edges.fasta.hmmout"))

preprocessTargets.append(os.path.join(OUTDIR, "all.hmmVOG.tbl"))


"""MISC"""
COVERAGE_FILE = os.path.join(OUTDIR, 'coverage.tsv')
VOG_ANNOT = os.path.join(OUTDIR, 'all.hmmVOG.tbl')
SMG_FILE = os.path.join(OUTDIR, 'edges.fasta.hmmout')
GRAPH_FILE = INPUT


"""RENEO TARGETS"""
reneoTargets.append(os.path.join(OUTDIR, "resolved_genome_info.txt"))
reneoTargets.append(os.path.join(OUTDIR, "resolved_component_info.txt"))
reneoTargets.append(os.path.join(OUTDIR, "component_vogs.txt"))
