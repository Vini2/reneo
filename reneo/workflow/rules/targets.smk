
preprocessTargets = []
reneoTargets = []
postprocessTargets = []


"""MISC"""
COVERAGE_FILE = os.path.join(RESDIR, 'reneo.coverage.tsv')
VOG_ANNOT = os.path.join(OUTDIR, 'all.hmmVOG.tbl')
SMG_FILE = os.path.join(OUTDIR, 'edges.fasta.hmmout')
GRAPH_FILE = INPUT
PICKLE_FILE = os.path.join(OUTDIR, "PE_junctions.pkl")


"""PREPROCESSING TARGETS"""
EDGES_FILE = os.path.join(OUTDIR, "edges.fasta")
preprocessTargets.append(EDGES_FILE)

BAM_PATH = os.path.join(OUTDIR, 'temp')
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam"), sample=SAMPLE_NAMES))
preprocessTargets.append(expand(os.path.join(BAM_PATH, "{sample}.bam.bai"), sample=SAMPLE_NAMES))

COVERAGE_PATH = os.path.join(OUTDIR, 'coverage_rpkm/')
preprocessTargets.append(COVERAGE_FILE)
if config["hmmsearch"]:
    preprocessTargets.append(SMG_FILE)
    preprocessTargets.append(VOG_ANNOT)


"""PHABLES TARGETS"""
RESOLVED_GENOMES = os.path.join(RESDIR, "resolved_paths.fasta")
reneoTargets.append(RESOLVED_GENOMES)

RESOLVED_GENOME_INFO = os.path.join(RESDIR, "resolved_genome_info.txt")
reneoTargets.append(RESOLVED_GENOME_INFO)

RESOLVED_COMP_INFO = os.path.join(RESDIR, "resolved_component_info.txt")
reneoTargets.append(RESOLVED_COMP_INFO)

COMP_VOGS = os.path.join(RESDIR, "component_vogs.txt")
reneoTargets.append(COMP_VOGS)


"""POSTPROCESSING TARGETS"""
GENOME_KOVERAGE_RES = os.path.join(RESDIR, "sample_coverage.tsv")
postprocessTargets.append(GENOME_KOVERAGE_RES)

GENOME_READ_COUNTS = os.path.join(RESDIR, "sample_genome_read_counts.tsv")
postprocessTargets.append(GENOME_READ_COUNTS)