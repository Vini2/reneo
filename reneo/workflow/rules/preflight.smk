"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

from metasnek import fastq_finder


"""
Setting the directory variables
"""
OUTDIR = config['output']
TMPDIR = os.path.join(OUTDIR, "temp")
RESDIR = os.path.join(config['output'], "results")
LOGSDIR = os.path.join(OUTDIR, 'logs')
BENCH = os.path.join(OUTDIR, 'benchmarks')


"""
Parse the samples
"""
SAMPLE_READS = dict(sorted(fastq_finder.parse_samples_to_dictionary(config['reads']).items()))
SAMPLE_NAMES = list(SAMPLE_READS.keys())


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(LOGSDIR):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(LOGSDIR))
        for logfile in oldLogs:
            os.unlink(os.path.join(LOGSDIR, logfile))


onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nReneo ran successfully!\n\n')


onerror:
    """Print an error message"""
    sys.stderr.write('\n\nReneo run failed\n\n')
