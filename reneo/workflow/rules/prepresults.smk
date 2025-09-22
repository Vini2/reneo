"""
Run prepresult to format the initial binning result from an existing binning tool.
The folder with the bins should be provided as inputs.
"""

rule run_prepresults:
    input:
        BINS
    output:
        BINS_FILE
    params:
        bins = BINS,
        output = BINS_FILE
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem_mb = config["resources"]["ram"]["mem"],
        mem = str(config["resources"]["ram"]["mem"]) + "MB",
        time = config["resources"]["ram"]["time"]
    log:
        stderr = os.path.join(LOGSDIR, "prep_results.err"),
        stdout = os.path.join(LOGSDIR, "prep_results.out"),
    conda:
        os.path.join("..", "envs", "reneo.yaml")
    script:
        os.path.join('..', 'scripts', 'prep_results.py')
