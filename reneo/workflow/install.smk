"""
Reneo: Unraveling Viral Genomes from Metagenomes.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")

include: os.path.join("rules", "databases.smk")


"""TARGETS"""
db_files = []

db_files.append(os.path.join(DBPATH, config['smg_hmm_file']))
db_files.append(os.path.join(DBPATH, config['vog_db_file']))


"""RUN SNAKEMAKE"""
rule all:
    input:
        db_files


"""RULES"""
rule smg_hmm_download:
    params:
        url=os.path.join(config['smg_hmm'])
    output:
        os.path.join(DBPATH, config['smg_hmm_file'])
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lo {output} {params.url}
        """

rule vog_db_download:
    params:
        url=os.path.join(config['vog_db']),
        file=os.path.join(DBPATH, config['vog_db_tar']),
    output:
        os.path.join(DBPATH, config['vog_db_file'])
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lo {params.file} {params.url}
            tar -xOf {params.file} > {output}
            rm -rf {params.file}
        """