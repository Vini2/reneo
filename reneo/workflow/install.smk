"""
Reneo: Unraveling Viral Genomes from Metagenomes using Assembly Graphs.

2023, Vijini Mallawaarachchi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")

include: "rules/00_database_preflight.smk"


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
        vog_db_path = directory(os.path.join(DBPATH, "VOG"))
    output:
        os.path.join(DBPATH, config['vog_db_file'])
    conda:
        os.path.join("envs", "curl.yaml")
    shell:
        """
            curl -Lo {params.file} {params.url}
            mkdir {params.vog_db_path}
            tar -xf {params.file} -C {params.vog_db_path}
            cat {params.vog_db_path}/* > {output}
            rm -rf {params.file}
            rm -rf {params.vog_db_path}
        """