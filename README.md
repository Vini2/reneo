<p align="center">
  <img src="https://raw.githubusercontent.com/Vini2/reneo/develop/reneo_logo.png" width="200" title="reneo logo" alt="reneo logo">
</p>

Reneo: Unraveling Viral Genomes from Metagenomes
===============

[![DOI](https://zenodo.org/badge/619432085.svg)](https://zenodo.org/badge/latestdoi/619432085)
![GitHub](https://img.shields.io/github/license/vini2/reneo)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/reneo/README.html)
[![Conda](https://img.shields.io/conda/v/bioconda/reneo)](https://anaconda.org/bioconda/reneo)
[![Conda](https://img.shields.io/conda/dn/bioconda/reneo)](https://anaconda.org/bioconda/reneo)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub last commit (by committer)](https://img.shields.io/github/last-commit/Vini2/reneo?color=8a35da)
[![CI](https://github.com/Vini2/reneo/actions/workflows/test_reneo.yml/badge.svg)](https://github.com/Vini2/reneo/actions/workflows/test_reneo.yml)
[![CodeQL](https://github.com/Vini2/reneo/actions/workflows/codeql.yml/badge.svg)](https://github.com/Vini2/reneo/actions/workflows/codeql.yml)


[Reneo](https://en.wiktionary.org/wiki/reneo) means to *unravel* or *untangle* in Latin. Reneo is a software developed to unravel or untangle high-quality genomes from viral communities (including both prokaryotic and eukaryotic viruses) found within metagenomes using assembly graphs. Reneo identifies viral components in the metagenomic assembly using virus orthologous groups from [VOGDB](https://vogdb.org/), models as flow networks and solves a minimum flow decomposition (MFD) problem to resolve genomic paths. Reneo was motivated based on a bacteriophage recovery tool named [Phables](https://github.com/Vini2/phables), specifically to extend the capabilities of Phables to all viruses.

**NEW:** Reneo is available on bioconda at [https://anaconda.org/bioconda/reneo](https://anaconda.org/bioconda/reneo)

## Setting up Reneo

### Installing Reneo using conda (recommended)

You can install Reneo from bioconda at [https://anaconda.org/bioconda/reneo](https://anaconda.org/bioconda/reneo). Make sure you have [`conda`](https://docs.conda.io/en/latest/) installed.

```bash
# create conda environment and install reneo
conda create -n reneo -c conda-forge -c anaconda -c bioconda reneo

# activate environment
conda activate reneo
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

### Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `reneo` conda environment does not include Gurobi. You have to install Gurobi using the following command.

```bash
conda install -c gurobi gurobi
```

To handle large models without any model size limitations, once you have installed Gurobi, you have to activate the (academic) license and add the key using the following command. You only have to do this once.

```bash
grbgetkey <KEY>
```

You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/). 


## Quick Start Guide

### Setting up databases

Run the following command to download and set up the databases used in Reneo.

```bash
reneo install
```

### Testing the setup

After setting up, run the following command to print out the Reneo help message.

```bash
reneo --help
```

You can simulate a Reneo run using the following command.

```bash
reneo simulate
```

You can also run Reneo with the test dataset provided.

```bash
reneo test
```

### Running Reneo

```bash
# Run Reneo
# locally: using 16 threads (default is 8 threads)
reneo run --input assembly_graph.gfa --reads fastq/ --threads 16
```


##  Issues and Questions

Reneo is still under testing. Please report any issues and suggestions under [Reneo Issues](https://github.com/Vini2/reneo/issues).


## Acknowledgement

Reneo uses the [Gurobi](https://www.gurobi.com/) implementation of [MFD-ILP](https://github.com/algbio/MFD-ILP) and code snippets from [Phables](https://github.com/Vini2/phables/). The Reneo logo was designed by [Laura Inglis](https://fame.flinders.edu.au/people/2021/01/01/laura-inglis).


## Citation

Reneo is a work in progress and the manuscript is currently in preparation. In the meantime, please cite Reneo as

```
V Mallawaarachchi, MJ Roach, LK Inglis and RA Edwards (2023). Reneo: Unraveling Viral Genomes from Metagenomes. Available at https://github.com/Vini2/reneo DOI: 10.5281/zenodo.8263066
```
