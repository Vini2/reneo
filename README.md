# Reneo: Unraveling Viral Genomes from Metagenomes using Assembly Graphs

[Reneo](https://en.wiktionary.org/wiki/reneo) means to *unravel* or *untangle* in latin. Reneo the software was developed to unravel or untangle high-quality genomes from viral communities found within metagenomes using assembly graphs. It models cyclic viral components in the metagenomic assembly as flow networks, models as a minimum flow decomposition problem and resolves genomic paths corresponding to flow paths determined.

## Setting up Reneo

### Step 1: Downloading Reneo

```bash
git clone https://github.com/Vini2/reneo.git
```

### Step 2: Setting up Reneo environment

```bash
# Move into reneo folder
cd reneo

# Setup conda environment
conda env create -f build/environment.yml

# Activate reneo environment
conda activate reneo

# Setup reneo
pip install -e .
```

Now you can go to [Setting up Gurobi](#setting-up-gurobi) to configure Gurobi.

### Setting up Gurobi

The MFD implementation uses the linear programming solver [Gurobi](https://www.gurobi.com/). The `phables` conda environment and pip setup does not include Gurobi. You have to install Gurobi using one of the following commands depending on your package manager.

```bash
# conda
conda install -c gurobi gurobi

# pip
pip install gurobipy
```

To handle large models without any model size limitations, once you have installed Gurobi, you have to activate the (academic) license and add the key using the following command. You only have to do this once.

```bash
grbgetkey <KEY>
```

You can refer to further instructions at [https://www.gurobi.com/academia/academic-program-and-licenses/](https://www.gurobi.com/academia/academic-program-and-licenses/). 


## Quick Start Guide

### Setting up databases

```bash
reneo install
```

### Testing the setup

After setting up, run the following command to print out the Reneo help message.

```bash
reneo --help
```


## Example Usage

```bash
# Run Reneo
# locally: using 8 threads (default is 1 thread)
reneo run --input assembly_graph.gfa --reads fastq/ --threads 8
```