# Reneo: Unraveling Viral Genomes from Metagenomes

[Reneo](https://en.wiktionary.org/wiki/reneo) means to *unravel* or *untangle* in latin. Reneo the software was developed to unravel or untangle high-quality genomes from viral communities found within metagenomes using assembly graphs. It models cyclic viral components in the metagenomic assembly as flow networks, models as a minimum flow decomposition (MFD) problem and resolves genomic paths corresponding to flow paths determined.

## Setting up Reneo

### Downloading Reneo

Clone the Reneo GitHub repository to your machine using the following command.

```bash
git clone https://github.com/Vini2/reneo.git
```

Move into the `reneo` folder.


```bash
# Move into reneo folder
cd reneo
```

### Setting up Reneo environment

We recommend to use conda for the setup. Run the following commands to setup and activate the conda environment.

```bash
# Setup conda environment
conda env create -f build/environment.yml

# Activate reneo environment
conda activate reneo
```

Now run the following command to install Reneo to the created environment.

```bash
# Setup reneo
pip install -e .
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

```bash
reneo install
```

### Testing the setup

After setting up, run the following command to print out the Reneo help message.

```bash
reneo --help
```

### Running Reneo

```bash
# Run Reneo
# locally: using 8 threads (default is 1 thread)
reneo run --input assembly_graph.gfa --reads fastq/ --threads 8
```


##  Issues and Questions

Reneo is still under testing. Please report any issues and suggestions under [Reneo Issues](https://github.com/Vini2/reneo/issues).
