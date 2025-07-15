# alignSTL: aligning STL-seq data

This Snakemake pipeline is designed to align STL-seq data. It is structured similarly to [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR) and can be run in the same way as that pipeline is run:

``` bash
### 
# PREREQUISITES: INSTALL MAMBA or CONDA AND GIT (only need to do once per system)
###

# CREATE ENVIRONMENT (only need to do once per system)
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do once per dataset)
mkdir path/to/working/directory
cd path/to/working/directory

# DEPLOY PIPELINE TO YOUR WORKING DIRECTORY (only need to do once per dataset)
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/alignSTL.git . --branch main

###
# EDIT CONFIG FILE (need to do once for each new dataset)
###

# RUN PIPELINE

# See [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details on all of the configurable parameters
snakemake --cores all --use-conda --rerun-triggers mtime --keep-going
```

## Details

Two general alignment strategies are implemented (the align_target parameter in the config sets which of these gets used):

1. Alignment to a provided FASTA file (e.g., a genome)
2. Alignment to a provided genome FASTA file followed by calling TSS's, creating a TSSome from the called TSS's, and then aligning to the called TSSome

TSS calling is done with [TSScall](https://github.com/lavenderca/TSScall), and is performed in both cases so that the output of this pipeline can be provied as input to fastq2EZbakR (if performing alignment to a provided genome FASTA file, you will need the TSS GTF created by alignSTL to run fastq2EZbakR). 

In either case, there are currently 3 aligners implemented in alignSTL (the s4U_aligner parameter in the config sets which of these gets used):

1. [NextGenMap](https://github.com/Cibiv/NextGenMap): implements a unique scoring function that turns off T-to-C mismatch penalization, perfect for STL-seq data.
2. [Bismark](https://github.com/FelixKrueger/Bismark) + bowtie2: alignment to a 3-base genome, an extreme form of T-to-C mismatch penalization removal
3. [Bowtie2](https://github.com/BenLangmead/bowtie2): Standard 4-base genome alignment


