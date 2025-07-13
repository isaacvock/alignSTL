# import basic packages
import pandas as pd
from snakemake.utils import validate
import glob
import os


### GENERAL HELPER FUNCTIONS/VARIABLES USED IN ALL CASES

# Sample names to help expanding lists of all bam files
# and to aid in defining wildcards
SAMP_NAMES = list(config.get("samples").keys())

CTL_SAMPLES = list(config.get("ctl_samples"))

# Which alignment index/reference genome to use?
if config.get("align_target", "TSSome") == "TSSome":

    ALIGN_ALL_INDEX = config.get("bowtie2_TSS_index")
    ALIGN_ALL_REF = "results/make_TSSome/TSSome.fasta"

else:

    ALIGN_ALL_INDEX = config.get("bowtie2_index")
    ALIGN_ALL_REF = config.get("genome_fasta")


### What is desired final output?

def get_target():
    target = []

    target.append("results/multiqc/multiqc_report.html")

    target.append(
        expand(
            "results/align_all/{SAMPLE}.bam",
            SAMPLE = SAMP_NAMES
        )
    )


### FastQC input info

def get_fastqc_read(wildcards):
    
    return expand(
        "results/trimmed/{SID}.{READ}.fastq",
        SID=wildcards.sample,
        READ=wildcards.read,
    )

