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


def get_input_fastqs(wildcards):

    fastq_path = config["samples"].get(wildcards.sample, None)
    if fastq_path is None:
        raise ValueError(f"No path found for sample {wildcards.sample}")

    print(f"path is: {fastq_path}")
    current_directory = os.getcwd()
    print("Current Working Directory:", current_directory)
    all_files = os.listdir(fastq_path)
    print(f"All files in directory: {all_files}")

    fastq_files = []

    # List all entries in the directory and store full paths in the list
    for entry in os.listdir(fastq_path):
        full_path = os.path.join(fastq_path, entry)
        fastq_files.append(full_path)

    # Filter files to only include those that end with .fastq or .fastq.gz
    # fastq_files = [f"data/fastq/{wildcards.sample}/{wildcards.sample}.{n}.fastq.gz" for n in ("1", "2")] #sorted([os.path.join(directory, f) for f in all_files if f.endswith('.fastq') or f.endswith('.fastq.gz')])
    # fastq_files = sorted(glob.glob(f"{fastq_path}/*.fastq*"))
    print(f"Found files: {sorted(fastq_files)}")
    return sorted(fastq_files)
