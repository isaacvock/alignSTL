##### PURPOSE OF THIS SCRIPT
### 1) Trim fastq files
### 2) Run FastQC


if config.get("PE_input", True):

    ### Trim adaptors
    rule cutadapt:
        input:
            get_input_fastqs,
        output:
            fastq1="results/trimmed/{sample}.1.fastq",
            fastq2="results/trimmed/{sample}.2.fastq",
            qc="results/trimmed/{sample}.qc.txt",
        params:
            adapters=config["adapters"],
            extra=config["cutadapt_extra"],
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8
        wrapper:
            "v7.2.0/bio/cutadapt/pe"

    ### Run fastQC
    rule fastqc:
        input:
            get_fastqc_read,
        output:
            html="results/fastqc/{sample}_r{read}_fastqc.html",
            zip="results/fastqc/{sample}_r{read}_fastqc.zip",
        log:
            "logs/fastqc/{sample}_r{read}.log",
        params:
            extra=config.get("fastqc_params"),
        resources:
            mem_mb=9000,
        threads: 4
        wrapper:
            "v2.2.1/bio/fastqc"

else:

    ### Trim adaptors
    rule cutadapt:
        input:
            get_input_fastqs,
        output:
            fastq="results/trimmed/{sample}.1.fastq",
            qc="results/trimmed/{sample}.qc.txt",
        params:
            adapters=config["adapters"],
            extra=config["cutadapt_extra"],
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8
        wrapper:
            "v7.2.0/bio/cutadapt/se"

    ### Run fastQC
    rule fastqc:
        input:
            get_fastqc_read,
        output:
            html="results/fastqc/{sample}_r1_fastqc.html",
            zip="results/fastqc/{sample}_r1_fastqc.zip",
        log:
            "logs/fastqc/{sample}.log",
        params:
            extra=config.get("fastqc_params"),
        resources:
            mem_mb=9000,
        threads: 4
        wrapper:
            "v2.2.1/bio/fastqc"
