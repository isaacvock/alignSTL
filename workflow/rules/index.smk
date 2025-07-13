##### PURPOSE OF THIS SCRIPT
### 1) Build bowtie2 alignment index
### 2) Build NextGenMap alignment index
### 3) Build grandRescue/gedi indices


### bowtie2 build index

rule bowtie2_build_genome:
    input:
        ref=config.get("genome_fasta"),
    output:
        multiext(
            config.get("bowtie2_index"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra=config.get("bowtie2_build_genome_extra"),
        index=config.get("bowtie2_index"),
    threads: 20
    conda:
        "../envs/bowtie2.yml"
    shell:
        """
        bowtie2-build --threads {threads} {params.extra} {input.ref} {params.index} &> {log}
        """


rule bowtie2_build_TSSome:
    input:
        ref="results/make_TSSome/TSSome.fasta",
    output:
        multiext(
            config.get("bowtie2_TSS_index"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra=config.get("bowtie2_build_genome_extra"),
        index=config.get("bowtie2_TSS_index"),
    threads: 20
    conda:
        "../envs/bowtie2.yml"
    shell:
        """
        bowtie2-build --threads {threads} {params.extra} {input.ref} {params.index} &> {log}
        """


### bismark

rule copy_genome:
    input:
        config.get("genome_fasta"),
    output:
        "bismark_genome/genome.fasta",
    log:
        "logs/copy_genome/copy.log"
    shell:
        "cp {input} {output} &> {log}"

rule bismark_prepare_genome:
    input:
        genome="bismark_genome/genome.fasta"
    output:
        ""
    shell:
        """
        bismark_genome_preparation --path_to_aligner ./ bismark_genome/
        """
