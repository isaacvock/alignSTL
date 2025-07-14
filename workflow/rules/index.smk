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
        ALIGN_ALL_REF,  
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
        TC=multiext("bismark_genome/Bisulfite_Genome/CT_conversion/BS_CT"
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2"),
        AG=multiext("bismark_genome/Bisulfite_Genome/GA_conversion/BS_GS"
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2"),
    params:
        extra=config.get("bismark_prep_genome_extra", ""),
    log:
        "logs/bismark_prepare_genome/bismark.log"
    shell:
        """
        bismark_genome_preparation --slam {params.extra} bismark_genome/ &> {log}
        """
