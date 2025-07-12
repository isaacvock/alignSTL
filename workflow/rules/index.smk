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
    threads: 20
    wrapper:
        "v7.2.0/bio/bowtie2/build"


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
        extra=config.get("bowtie2_build_TSSome_extra"),
    threads: 20
    wrapper:
        "v7.2.0/bio/bowtie2/build"