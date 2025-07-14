##### PURPOSE OF THIS SCRIPT
### 0) Might need to edit FASTA file, or at least validate it, to deal with grandRescue bugs
### 1) Align to the called TSSome with grandRescue, NextGenMap, or bowtie2


if config["s4U_aligner"] == "ngm":

    rule align_all:
        input:
            read="results/trimmed/{sample}.1.fastq",
            ref=ALIGN_ALL_REF,
        output:
            "results/align_all/{sample}.bam",
        threads: 20
        log:
            "logs/align_all/{sample}.log"
        conda:
            "../envs/ngm.yml"
        shell:
            """
            set -euo pipefail
            ngm -q {input.read} -r {input.ref} -t {threads} --slam-seq 2 2>> {log} | samtools view -Sb - 2>> {log} > {output}
            """


elif config["s4U_aligner"] == "grandRescue":

    ### NOT COMPLETE; MAYBE NEVER AS grandRescue IS VERY BUGGY

    dummy = 2 + 2
    # rule prep_pseduotranscriptome:
    #     input:
    #         fasta=ALIGN_ALL_REF,
    #         gtf=ALIGN_ALL_GTF,
    #     output:
    #         ### TO-DO: figure out what files get created
    #         multiext(
    #             "gedi_genome/pseudoTranscriptome/",
    #             ".gtf",
    #             ".index",
    #             "_T2C.fa" 
    #         )
    #     params:
    #         gedi_path=config["gedi_path"],
    #         extra=config["CreatePseudo_extra"],
    #     threads: 1
    #     shell:
    #         """
    #         {params.gedi_path} -e CreatePseudo -fasta {input.fasta} -gtf {input.gtf} -o gedi_genome/pseudoTranscriptome/
    #         """

    # rule prep_gediGenome:
    #     input:
    #         fasta=ALIGN_ALL_REF,
    #         gtf=ALIGN_ALL_GTF,
    #     output:
    #         ### TO-DO: Figure out what gets created
    #     params:
    #         name=config["gedi_genome_name"],
    #     shell:
    #         """
    #         gedi -e IndexGenome -s {input.fasta} -a {input.gtf} -n {params.name}
    #         """

    # rule prep_gediPseudoGenome:
    #     input:
    #         fasta=ALIGN_ALL_PSEUDOREF,
    #         gtf=ALIGN_ALL_PSEUDOGTF,
    #     output:
    #         ### TO-DO: Figure out what gets created
    #     params:
    #         name=config["gedi_genome_pseudoname"],
    #     shell:
    #         """
    #         gedi -e IndexGenome -s {input.fasta} -a {input.gtf} -n {params.name}
    #         """


elif config["s4U_aligner"] == "bismark":

    rule align_all:
        input:
            sample=["results/trimmed/{sample}.1.fastq"],
        output:
            "results/align_all/{sample}.bam",
        log:
            "logs/align_all/{sample}.log",
        params:
            extra=config.get("bismark_align_extra", ""),
        shell:
            """
            bismark -p {threads} --slam {params.extra} --genome bismark_genome/ {input.sample} -o /results/align_all/ --basename {wildcards.sample}
            """

elif config["s4U_aligner"] == "bowtie2":

    rule align_all:
        input:
            sample=["results/trimmed/{sample}.1.fastq"],
            idx=multiext(
                ALIGN_ALL_INDEX,
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        output:
            "results/align_all/{sample}.bam",
        log:
            "logs/align_all/{sample}.log",
        params:
            extra=config.get("align_ctl_extra"),
            index=config.get("bowtie2_index")
        threads: 20
        conda:
            "../envs/bowtie2.yml"
        shell:
            """
            set -euo pipefail
            bowtie2 --threads {threads} \
                -U {input.sample} \
                -x {params.index} {params.extra} 2>> {log} \
                | samtools view -@ {threads} -h -b -o {output} - 2>> {log}
            """
