##### PURPOSE OF THIS SCRIPT
### 0) Might need to edit FASTA file, or at least validate it, to deal with grandRescue bugs
### 1) Align to the called TSSome with grandRescue, NextGenMap, or bowtie2


if config["s4U_aligner"] == "ngm":

    rule align_all:
        input:
            read="results/trimmed/{sample}.1.fastq",
            ref=ALIGN_ALL_REF,
        output:
            sam="results/align_all/{sample}.sam",
            bam="results/align_all/{sample}.bam",
        threads: 20
        conda:
            "../envs/ngm.yml"
        shell:
            """
            ngm -q {input.read} -r {input.ref} -t {threads} -o {output.sam} --slam-seq 2
            samtools view -b -@ {threads} -o {output.bam}
            """


else if config["s4U_aligner"] == "grandRescue":

    ### NOT COMPLETE; MAYBE NEVER AS grandRescue IS VERY BUGGY

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


else if config["s4U_aligner"] == "bismark"

    # rule align_all:
    #     input:

else if config["s4U_aligner"] == "bowtie2"

    rule align_all:
        input:
            sample=["results/trimmed/{sample}.1.fastq"]
            idx=muitiext(
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
            extra=config["align_all_extra"]
        threads: 20
        wrapper:
            "v7.2.0/bio/bowtie2/align"