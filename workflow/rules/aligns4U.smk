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
            "logs/align_all/{sample}.log",
        params:
            extra=config.get("NGM_extra", ""),
        conda:
            "../envs/ngm.yml"
        shell:
            """
            set -euo pipefail
            ngm -q {input.read} -r {input.ref} -t {threads} --slam-seq 2 2>> {log} | samtools view -Sb - 2>> {log} > {output}
            """

elif config["s4U_aligner"] == "bismark":

    rule align_all:
        input:
            sample=["results/trimmed/{sample}.1.fastq"],
            TC=multiext(
                "bismark_genome/Bisulfite_Genome/CT_conversion/BS_CT" ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
            AG=multiext(
                "bismark_genome/Bisulfite_Genome/GA_conversion/BS_GS" ".1.bt2",
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
            extra=config.get("bismark_align_extra", ""),
        conda:
            "../envs/bismark.yml"
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
            extra=config.get("bowtie2_align_extra", ""),
            index=ALIGN_ALL_INDEX,
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
