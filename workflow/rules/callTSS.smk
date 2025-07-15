##### PURPOSE OF THIS SCRIPT
### 1) Align ctl data with bowtie2
### 2) Call TSS with Adelman lab script
### 3) Create TSS FASTA and GTF


### Align samples used for TSScall


rule align_ctl:
    input:
        sample=["results/trimmed/{ctl}.1.fastq"],
        idx=multiext(
            config.get("bowtie2_index"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "results/align_ctl/{ctl}.bam",
    log:
        "logs/align_ctl/{ctl}.log",
    params:
        extra=config.get("align_ctl_extra"),
        index=lambda wildcards, input: input[1].rsplit(".1.bt2", 1)[0],
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


### Merge and filter bam file for TSScall


rule merge_ctl:
    input:
        expand("results/align_ctl/{CTL}.bam", CTL=CTL_SAMPLES),
    output:
        "results/merge_ctl/merged.bam",
    threads: 8
    log:
        "logs/merge_ctl/merge_ctl.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} &> {log}
        """


rule filter_and_sort_ctl:
    input:
        "results/merge_ctl/merged.bam",
    output:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    log:
        "logs/filter_and_sort_ctl/filter_and_sort_ctl.log",
    threads: 8
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -@ {threads} -b -h -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x800 {input} 2>> {log} | \
            samtools sort -@ {threads} - 2>> {log} > {output}
        """


### Make bedgraph files for TSScall
rule make_forward_bedgraph_ctl:
    input:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    output:
        "results/make_bedgraph_ctl/forward.bedgraph",
    threads: 1
    conda:
        "../envs/bedtools.yml"
    params:
        shellscript=workflow.source_path("../scripts/bedtools.sh"),
    log:
        "logs/make_forward_bedgraph_ctl/make_forward_bedgraph_ctl.log",
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {input} {output} &> {log}
        """


rule make_reverse_bedgraph_ctl:
    input:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    output:
        "results/make_bedgraph_ctl/reverse.bedgraph",
    threads: 1
    conda:
        "../envs/bedtools.yml"
    params:
        shellscript=workflow.source_path("../scripts/bedtools.sh"),
    log:
        "logs/make_reverse_bedgraph_ctl/make_reverse_bedgraph_ctl.log",
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {input} {output} &> {log}
        """


### Get chromosome size information from bam file
rule get_chr_sizes:
    input:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    output:
        "results/get_chr_sizes/genome.chrom.sizes",
    log:
        "logs/get_chr_sizes/get_chr_sizes.log",
    conda:
        "../envs/chrom.yml"
    params:
        shellscript=workflow.source_path("../scripts/chrom.sh"),
    threads: 1
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {input} {output} &> {log}
        """


### Use TSScall script from Adelman lab to ID TSS
rule callTSS:
    input:
        fwd="results/make_bedgraph_ctl/forward.bedgraph",
        rev="results/make_bedgraph_ctl/reverse.bedgraph",
        chrom="results/get_chr_sizes/genome.chrom.sizes",
        anno=config.get("genome_gtf"),
    output:
        bed="results/callTSS/callTSS.tss.bed",
        details="results/callTSS/callTSS.tss_detail.txt",
        clusters="results/callTSS/callTSS.tss_clusters.bed",
    params:
        pyscript=workflow.source_path("../scripts/TSScall.py"),
        fdr=config.get("TSScall_fdr"),
        utss_filter_size=config.get("TSScall_utss_filter_size"),
        utss_search_window=config.get("TSScall_utss_search_window"),
        bidirectional_thresh=config.get("TSScall_bidirectional_thresh"),
        cluster_thresh=config.get("TSScall_cluster_thresh"),
        call_method=config.get("TSScall_call_method"),
        bin_winner_size=config.get("TSScall_bin_winner_size"),
        a_search_win=config.get("annotation_search_window", 1000),
        a_join_dist=config.get("annotation_join_distance", 200),
    log:
        "logs/callTSS/callTSS.log",
    conda:
        "../envs/TSScall.yml"
    threads: 1
    shell:
        """
        chmod +x {params.pyscript}
        python {params.pyscript} \
            --fdr {params.fdr} \
            --utss_filter_size {params.utss_filter_size} \
            --utss_search_window {params.utss_search_window} \
            --bidirectional_threshold {params.bidirectional_thresh} \
            --cluster_threshold {params.cluster_thresh} \
            --cluster_bed {output.clusters} \
            --detail_file {output.details} \
            --annotation_file {input.anno} \
            --call_method {params.call_method} \
            --bin_winner_size {params.bin_winner_size} \
            --annotation_search_window {params.a_search_win} \
            --annotation_join_distance {params.a_join_dist} \
            {input.fwd} {input.rev} {input.chrom} {output.bed} \
            &> {log}
        """


### Use my custom script to get TSS FASTA
rule make_TSSome:
    input:
        bed="results/callTSS/callTSS.tss.bed",
        fasta=config.get("genome_fasta"),
        details="results/callTSS/callTSS.tss_detail.txt",
    output:
        fasta="results/make_TSSome/TSSome.fasta",
    params:
        Rscript=workflow.source_path("../scripts/createTSSome.R"),
        extra=config.get("createTSSome_extra"),
    log:
        "logs/make_TSSome/makeTSSome.log",
    conda:
        "../envs/TSSome.yml"
    threads: 1
    shell:
        r"""
        chmod +x {params.Rscript}
        {params.Rscript} \
            --fasta {input.fasta} \
            --bed {input.bed} \
            --output_fasta {output.fasta} \
            --details {input.details} \
            {params.extra} &> {log}
        """
