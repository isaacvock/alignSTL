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
        "results/align_ctl/{ctl}.bam"
    log:
        "logs/align_ctl/{ctl}.log"
    params:
        extra=config.get("align_ctl_extra")
    threads: 20
    wrapper:
        "v7.2.0/bio/bowtie2/align"


### Merge and filter bam file for TSScall

rule merge_ctl:
    input:
        expand("results/align_ctl/{CTL}.bam", CTL = CTL_SAMPLES),
    output:
        "results/merge_ctl/merged.bam"
    threads: 8
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools merge -@ {threads} -o {output} {input}
        """


rule filter_and_sort_ctl:
    input:
        "results/merge_ctl/merged.bam"
    output:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    log:
        "logs/filter_and_sort_ctl/filter_and_sort_ctl.log",
    threads: 8
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -@ {threads} -b -h -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x800 {input} | \
            samtools sort -@ {threads} - > {output}
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
    shell:
        """
        # Make bedgraph
        bedtools genomecov -ibam {input} -bg -5 -strand + | awk '{ printf "%s \t %d \t %d \t %d\n", $1,$2,$3,$4 }' > {output}
        """

rule make_reverse_bedgraph_ctl:
    input:
        "results/filter_and_sort_ctl/filter_and_sort.bam",
    output:
        "results/make_bedgraph_ctl/reverse.bedgraph"
    threads: 1
    conda:
        "../envs/bedtools.yml"  
    shell:
        """
        # Make bedgraph
        bedtools genomecov -ibam {input} -bg -5 -strand - | awk '{ printf "%s \t %d \t %d \t %d\n", $1,$2,$3,$4 }' > {output.reverse}
        """

### Get chromosome size information from bam file
rule get_chr_sizes:
    input:
        expand("results/sorted_bam/{sample_one}.bam", sample_one=CTL_SAMPLES[0]),
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
    log:
        "logs/callTSS/callTSS.log"
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
            --detail_file {output.detail} \
            --annotation_file {input.anno} \
            --call_method {params.call_method} \
            --bin_winner_size {params.bin_winner_size} \
            {input.fwd} {input.rev} {input.chrom} {output.bed} \
            &> {log}
        """


### Use my custom script to get TSS FASTA
rule make_TSSome:
    input:
        bed="results/callTSS/callTSS.tss.bed",
        fasta=config.get("genome_fasta"),
    output:
        fasta="results/make_TSSome/TSSome.fasta",
    params:
        Rscript=workflow.source_path("../scripts/createTSSome.R"),
        extra=config.get("createTSSome_extra")
    log:
        "logs/make_TSSome/makeTSSome.log"
    conda:
        "../envs/TSSome.yml"
    threads: 1
    shell:
        """
        chmod +x {params.Rscript}
        R {params.Rscript} \
            --fasta {input.fasta} \
            --bed {input.bed} \
            --output_fasta {output.fasta} \
            {params.extra} &> {log}
        """
