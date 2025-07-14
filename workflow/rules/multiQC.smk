### Compile all QC into one report
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_r{read}.{ext}",
            sample=SAMP_NAMES,
            read=READ_NAMES,
            ext=["html", "zip"],
        ),
    output:
        "results/multiqc/multiqc_report.html",
    params:
        # 1. raise / disable the size filter
        extra="--file-size-limit 0",          # 0 = no limit
        # 2. hand actual filenames directly to MultiQC
        use_input_files_only=True,
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml",
    script:
        "../scripts/multiQC.py"
