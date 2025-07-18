### Compile all QC into one report
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_r{read}_fastqc.{ext}",
            sample=SAMP_NAMES,
            read=READ_NAMES,
            ext=["html", "zip"],
        ),
    output:
        "results/multiqc/multiqc_report.html",
    params:
        extra=config.get("multiqc_extra", ""),
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    script:
        "../scripts/multiQC.py"
