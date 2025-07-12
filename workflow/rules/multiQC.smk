### Compile all QC into one report
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_r{read}.{ext}",
            sample=SAMPLE_NAMES,
            read=["1", "2"],
            ext=["html", "zip"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra=config.get("multiqc_extra", ""),
    log:
        "results/multiqc/multiqc.log",
    wrapper:
        "v6.0.0/bio/multiqc"