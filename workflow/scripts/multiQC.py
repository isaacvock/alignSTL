#!/usr/bin/env python
"""Snakemake wrapper for running MultiQC.

This version fixes two issues in the stock wrapper:

1. `input_data` is now rendered as a *space-separated, shell-quoted*
   list of paths instead of a Python set literal, so MultiQC actually
   receives the files / directories it needs.

2. Every pathname handed to the shell command is wrapped in
   `shlex.quote()` to survive spaces or other special characters.

The wrapper remains fully generic: any additional command-line
switches can still be injected with the Snakemake parameter
`params.extra`, and you can switch between feeding MultiQC whole
directories or individual files with the boolean
`params.use_input_files_only`.
"""

__author__      = "Julian de Ruiter (edited by ChatGPT)"
__copyright__   = "Copyright 2017, Julian de Ruiter"
__email__       = "julianderuiter@gmail.com"
__license__     = "MIT"

from pathlib import Path
from shlex   import quote
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import is_arg

# ----------------------------------------------------------------------
# Gather Snakemake inputs
# ----------------------------------------------------------------------
extra = snakemake.params.get("extra", "")
log   = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Configuration files supplied as normal Snakemake inputs
mqc_config = snakemake.input.get("config", "")
if isinstance(mqc_config, list):
    for fp in mqc_config:
        extra += f" --config {quote(fp)}"
    mqc_config = set(mqc_config)
elif mqc_config:
    extra += f" --config {quote(mqc_config)}"
    mqc_config = {mqc_config}
else:
    mqc_config = set()

# ----------------------------------------------------------------------
# Decide which paths to give MultiQC
# ----------------------------------------------------------------------
use_input_files_only = snakemake.params.get("use_input_files_only", False)

if use_input_files_only:
    # Hand MultiQC each file explicitly
    input_paths = {Path(fp) for fp in snakemake.input if fp not in mqc_config}
else:
    # Hand MultiQC the *parent directories* to let it discover files
    input_paths = {Path(fp).parent for fp in snakemake.input if fp not in mqc_config}

# Convert the Python set to a proper command-line argument string
input_data = " ".join(quote(str(p)) for p in sorted(input_paths))

# ----------------------------------------------------------------------
# Adjust command-line flags based on expected outputs
# ----------------------------------------------------------------------
no_report = True
for output in snakemake.output:
    if output.endswith(".html"):
        no_report = False
    if output.endswith("_data"):
        extra += " --data-dir"
    if output.endswith(".zip"):
        extra += " --zip-data-dir"

if no_report:
    extra += " --no-report"
if (
    not is_arg("--data-dir", extra)
    and not is_arg("-z", extra)
    and not is_arg("--zip-data-dir", extra)
):
    extra += " --no-data-dir"

# ----------------------------------------------------------------------
# Derive output directory and file name from Snakemake's first output
# ----------------------------------------------------------------------
out_dir   = Path(snakemake.output[0]).parent
file_name = Path(snakemake.output[0]).with_suffix("").name

# ----------------------------------------------------------------------
# Launch MultiQC
# ----------------------------------------------------------------------
shell(
    "multiqc"
    " {extra}"
    " --outdir {quote(str(out_dir))}"
    " --filename {quote(file_name)}"
    " {input_data}"
    " {log}"
)
