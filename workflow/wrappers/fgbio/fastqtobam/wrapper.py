#!/usr/bin/env python

__author__ = "Sarath Murugan"
__copyright__ = "Copyright 2024, Sarath Murugan"
__maintainer__ = "Sarath Murugan"
__email__ = "sarath.murugan@outlook.com"

import os
import snakemake.shell as shell

extra = snakemake.params.get("extra", "")
java_options = snakemake.params.get("java_options", "")
sample = snakemake.params.get("sample", "")
tmpdir = snakemake.params.get("tmpdir", "")
read_struct = snakemake.params.get("read_struct", "")
library = snakemake.params.get("library", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "fgbio {java_options} "
    "   --tmp-dir {tmpdir} FastqToBam "
    "    -i {snakemake.input.fq1}     "
    "       {snakemake.input.fq2}      "
    "    -o {snakemake.output.bam}     "
    "    --sample {sample}             "
    "    --library {library}           "
    "    -r {read_struct}              "
    "    -s true  {log}  &&  "
    "rm -rf {tmpdir} "
)