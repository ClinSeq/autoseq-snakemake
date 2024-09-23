#!/usr/bin/env python

__author__ = "Sarath Murugan"
__copyright__ = "Copyright 2024, Sarath Murugan"
__maintainer__ = "Sarath Murugan"
__email__ = "sarath.murugan@outlook.com"

import os
import snakemake.shell as shell

extra = snakemake.params.get("extra", "")
tmpdir = snakemake.params.get("tmpdir", "")
java_options = snakemake.params.get("java_options", "")

## Input reference files
reference_genome = snakemake.input.reference_genome
target_region = snakemake.input.target_region
known_1kg = snakemake.input.known_1kg
known_mills_1kg = snakemake.input.known_mills_gs

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


shell(
    "gatk3 {java_options} "
    " -Djava.io.tmpdir={tmpdir} "
    " -T RealignerTargetCreator "
    " -R {reference_genome} "
    " -known {known_1kg} "
    " {extra} "
    " -L {target_region} "
    " -known {known_mills_gs} "
    " -I {snakemake.input.bam} "
    " -o {snakemake.output.target_intervals} 2> {log} "
    " && rm -rf {tmpdir} "
)
