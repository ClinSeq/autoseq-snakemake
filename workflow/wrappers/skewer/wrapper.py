#!/usr/bin/env python

__author__ = "Sarath Murugan"
__copyright__ = "Copyright 2024, Sarath Murugan"
__maintainer__ = "Sarath Murugan"
__email__ = "sarath.murugan@outlook.com"

import os
import snakemake.shell as shell

extra = snakemake.params.get("extra", "")
tmpdir = snakemake.params.get("tmpdir", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """
    mkdir {tmpdir} 
    skewer -z -t {snakemake.threads} --quiet \
           -o {tmpdir}/skewer  \
           {snakemake.input.fq1} \
           {snakemake.input.fq2} 
    cp {tmpdir}/skewer-trimmed-pair1.fastq.gz {snakemake.output.out_fq1} 
    cp {tmpdir}/skewer-trimmed-pair2.fastq.gz {snakemake.output.out_fq2}
    rm -rf {tmpdir} 2> {log}
    """
)
