#!/usr/bin/env python

__author__ = "Sarath Murugan"
__copyright__ = "Copyright 2024, Sarath Murugan"
__maintainer__ = "Sarath Murugan"
__email__ = "sarath.murugan@outlook.com"

import os
import snakemake.shell as shell

extra = snakemake.params.get("extra", "")
tmpdir = snakemake.params.get("tmpdir", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

bwa_index = snakemake.input.bwa_index


if not isinstance(bwa_index, str):
    raise ValueError(
        "BWA Index must be given in the Inputs ! "
    )


shell(
    "picard SamToFastq I={snakemake.input.bam} F=/dev/stdout "
        "   INTERLEAVE=true TMP_DIR={tmpdir} "
        " | bwa mem -p -t {snakemake.threads} {bwa_index} /dev/stdin " 
        " | picard -Djava.io.tmpdir={tmpdir}  "
        " {java_options} MergeBamAlignment "
        " UNMAPPED={snakemake.input.bam} "
        " ALIGNED=/dev/stdin "
        " O={snakemake.output.bam}"
        " R={bwa_index} "
        " SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true "
        " MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true "
        " TMP_DIR={tmpdir} 2> {log} && rm -rf {tmpdir} "
)


