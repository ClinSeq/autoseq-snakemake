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
