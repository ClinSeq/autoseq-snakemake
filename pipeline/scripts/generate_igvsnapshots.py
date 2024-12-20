#!/usr/bin/env python3

"""
    generate_igvsnapshots.py is a script designed to create IGV snapshots using a customized XML file.
    This tool allows for a quick review of variants, facilitating manual curation.

"""

__author__      = "Sarath Kumar Murugan"
__copyright__   = "Copyright 2024, Sarath Murugan"
__email__       = "sarath.murugan@outlook.com"

import os
import argparse
import logging
import subprocess

_default_ops = ["new", "genome hg19"]


def is_indel(variant):
    """
    Check if a variant is an indel.
    """
    ref = variant[3]
    alt = variant[4]

    if len(ref) > 1:
        return True
    
    for _alt in alt:
        if len(_alt) != len(ref):
            return True

    return False


def adj_base(variant):
    """
    
    """
    ref = variant[3]
    alt = variant[4]
    end = int(variant[2])

    for _alt in alt:
        if _alt is None:
            return None
        if len(_alt) < len(ref):
            return round(end + len(ref)/2)
        else:
            return end
    return None


def create_batch_script(vars, xml, outdir):
    """
    """
    cmds = list()
    cmds.extend(_default_ops)
    cmds.append(f"snapshotDirectory {outdir}")
    cmds.append(f"load {xml}")
    with open(vars, 'r') as fh:
        header = fh.readline()
        for line in fh:
            line = line.strip().split('\t')
            _chr = line[0] 
            _start = line[1]
            _end = line[2]
            if is_indel(line):
                # small deletion 
                cmds.append(f"goto chr{_chr}:{str(adj_base(line))}")
            else:
                cmds.append(f"goto chr{_chr}:{_end}")
            cmds.append("sort base")
            cmds.append(f"snapshot chr{_chr}_{_start}-{_end}.png")
    
    cmds.append("exit")

    batch_script_path = outdir + "/igv_batch.sh"
    with open(batch_script_path, 'w') as fh:
        fh.write("\n".join(cmds))

    return(batch_script_path)


def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
            format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.info("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})


def main():
    parser = argparse.ArgumentParser(description=
        'Generate IGVsnapshots for single nucleotide variants')
    parser.add_argument('--xml', required=True, help="IGV session xml file ")
    parser.add_argument('-i', '--input', required=True, help="SNVs txt file as input")
    parser.add_argument('-p', '--preferences', required=True, help="Preferences for IGV batch")
    parser.add_argument('--output', help="output directory to store all images")
    args = parser.parse_args()
    
    setup_logging()

    if not os.path.exists(args.input) or not os.path.exists(args.xml) :
        logging.error(FileNotFoundError)
        exit
    
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
    
    batch_script = create_batch_script(args.input, args.xml, args.output)
    igv_cmd = f"igv.sh -b {batch_script} --preferences {args.preferences} " 

    try:
        subprocess.run(igv_cmd, shell=True)
    except Exception as err:
        logging.error(err)


if __name__ == "__main__":
    main()