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
    Adjust bases in small indels
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


def create_batch_script_snv(vars, xml, outdir):
    """
    Create igv batch script for SNVs
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


def adj_window(start, end, d):
    """
    Adjust window size based on the distance between the variant and the end of the window
    """
    diff = end - start
    center = start + round(diff/2)
    if d == '+':
        return center + 150
    else:
        return center - 150



def create_batch_script_sv(vars, xml, outdir):
    """
    Create igv batch script for SVs
    """
    cmds = list()
    cmds.extend(_default_ops)
    cmds.append(f"snapshotDirectory {outdir}")
    cmds.append(f"load {xml}")
    cmds.append("squish tb_simpleRepeat ")
    with open(vars, 'r') as fh:
        header = fh.readline()
        for line in fh:
            data = line.strip().split('\t')
            # skip GSR not in curator
            if data[12] == "NO":
                continue
            chr_a = data[0] 
            start_a = int(data[1])
            end_a = int(data[2])
            chr_b = data[3] 
            start_b = int(data[4])
            end_b = int(data[5])
            if data[10] == "gridss":
                bp1_start = end_a - 150
                bp1_end = end_a + 150
                bp2_start = end_b - 150
                bp2_end = end_b + 150
            else:
                bp1_start = start_a if (end_a - start_a < 300) else adj_window(start_a, end_a, '-')
                bp1_end = end_a if (end_a - start_a < 300) else adj_window(start_a, end_a, '+')
                bp2_start = end_b if (end_b - start_b < 300) else adj_window(start_b, end_b, '-')
                bp2_end = end_b if (end_b - start_b < 300) else adj_window(start_b, end_b, '+')
            
            cmds.append(f"goto chr{chr_a}:{bp1_start}-{bp1_end}")
            cmds.append("sort base")
            cmds.append(f"snapshot chr{chr_a}_{start_a}-{end_a}.png")
            cmds.append(f"goto chr{chr_b}:{bp2_start}-{bp2_end}")
            cmds.append("sort base")
            cmds.append(f"snapshot chr{chr_b}_{start_b}-{end_b}.png")
    
    cmds.append("exit")

    batch_script_path = outdir + "/igv_batch_sv.sh"
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
    parser.add_argument('-i', '--input', required=True, help="SNVs or SVs txt file as input")
    parser.add_argument('-t', '--type', choices = ['snv', 'sv'], help="Variant type SNV or SV")
    parser.add_argument('-p', '--preferences', required=True, help="Preferences for IGV batch")
    parser.add_argument('--output', help="output directory to store all images")
    args = parser.parse_args()
    
    setup_logging()

    if not os.path.exists(args.input) or not os.path.exists(args.xml) :
        logging.error(FileNotFoundError)
        exit
    
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
    
    if args.type == 'snv':
        batch_script = create_batch_script_snv(args.input, args.xml, args.output)
    else:
        batch_script = create_batch_script_sv(args.input, args.xml, args.output)

    igv_cmd = f"igv.sh -b {batch_script} --preferences {args.preferences} " 

    try:
        subprocess.run(igv_cmd, shell=True)
    except Exception as err:
        logging.error(err)


if __name__ == "__main__":
    main()