import os
import sys
import json
import yaml
from rich.console import Console
from rich.table import Column, Table
import click
import subprocess

from loguru import logger as Log
import pipeline
from pipeline.utils.utils import make_paths_absolute, Pipeline
from pipeline.utils.clinseq_barcodes import data_available_for_clinseq_barcode, \
    extract_clinseq_barcodes, validate_clinseq_barcodes, convert_barcodes_to_sampledict, \
    check_sampledata, normpath


def console_autoseq():
    console = Console()
    console.print("[magenta]     _         _       ____             ")
    console.print("[magenta]    / \  _   _| |_ ___/ ___|  ___  __ _ ")
    console.print("[magenta]   / _ \| | | | __/ _ \___ \ / _ \/ _` |")
    console.print("[magenta]  / ___ \ |_| | || (_) |__) |  __/ (_| |")
    console.print("[magenta] /_/   \_\__,_|\__\___/____/ \___|\__, |")
    console.print("[magenta]                                     |_| :snake:")
    console.print("                         version: {}".format(pipeline.__version__))
    
    console.print("\n")

    # run autoseq-cli
    cli()


@click.group()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option("-v", "--verbose", is_flag=True, default=False, help="Print verbose output to the console.")
@click.pass_context
def cli(context, loglevel, verbose):
    """
    Autoseq - pipeline

    Autoseq consists of a custom-pipeline with additional support modules aimed 
    primarily for the analysis of data from high-throughput sequencing of liquid biopsies.
    """
    Log.remove()
    Log.add(sys.stdout, colorize=True, 
            format="<green>{time:YYYY-MM-DD at HH:mm:ss}</green> | {level} | <level>{message}</level>", 
            level=loglevel)


@cli.command()
@click.option("--outdir", help="output directory")
@click.argument('barcodes-file', type=str)
@click.pass_context
def config(context, barcodes_file, outdir):
    """
    Create sample json for given clinseq barcodes
    """
    clinseq_barcodes = extract_clinseq_barcodes(barcodes_file)
    validate_clinseq_barcodes(clinseq_barcodes)

    sample_dict = convert_barcodes_to_sampledict(clinseq_barcodes)

    for sdid in sample_dict:    
        fn = "{}/{}.json".format(outdir, sdid)
        with open(fn, 'w') as f:
            json.dump(sample_dict[sdid], f, sort_keys=True, indent=4)
            Log.info(f"Autoseq samples config file created - {fn}")


@cli.command()
@click.pass_context
def list(context):
    """
    List autoseq available pipelines with version
    """
    console = Console()

    pipelines = Table(show_header=True, header_style="bold magenta")
    pipelines.add_column("No.")
    pipelines.add_column("Pipeline Name")
    pipelines.add_column("Type")
    pipelines.add_column("Last update")
    
    pipelines.add_row("1", "Autoseq", "Targeted Re-sequencing", "May 25 2021")

    console.print(pipelines)
    

@cli.command()
@click.option("--ref", '-r', help="json file with reference files to use", 
            type=click.Path(exists=True))
@click.option("--samples", help="json file contains list of samples")
@click.option("--outdir", default=os.getcwd() ,help="output directory")
@click.option("--libdir", help="directory to search libraries")
@click.option("--configfile", help="configuration file for params")
@click.option("--scratch", default="/tmp", help="path to /tmp/scratch")
@click.option("--dryrun/--run", default=False, help=" --dryrun for testing snakemake workflow")
@click.option("--umi", is_flag=True, help="To process the data with UMI- Unique Molecular Identifier")
@click.option("--profile", default='shell', help="job schedulers eg. SLURM")
@click.option("--pipeline", default='autoseq', help="Pipeline to be launched")
@click.option("-n", "--normal-bam", default=None, help="Normal bam files dir, Applicable only to tumor only pipeline")
@click.option("--use-singularity", is_flag=True, help="To use singularity")
@click.option("--singularity", help="Path to singularity image")
@click.option("--smk-opt", help="snakemake option")
@click.option("--cores", help="max number of cores")
@click.pass_context
def launch(context, ref, samples, outdir, libdir, 
            configfile, scratch, dryrun, umi, profile, 
            pipeline, normal_bam, use_singularity, 
            singularity, cores, smk_opt):
    """
    launch the respective pipeline with samples json 
    """
    sample_json = json.load(open(samples))
    sdid = sample_json['sdid']

    Log.info(f"Looking for fastq files {sdid} in {libdir}")
    sampledata, all_clinseq_barcodes = check_sampledata(libdir, sample_json)
    Log.info(f"Libraries {all_clinseq_barcodes} have data. Using it.")

    if not configfile:
        tool_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        configfile = os.path.join(tool_dir, 'config.yml')
    
    config_dict = yaml.load(open(configfile), Loader=yaml.FullLoader)

    normal_barcode = [i for i in all_clinseq_barcodes if '-N-' in i]
    tumor_barcode = [i for i in all_clinseq_barcodes if '-T-' in i or '-CFDNA-' in i]

    sample_str = "_".join(tumor_barcode + normal_barcode)
    outdir = os.path.join(outdir, sampledata['sdid'], sample_str)

    if use_singularity and not singularity:
        Log.error('Singularity file path is not available, use --singularity for file path')
        raise click.Abort()
    
    if use_singularity:
        if not os.path.exists(os.path.join(singularity, "autoseq-smk.sif")) or \
            not os.path.exists(os.path.join(singularity, "gridss.sif")):
            Log.error('Singularity file does not exist !!')
            raise click.Abort()

    # config dictionary update
    config_dict['samples'] = normpath(samples)
    config_dict['reference'] = normpath(ref)
    config_dict['outdir'] = normpath(outdir)
    config_dict['libdir'] = normpath(libdir)
    config_dict['umi'] = umi
    config_dict['container']['base'] = os.path.join(singularity, "autoseq-smk.sif") if use_singularity else ' '
    config_dict['container']['gridss'] = os.path.join(singularity, "gridss.sif") if use_singularity else ' '
    config_dict['container']['franken'] = os.path.join(singularity, "franken.sif") if use_singularity else ' '
    
    # pipeline based args
    if pipeline == 'tumor_only' and normal_bam:
        nClip_bam = os.path.join(normal_bam, normal_barcode[0] + "_clipoverlap.bam")
        nNodups_bam = os.path.join(normal_bam , normal_barcode[0] + "_nodups.bam")
        for bam in [nClip_bam, nNodups_bam]:
            if os.path.isfile(bam):
                Log.info(f"Normal sample bam file - {bam}")
            else:
                Log.error(f"{bam} does not exist")
                raise click.Abort()

        config_dict['normal_bams'] = [nClip_bam, nNodups_bam]
        

    out_configpath = os.path.join(normpath(outdir), f"config_{sample_str}.yml")
    jobdb = os.path.join(normpath(outdir), f"{sample_str}.jobdb")

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    with open(out_configpath, 'w') as cf:
        yaml.safe_dump(config_dict, cf, default_flow_style=False)
    
    bind_paths = set()
    if use_singularity:
        for path in ('samples', 'reference', 'outdir', 'libdir'):
            bind_paths.add(os.path.dirname(config_dict[path]))
        
        bind_paths.add(os.path.dirname(os.path.dirname(config_dict['reference'])))


    if pipeline == "tumor_only":
        snakefile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tumor_only/Snakefile')
    else:
        snakefile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'autoseq/Snakefile')

    autoseq = Pipeline(snakefile = snakefile, 
                      config = out_configpath, 
                      sdid = sdid, 
                      workdir = outdir, 
                      dryrun = dryrun, 
                      profile = profile,
                      jobdb = jobdb,
                      smk_option = smk_opt,
                      use_singularity = use_singularity,
                      bind_paths = bind_paths,
                      cores = cores)
    
    cmd = autoseq.build_cmd()

    Log.info(f"Launching autoseq - {pipeline} pipeline ...")
    # print(cmd) 
    try:
        subprocess.run(cmd, shell=True)
    except Exception as err:
        Log.error(err)


if __name__ == "__main__":
    console_autoseq()
