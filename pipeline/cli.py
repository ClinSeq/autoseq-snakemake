import os
import json
import yaml
import click
import subprocess


from pipeline.rules.utils import make_paths_absolute, Pipeline
from pipeline.rules.clinseq_barcodes import data_available_for_clinseq_barcode, \
    extract_clinseq_barcodes, validate_clinseq_barcodes, convert_barcodes_to_sampledict, \
    check_sampledata, normpath


@click.group()
@click.pass_context
def cli(context):
    pass


@cli.command()
@click.option("--outdir", help="output directory")
@click.argument('barcodes-file', type=str)
@click.pass_context
def config(context, barcodes, outdir):
    clinseq_barcodes = extract_clinseq_barcodes(barcodes)
    validate_clinseq_barcodes(clinseq_barcodes)

    sample_dict = convert_barcodes_to_sampledict(clinseq_barcodes)

    for sdid in sample_dict:    
        fn = "{}/{}.json".format(outdir, sdid)
        with open(fn, 'w') as f:
            json.dump(sample_dict[sdid], f, sort_keys=True, indent=4)


@cli.command()
@click.option("--ref", '-r', help="json file with reference files to use", 
            type=click.Path(exists=True))
@click.option("--samples", help="json file contains list of samples")
@click.option("--outdir", default=os.getcwd() ,help="output directory")
@click.option("--libdir", help="directory to search libraries")
@click.option("--configfile", help="configuration file for params")
@click.option("--scratch", default="/tmp", help="path to /tmp/scratch")
@click.option("--dryrun/--run", default=False)
@click.option("--cores", help="max number of cores")
@click.pass_context
def launch(context, ref, samples, outdir, libdir, configfile, scratch, dryrun, cores):
    # samples
    sample_json = json.load(open(samples))

    # check sample data
    sampledata, all_clinseq_barcodes = check_sampledata(libdir, sample_json)

    if not configfile:
        tool_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        configfile = os.path.join(tool_dir, 'config.yml')
    
    config_dict = yaml.load(open(configfile), Loader=yaml.FullLoader)

    sample_str = "_".join(all_clinseq_barcodes)
    outdir = os.path.join(outdir, sampledata['sdid'], sample_str)

    # update config dict
    config_dict['samples'] = normpath(samples)
    config_dict['reference'] = normpath(ref)
    config_dict['outdir'] = normpath(outdir)
    config_dict['libdir'] = normpath(libdir)
    out_configpath = os.path.join(normpath(outdir), 'config.yml')

    # create output dir
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    
    # write config file inside outdir
    if not os.path.exists(out_configpath):
        with open(out_configpath, 'w') as cf:
            yaml.safe_dump(config_dict, cf, default_flow_style=False)

    snakefile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Snakefile')

    # building snakemake pipeline 
    liqbio = Pipeline(snakefile, out_configpath, outdir, dryrun=dryrun)
    cmd = liqbio.build_cmd()

    try:
        subprocess.run(cmd, shell=True)
    except Exception as err:
        click.echo(err)





    
    
    
    
