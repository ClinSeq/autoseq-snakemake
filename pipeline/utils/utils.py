import os, re
from pipeline.utils.clinseq_barcodes import parse_prep_id, compose_sample_str, \
    extract_unique_capture, find_fastqs


def get_containers(_path):
    """
    """
    containers = {
        "base": os.path.join(_path, "autoseq-base.sif"),
        "franken": os.path.join(_path, "autoseq-franken.sif"),
        "gatk3": os.path.join(_path, "autoseq-gatk3.sif"),
        "gridss": os.path.join(_path, "autoseq-gridss.sif"),
        "jumble": os.path.join(_path, "autoseq-jumble.sif"),
        "purecn": os.path.join(_path, "autoseq-purecn.sif"),
        "ensemblvep": os.path.join(_path, "autoseq-ensemblvep.sif"),
        "somaticseq": os.path.join(_path, "autoseq-somaticseq.sif"),
        "svcaller": os.path.join(_path, "autoseq-svcaller.sif")
    }

    for k, v in containers.items():
        if not os.path.exists(v):
            raise ValueError("Invalid container PATH to " + v)
    
    return containers


def get_scheduler(scheduler, filetype):
    """
    In cluster environment, to get sheduler script and config file

    params: scheduler type, filetype
    return: file (script or config) path
    """
    tool_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    submit_fp = ""

    if scheduler == 'slurm' and filetype == 'pyscript':
        submit_fp = "scheduler/slurm_submit.py"
    
    if scheduler == 'slurm' and filetype == 'config':
        submit_fp = "scheduler/cluster_config.json"

    script  = os.path.join(tool_dir, submit_fp)

    if not os.path.isfile(script):
        raise FileNotFoundError(script)
    
    return script


class Pipeline:
    """
    Class pipeline to build snakmake command based on given args.

    """
    def __init__(self, snakefile, config, cluster_config, sdid, project_id, workdir, dryrun, 
                profile, jobdb, smk_option, use_singularity, bind_paths, cores='4'):
        self.snakefile = snakefile
        self.cores = cores
        self.configfile = config
        self.cluster_config = cluster_config
        self.sdid = sdid
        self.project_id = project_id
        self.workdir = workdir
        self.profile = profile
        self.jobdb = jobdb
        self.dryrun = dryrun
        self.smk_option = smk_option
        self.use_singularity = use_singularity
        self.bind_paths = bind_paths
    
    def build_cmd(self):
        dryrun = ''
        profile_cmd = ''
        slurm_cmd = ''
        smk_opt = ''
        singularity_cmd = ''
        
        if self.dryrun:
            dryrun = "-n"

        if self.smk_option:
            smk_opt = self.smk_option

        if self.profile == 'slurm':
            if self.cluster_config:
                cluster_config = self.cluster_config
            else:
                cluster_config = get_scheduler(self.profile, 'config')
            
            slurm_submit = get_scheduler(self.profile, 'pyscript')
            slurm_cmd = " --notemp --immediate-submit -j 500 "
            slurm_cmd += " --jobname smk.{{rulename}}.{}-{}.{{jobid}}.sh ".format(self.project_id, self.sdid)
            slurm_cmd += " --cluster-config {} ".format(cluster_config)
            slurm_cmd += (" --cluster '{} "
                          " --jobdb {} "
                          " {{dependencies}} '".format(slurm_submit, self.jobdb))

        if self.use_singularity:
            singularity_cmd = " --use-singularity "
            singularity_cmd += " --singularity-args ' "

            for path in self.bind_paths:
                singularity_cmd += f" --bind {path}:{path} "
            
            singularity_cmd += "'"


        cmd = ("snakemake -p --snakefile {} "
               " --cores {} "
               " --directory {} "
               " --configfile {} "
               " {} {} {} {} {} ").format(self.snakefile,
                              self.cores,
                              self.workdir,
                              self.configfile,
                              dryrun,
                              profile_cmd,
                              singularity_cmd,
                              slurm_cmd,
                              smk_opt)
        return cmd


class SinglePanelResults():
    def __init__(self):
        self.bamfile = None
        self.umibam = None
        
        # CNV kit outputs:
        self.cnr = None
        self.cns = None
        self.seg = None

        # Coverage QC call:
        self.cov_qc_call = None

        # Structural variants, organised as a dictionary with event type as key,
        # and their effects:
        self.svs = {}
        self.sv_effects = None

        # FIXME: Msings should never be run for normal samples => OO progr. fail. Refactor.
        # Msings output:
        self.msings_output = None


def get_fqwildcards(sample_barcode, libdir):
    """
    function to extract fastq prefix and suffix

    param: sample barcode
    param: library directory
    return: fastq prefix, suffix for R1 and R2
    """
    fq1_files, fq2_files = find_fastqs(sample_barcode, libdir)
    fq1_abs = [os.path.basename(x) for x in fq1_files]
    fq2_abs = [os.path.basename(x) for x in fq2_files]
    fq_prefix = list()

    regex_fq1 = r'(.+)(_1.fastq.gz|_1.fq.gz|R1_\d{3}.fastq.gz)'
    regex_fq2 = r'(.+)(_2.fastq.gz|_2.fq.gz|R2_\d{3}.fastq.gz)'
    s1 = ''
    
    for fq in fq1_abs:
        _fq_ = [i for i in re.split(regex_fq1, fq) if i != '']
        fq_prefix.append(_fq_[0])
        s1 = _fq_[1]

    _fq_ = [i for i in re.split(regex_fq2, fq2_abs[0]) if i != '']
    s2 = _fq_[1]

    return fq_prefix, s1, s2


def get_capture_bam(unique_capture, bamfiles):
    """
    return bamfiles for given unique capture
    """
    sample_str = compose_sample_str(unique_capture)

    for bam in bamfiles:
        filename = os.path.basename(bam)
        if sample_str in filename:
            return bam
    
    return False


def get_cnvkitref(wildcards, reference):
    """
    return cnvkit reference file
    """
    unique_capture = extract_unique_capture(wildcards.sample)
    capture_name = get_capture_name(unique_capture.capture_kit_id)
    library_name = get_prep_kit_name(unique_capture.library_kit_id)
    sample_type = unique_capture.sample_type

    cnvkit_ref = None
    if 'cnvkit-ref' in reference['targets'][capture_name]:
        cnvkit_ref = list(list(reference['targets'][capture_name]['cnvkit-ref'].values())[0].values())[0]
    
    try:
        cnvkit_ref = reference['targets'][capture_name]['cnvkit-ref'][library_name][sample_type]
    except KeyError:
        pass

    return cnvkit_ref


def get_jumbleref(wildcards, reference):
    """
    return jumble reference file
    """
    unique_capture = extract_unique_capture(wildcards.sample)
    capture_name = get_capture_name(unique_capture.capture_kit_id)

    jumble_ref = None
    if 'jumble-ref' in reference['targets'][capture_name]:
        jumble_ref = reference['targets'][capture_name]['jumble-ref']

    return jumble_ref


def get_capture_svs(wildcards, outdir):
    """
    return gtfs dictionary for given sample
    """
    events = ["DEL", "DUP", "INV", "TRA"]
    gtfs = dict()
    for event in events:
        gtfs[event] = outdir + "/svs/svcaller/{}-{}.gtf".format(wildcards.sample, event)

    return gtfs


def get_svcaller_mut(wildcards, outdir):
    """
    return svcaller mut file name
    """
    sample_str = compose_sample_str(extract_unique_capture(wildcards.sample))
    
    return "{}/svs/igv/{}_svcaller.mut".format(outdir, sample_str)


def get_readgroup(wildcards):
    """
    return readgroup for alignments
    """
    try:
        sample = wildcards.sample
    except AttributeError:
        sample = wildcards

    library_id = parse_prep_id(sample)
    sample_string = compose_sample_str(extract_unique_capture(sample))

    readgroup = "\"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:ILLUMINA\"".format(\
        rg_id=sample, rg_sm=sample_string, rg_lb=library_id)
    
    return readgroup


def get_targets(wildcards, reference, key):
    """
    return bed file corresponds to capture id
    """
    unique_capture = extract_unique_capture(wildcards.sample)
    targets = get_capture_name(unique_capture.capture_kit_id)

    if unique_capture.capture_kit_id in ["P2", "S2", "B2"]:
        return reference['small-design'][targets][key]
    
    return reference['targets'][targets][key]


def get_target_name(wildcards):
    """
    return capture id
    """
    unique_capture = extract_unique_capture(wildcards.sample)
    targets = get_capture_name(unique_capture.capture_kit_id)
    
    return targets


def get_chromosomes(targets):
    """
    extract chromosomes from target bed file 
    """
    chromos = set()
    with open(targets, 'r') as bedfile:
        for line in bedfile.readlines():
            chromos.add(line.split('\t')[0])
    
    return chromos
            

def get_target_region(wildcards, chrsizes):
    """
    utility function to pass target region param to indelrealigner
    
    return: target_region eg: 1:1-122121212
    """
    
    chromo = wildcards.chr

    if chromo in chrsizes:
        return ":".join([chromo, chrsizes[chromo]])

    raise KeyError(chromo)


def get_capture_name(capture_kit_code):
    """
    Convert a two-letter capture kit code to the corresponding capture kit name.

    :param capture_kit_code: The two-letter capture kit code.
    :return: The capture-kit name.
    """
    
    # FIXME: Move this information to a config JSON file.
    capture_kit_loopkup = {"CS": "clinseq_v3_targets",
                            "CZ": "clinseq_v4",
                            "EX": "EXOMEV3",
                            "EO": "EXOMEV1",
                            "RF": "fusion_v1",
                            "CC": "core_design",
                            "CD": "discovery_coho",
                            "CB": "big_design",
                            "AL": "alascca_targets",
                            "TT": "test-regions",
                            "CP": "progression",
                            "CM": "monitor",
                            "PC": "probio_comprehensive",
                            "PB": "probio_biomarker_signature",
                            "PA": "pancancer",
                            "C2": "probio_comprehensive2",
                            "C3": "probio_comprehensive3",
                            "C4": "probio_comprehensive4",
                            "PN": "pancancer2",
                            "PE": "pancancer2_enzymatic",
                            "P2": "probio_biomarkersignature2",
                            "S2": "probio_biomarkersignature2",
                            "B2": "probio_biomarkersignature2",
                            "PS": "probio_snvindel",
                            "S3": "probio_biomarkersignature2" # S3 pointed to S2 capture files
                            }

    if capture_kit_code == 'WG':
        return 'lowpass_wgs'

    else:
        return capture_kit_loopkup[capture_kit_code]


def get_prep_kit_name(prep_kit_code):
    """
    Convert a two-letter library kit code to the corresponding library kit name.

    :param prep_kit_code: Two-letter library prep code. 
    :return: The library prep kit name.
    """

    # FIXME: Move this information to a config JSON file.
    prep_kit_lookup = {"BN": "BIOO_NEXTFLEX",
                        "KH": "KAPA_HYPERPREP",
                        "TD": "THRUPLEX_DNASEQ",
                        "TP": "THRUPLEX_PLASMASEQ",
                        "TF": "THRUPLEX_FD",
                        "TS": "TRUSEQ_RNA",
                        "NN": "NEBNEXT_RNA",
                        "VI": "VILO_RNA"}

    return prep_kit_lookup[prep_kit_code]


def make_paths_absolute(input_dict, base_path):
    """Processes the input dictionary, converting relative file paths to absolute
    file paths throughout the dictionary structure.

    Specifically, for each value in the dictionary:
    - If it is also a dictionary, then recursively apply this function,
    replacing the initial dictionary.
    - Otherwise:
    -- If the value is a non-null string that is not already an absolute path,
    then try prepending the specified base_path and see if the resulting file name
    exists, and in that case then replace the string with the resulting absolute path.
    """

    for curr_key, curr_value in input_dict.items():
        if isinstance(curr_value, dict):
            input_dict[curr_key] = make_paths_absolute(curr_value, base_path)
        else:
            converted_value = convert_to_absolute_path(curr_value, base_path)
            input_dict[curr_key] = converted_value

    return input_dict


def convert_to_absolute_path(possible_relative_path, base_path):
    """
    Convert the input potential relative file path to an absolute path by
    prepending the specified base_path, but only if the resulting absolute path points
    to a pre-existing file or directory.

    If the base_path cannot be prepended, then simply return the original input value.

    :param possible_relative_path: A string potentially indicating a relative file/directory path. 
    :param base_path: The base path to prepend.
    :return: Modified path string.
    """

    converted_value = possible_relative_path
    try:
        if not os.path.isabs(possible_relative_path):
            joined_path = os.path.join(base_path, possible_relative_path)
            if os.path.isfile(joined_path) or os.path.isdir(joined_path):
                converted_value = joined_path

    except Exception:
        pass

    return converted_value
 
