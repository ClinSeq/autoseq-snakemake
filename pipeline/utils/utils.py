import os, re
import glob
import sys
import subprocess
import requests
import urllib3
from urllib3.exceptions import InsecureRequestWarning
from pipeline.settings import PROJECTS_CODE_TO_NAME
from pipeline.utils.clinseq_barcodes import parse_prep_id, compose_sample_str, \
    extract_unique_capture, find_fastqs, compose_lib_capture_str



def check_project_exists(project_id, GET_PROJECT_LIST_API):
    """
    Check if the project exists in the project list.
    """
    urllib3.disable_warnings(InsecureRequestWarning)
    response = requests.get(GET_PROJECT_LIST_API, verify = False)
    response.raise_for_status()
    all_projects = response.json()

    for proj in all_projects['data']:
        if proj['project_name'] == project_id:
            return True
    return False


def create_project_info(project_id, CREATE_PROJECT_INFO_API):
    """
    Create a new project info entry if the project does not exist.
    """
    if project_id not in PROJECTS_CODE_TO_NAME:
        return {"status": "error", "message": f"Project ID '{project_id}' is not recognized."}

    params = {
        "project_name": PROJECTS_CODE_TO_NAME[project_id],
        "prefix_name": project_id,
        "nfs_path": "",
        "proj_status": 0,
        "pdf_report": 1,
        "mtbp_report": 1,
        "sort_order": 1 
    }

    urllib3.disable_warnings(InsecureRequestWarning)
    response = requests.post(CREATE_PROJECT_INFO_API, params = params, verify=False)
    response.raise_for_status()
    
    if response.status_code != 200:
        return {"status": "error", "message": "Failed to create project info."}
    
    return {"status": "success", "message": f"Project '{project_id}' created successfully."}


def update_sample_info(project_id, sdid, capture, status, pipeline):
    """
    Validate the project ID by checking if it exists in the project list.
    """
    BASE_URL = os.environ.get("CURATOR_BASE_URL")  # Replace with actual API base URL
    GET_PROJECT_LIST = f"{BASE_URL}/all_project_list"
    CREATE_SAMPLE_INFO = f"{BASE_URL}/create-sample-info"
    CREATE_PROJECT_INFO = f"{BASE_URL}/create-project-info"
    
    if project_id not in PROJECTS_CODE_TO_NAME:
        return {"status": "error", "message": f"Project ID '{project_id}' is not recognized."}

    proj_name = PROJECTS_CODE_TO_NAME[project_id]
    try:
        if not check_project_exists(proj_name, GET_PROJECT_LIST):
            create_response = create_project_info(project_id, CREATE_PROJECT_INFO)
            if create_response['status'] == 'error':
                return create_response
        
        params = {
            "project_name": proj_name,
            "sample_id": sdid,
            "capture_id": capture,
            "sample_status": status,
            "pipeline_type": 1 if pipeline == "tumor_only" else 0
        }

        urllib3.disable_warnings(InsecureRequestWarning)
        response = requests.post(CREATE_SAMPLE_INFO, params=params, verify=False)
        response.raise_for_status()

        if response.status_code == 200:
            return {"status": "success", "message": f"Sample info for '{sdid}' updated successfully."}
        else:
            return {"status": "error", "message": "Failed to update sample info."}
    
    except requests.exceptions.HTTPError as http_err:
        return {"status": "error", "message": f"HTTP error occurred: {http_err}"}
    except requests.exceptions.ConnectionError:
        return {"status": "error", "message": "Connection error. Check your network or API endpoint."}
    except requests.exceptions.Timeout:
        return {"status": "error", "message": "Request timed out."}
    except requests.exceptions.RequestException as e:
        return {"status": "error", "message": f"An error occurred: {e}"}
    except ValueError:
        return {"status": "error", "message": "Failed to decode JSON from response."}


def extract_bam(sample, libdir, umi = False):
    """
    for re-run pipeline, need to extract sample specific bam files
    from library directory

    params: 
    sample: clinseqbarcode
    libdir: library directory path
    umi   : umi based processing or not
    """
    project = sample.split("-")[0]

    if project == "AL":
        sample_capture_str = compose_lib_capture_str(extract_unique_capture(sample))
        pattern_bam = libdir + "/" + sample_capture_str + "*nodups.bam"
    
    if project in ["PB", "LB", "PSFF", "iPCM"]:
        pattern_bam = libdir + "/" + sample + "*nodups.bam"
        pattern_umibam = libdir + "/" + sample + "*clipoverlap.bam"
        pattern_umibai = libdir + "/" + sample + "*clipoverlap.bai"
    
    pattern_bai = pattern_bam +  ".bai"
    if umi:
        umi_bam  = glob.glob(pattern_umibam)
        umi_bai  = glob.glob(pattern_umibai)
        if len(umi_bam) == 1:
            return(umi_bam[0], umi_bai[0])
        else:
            raise ValueError("Invalid UMI bam search: " + sample)

    nodups_bam = glob.glob(pattern_bam)
    nodups_bai = glob.glob(pattern_bai)
    if len(nodups_bam) == 1:
        return(nodups_bam[0], nodups_bai[0])
    else:
        raise ValueError(" ".join(["Invalid bam : ", sample, nodups_bam]))


def get_fastqs(wildcards, libdir):
    """
    helper function to get all fq files 
    """
    r1, r2 = find_fastqs(wildcards.sample, libdir)
    return r1 + r2


def get_containers(_path):
    """
    fetch containers for each conda env 
    
    """
    containers = {
        "base": os.path.join(_path, "autoseq-base.sif"),
        "franken": os.path.join(_path, "autoseq-franken.sif"),
        "gatk3": os.path.join(_path, "autoseq-gatk3.sif"),
        "gridss": os.path.join(_path, "autoseq-gridss.sif"),
        "jumble": os.path.join(_path, "autoseq-jumble.sif"),
        "purecn": os.path.join(_path, "autoseq-purecn.sif"),
        "ensemblvep": os.path.join(_path, "autoseq-ensemblvep-v113.1.sif"),
        "somaticseq": os.path.join(_path, "autoseq-somaticseq.sif"),
        "svcaller": os.path.join(_path, "autoseq-svcaller.sif"),
        "mulled-v2": os.path.join(_path, "autoseq-mulled-v2.sif"),
        "hmftools-markdups": os.path.join(_path, "autoseq-hmftools-markdups.sif"),
        "hmftools-redux": os.path.join(_path, "autoseq-hmftools-redux_1.0.sif"),
        "autoseq-rnastar": os.path.join(_path, "autoseq-rnastar-2.7.3a.sif"),
        "igv": os.path.join(_path, "autoseq-igvbatch.sif"),
        "dpyd": os.path.join(_path, "autoseq-typeDPYD.sif")
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
                profile, profile_path, jobs, jobdb, smk_option, use_singularity, bind_paths,
                cores='4', account=None, qos=None, max_run_hours=168):
        self.snakefile = snakefile
        self.cores = cores
        self.configfile = config
        self.cluster_config = cluster_config
        self.sdid = sdid
        self.project_id = project_id
        self.workdir = workdir
        self.profile = profile
        self.profile_path = profile_path
        self.jobs = jobs
        self.jobdb = jobdb
        self.dryrun = dryrun
        self.smk_option = smk_option
        self.use_singularity = use_singularity
        self.bind_paths = bind_paths
        self.account = account
        self.qos = qos
        self.max_run_hours = max_run_hours

    def build_cmd(self):
        dryrun = ''
        cores_cmd = ''
        jobs_cmd = ''
        smk_opt = ''
        singularity_cmd = ''
        slurm_options = ''

        if self.dryrun:
            dryrun = "-n"

        if self.smk_option:
            smk_opt = self.smk_option

        if self.profile_path:
            jobs_cmd = f" --jobs {self.jobs} "
            slurm_options = (f" --slurm-jobname-prefix {self.project_id}-{self.sdid} "
                             f" --profile {self.profile_path} ")
        else:
            # Local execution.
            cores_cmd = f" --cores {self.cores} "

        if self.use_singularity:
            singularity_cmd = " --software-deployment-method apptainer "
            singularity_cmd += " --apptainer-args ' "

            for path in self.bind_paths:
                singularity_cmd += f" --bind {path}:{path} "

            singularity_cmd += "'"


        cmd = ("snakemake --notemp -p --rerun-triggers mtime "
               f"--snakefile {self.snakefile} "
               f" --directory {self.workdir} "
               f" --configfile {self.configfile} "
               f" {dryrun} {cores_cmd} "
               f" {jobs_cmd} {singularity_cmd} "
               f" {slurm_options} "
               f" {smk_opt}")
        
        return cmd

    def _build_sbatch_head_job(self):
        """
        Build the sbatch "head job" script that runs the snakemake orchestrator
        on a compute node (instead of the login node). Snakemake itself keeps
        submitting the per-rule child jobs to SLURM via the workflow profile.

        return: sbatch script as a string
        """
        lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name=autoseq.head_job.{self.sdid}.%j",
            f"#SBATCH --output={self.workdir}/autoseq.head_job.{self.sdid}.%j.out",
            f"#SBATCH --error={self.workdir}/autoseq.head_job.{self.sdid}.%j.err",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=500M",
            f"#SBATCH --time={self.max_run_hours}:00:00",
            "#SBATCH --cpus-per-task=1",
        ]

        if self.account:
            lines.append(f"#SBATCH --account={self.account}")
        if self.qos:
            lines.append(f"#SBATCH --qos={self.qos}")

        lines.append("")
        lines.append("set -eo pipefail")
        lines.append("")
        lines.append(self.build_cmd())
        lines.append("")

        return "\n".join(lines)

    def submit_job(self):
        """
        Write the head-job sbatch script into the workdir and submit it to SLURM.

        return: the submitted SLURM job id
        """
        script_path = os.path.join(self.workdir, f"autoseq_smk_head_job_{self.sdid}.sh")

        with open(script_path, 'w') as sf:
            sf.write(self._build_sbatch_head_job())

        result = subprocess.run(
            ["sbatch", "--parsable", script_path],
            check=True, capture_output=True, text=True)

        return result.stdout.strip()


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
    unique_capture = extract_unique_capture(wildcards.sample, validation=False)
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
    unique_capture = extract_unique_capture(wildcards.sample, validation = False)
    targets = get_capture_name(unique_capture.capture_kit_id)

    if unique_capture.capture_kit_id in ["P2", "S2", "B2"]:
        return reference['small-design'][targets][key]
    
    return reference['targets'][targets][key]


def get_target_name(wildcards):
    """
    return capture id
    """
    unique_capture = extract_unique_capture(wildcards.sample, validation = False)
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
                            "S3": "probio_biomarkersignature2", # S3 pointed to S2 capture files
                            "P3": "probio_biomarkersignature2",
                            "N3": "pancancer3",
                            "UG": "ghent_GMCK"
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
 
