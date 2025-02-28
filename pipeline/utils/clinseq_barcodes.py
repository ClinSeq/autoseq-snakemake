import collections, os, re
import errno
import sys
from loguru import logger as Log


# Fields defining a unique library capture item. Note: A library capture item can
# include an item where no capture was performed; i.e. capture_kit_id indicates
# whole-genome sequencing:
UniqueCapture = collections.namedtuple(
    'UniqueCapture',
    ['project',
    'sdid',
    'sample_type',
    'sample_id',
    'library_kit_id',
    'capture_kit_id']
)


def normpath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def find_fastqs(library, libdir):
    """Find fastq files for a given library id in a given direcory.

        Returns a tuple with two lists:
    (['foo_1.fastq.gz', 'bar_1.fastq.gz'], # read 1
     ['foo_2.fastq.gz', 'bar_2.fastq.gz'])

    Supports the following file naming convenstions:
    *_1.fastq.gz / *_2.fastq.gz
    *_1.fq.gz / *_2.fq.gz
    *R1_nnn.fastq.gz / *R2_nnn.fastq.gz

    :rtype: tuple[str,str]
    """
    if not library:
        return (None, None)
    regex_fq1 = r'(.+)(_1.fastq.gz|_1.fq.gz|R1_\d{3}.fastq.gz)'
    regex_fq2 = r'(.+)(_2.fastq.gz|_2.fq.gz|R2_\d{3}.fastq.gz)'

    d = normpath(os.path.join(libdir, library))

    fq1s = []
    fq2s = []

    for f in os.listdir(d):
        match1 = re.search(regex_fq1, f)
        if match1:
            fn = "".join(match1.groups())
            fq1s.append(os.path.join(libdir, library, fn))
        match2 = re.search(regex_fq2, f)
        if match2:
            fn = "".join(match2.groups())
            fq2s.append(os.path.join(libdir, library, fn))

    fq1s.sort()
    fq2s.sort()

    return fq1s, fq2s


def check_sampledata(libdir, sampledata):
    """
    Check that data is available for all clinseq barcodes given in the sample json

    :param libdir, sample json
    :return sample json, clinseq barcodes
    """
    clinseq_barcodes = []
    for sample_type in ['N', 'T', 'CFDNA']:
        clinseq_barcodes_with_data = []
        for clinseq_barcode in sampledata[sample_type]:
            if data_available_for_clinseq_barcode(libdir, clinseq_barcode):
                clinseq_barcodes_with_data.append(clinseq_barcode)
            else:
                Log.error("No fastq files found for {} in dir {}".format(clinseq_barcode, libdir))
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), clinseq_barcode)
        
        sampledata[sample_type] =  clinseq_barcodes_with_data
        clinseq_barcodes.extend(clinseq_barcodes_with_data)
    
    return sampledata, clinseq_barcodes


def data_available_for_clinseq_barcode(libdir, clinseq_barcode):
    """
    Check that data is available for the specified clinseq barcode in the specified library folder.

    :param libdir: Directory name where fastqs are organised.
    :param clinseq_barcode: A valid clinseq barcode string
    :return: True if data is available, False otherwise
    """

    if not clinseq_barcode_is_valid(clinseq_barcode):
        raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

    filedir = os.path.join(libdir, clinseq_barcode)
    if not os.path.exists(filedir):
        Log.error("Dir {} does not exists for {}. Not using library.".format(filedir, clinseq_barcode))
        return False
    if find_fastqs(clinseq_barcode, libdir) == (None, None):
        Log.error("No fastq files found for {} in dir {}".format(clinseq_barcode, filedir))
        return False
    
    return True


def convert_to_libdict(clinseq_barcode):
    """
    Converts a clinseq_barcode into a dictionary - introduced to maintain compatibility
    with the mongoDB "libraries" collection, when parse_orderform is used in autoseq-api.

    :param clinseq_barcode: A clinseq barcode string 
    :return: A dictionary containing the field types and values present in the barcode
    """

    library_id = clinseq_barcode

    if not clinseq_barcode_is_valid(clinseq_barcode):
        raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

    capture_id = parse_capture_id(clinseq_barcode)
    type = parse_sample_type(clinseq_barcode)
    sample_id = parse_sample_id(clinseq_barcode)
    project_id = parse_project(clinseq_barcode)
    sdid = parse_sdid(clinseq_barcode)
    prep_id = parse_prep_id(clinseq_barcode)

    return {"library_id": library_id, "capture_id": capture_id, "type": type, "sample_id": sample_id,
            "project_id": project_id, "sdid": sdid, "prep_id": prep_id}


def extract_unique_capture(clinseq_barcode, validation = True):
    """
    Parses the specified clinseq barcode and produces a corresponding
    UniqueCapture representing the unique library capture corresponding
    to this clinseq barcode.

    :param clinseq_barcode_tuple: A clinseq barcode string
    :return: UniqueCapture named tuple
    """
    if validation:
        if not clinseq_barcode_is_valid(clinseq_barcode):
            raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)

    project = parse_project(clinseq_barcode)
    sdid = parse_sdid(clinseq_barcode)
    sample_type = parse_sample_type(clinseq_barcode)
    sample_id = parse_sample_id(clinseq_barcode)
    prep_id = parse_prep_id(clinseq_barcode)
    capture_id = parse_capture_id(clinseq_barcode)

    return UniqueCapture(project,
                         sdid,
                         sample_type,
                         sample_id,
                         extract_kit_id(prep_id, validation),
                         extract_kit_id(capture_id, validation))


def parse_project(clinseq_barcode):
    """
    Extract the project string from the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The project field from the input string.
    """

    return clinseq_barcode.split("-")[0]


def compose_sample_str(capture):
    """
    Produce a string representing the unique sample for the specified library capture item.

    :param capture: Named tuple indicating a unique sample library capture. 
    :return: Dash-delimited string of the fields uniquely identifying the sample.
    """

    return "{}-{}-{}-{}".format(
        capture.project,
        capture.sdid,
        capture.sample_type,
        capture.sample_id
    )


def compose_lib_capture_str(capture):
    """
    Produce a string for a unique library capture item.

    :param capture: A named tuple identifying a unique sample library capture.
    :return: A dash-delimted string of the fields uniquely identifying the capture.
    """

    return "{}-{}-{}-{}-{}-{}".format(
        capture.project,
        capture.sdid,
        capture.sample_type,
        capture.sample_id,
        capture.library_kit_id,
        capture.capture_kit_id)


def parse_sample_type(clinseq_barcode):
    """
    Extract the sample type from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The sample type field from the input string.
    """

    return clinseq_barcode.split("-")[3]


def parse_sample_id(clinseq_barcode):
    """
    Extract the sample ID from the clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string.
    :return: The sample ID field from the input string.
    """

    return clinseq_barcode.split("-")[4]


def parse_sdid(clinseq_barcode):
    """
    Extract the SDID from the clinseq barcode, including the "P-" prefix.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: The SDID field from the input string.
    """

    return "-".join(clinseq_barcode.split("-")[1:3])


def parse_prep_id(clinseq_barcode):
    """
    Extract the library prep ID (the entire string - not just the kit ID) from
    the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: Library prep ID string.
    """

    return clinseq_barcode.split("-")[5]


def parse_capture_id(clinseq_barcode):
    """
    Extract the library capture ID (the entire string - not just the capture kit ID)
    from the specified clinseq barcode.

    :param clinseq_barcode: Dash-delimited clinseq barcode string. 
    :return: Library capture ID string.
    """

    return clinseq_barcode.split("-")[6]


def extract_kit_id(kit_string, validation = True):
    """
    Extract the kit type from a specified kit string (either library or capture kit).

    :param kit_string: A string indicating a kit type, where the first two letters indicate the kit ID. 
    :return: The kit ID, comprising the first two letters.
    """
    if validation:
        if len(kit_string) < 3:
            raise ValueError("Invalid kit string: " + kit_string)

    return kit_string[:2]


def project_valid(project_str):
    return project_str in ["AL", "LB", "OT", "PB", "PSFF", "UL", "iPCM", "CRCR", "SARC", "CPC", "BM", "KA", "UM", "COV"]


def sdid_valid(sdid_str):
    return re.match("^P-[a-zA-Z0-9]+$", sdid_str) is not None


def sample_type_valid(sample_type_str):
    return sample_type_str in ["N", "T", "CFDNA"]


def sample_id_valid(sample_id_str):
    return re.match("^[a-zA-Z0-9]+$", sample_id_str) is not None


def prep_id_valid(prep_id_str):
    return re.match("^[A-Z]{2}[0-9_]+$", prep_id_str) is not None


def capture_id_valid(capture_id_str):
    return (re.match("^[A-Z0-9]{2}[0-9_]+$", capture_id_str) is not None) or \
           (capture_id_str == "WGS")


def clinseq_barcode_is_valid(clinseq_barcode):
    """
    Test the structure of the specified clinseq barcode for validity.
    Barcode format is defined at https://github.com/clinseq/autoseq

    :param clinseq_barcode: The input clinseq barcode.
    :return: True if the barcode has valid structure, False otherwise.
    """

    fields = clinseq_barcode.split("-")
    if len(fields) != 7:
        return False

    project_str = fields[0]
    sdid_str = "-".join(fields[1:3])
    sample_type_str = fields[3]
    sample_id_str = fields[4]
    prep_id_str = fields[5]
    capture_id_str = fields[6]

    barcode_valid = \
        project_valid(project_str) and sdid_valid(sdid_str) and sample_type_valid(sample_type_str) and \
        sample_id_valid(sample_id_str) and prep_id_valid(prep_id_str) and capture_id_valid(capture_id_str)

    return barcode_valid


def extract_clinseq_barcodes(input_filename):
    """
    Extrat clinseq barcodes from the specified input file:

    :param input_filename: Either a .txt listing clinseq barcodes one per line,
    or a .xlsx order form file containing the barcodes.

    :return: A list of (not-yet validated) dash-delimited clinseq barcodes.
    """

    toks = input_filename.split(".")

    if toks[-1] == "txt":
        return list(set([line.strip() for line in open(input_filename).readlines()]))
    else:
        raise ValueError("Invalid clinseq barcodes file type: " + input_filename)


def validate_clinseq_barcodes(clinseq_barcodes):
    """Checks all the specified clinseq barcodes for validity and raises an
    Exception if any are not.
    
    :param clinseq_barcodes: List of clinseq barcode strings.
    """

    # Check the input clinseq barcodes for validity:
    for clinseq_barcode in clinseq_barcodes:
        if not clinseq_barcode_is_valid(clinseq_barcode):
            raise ValueError("Invalid clinseq barcode: " + clinseq_barcode)


def create_scaffold_sampledict(sdids):
    """
    Generate a scaffold sample dictionary, with SDIDs as keys and clinseq barcode information
    dictionaries as values.

    :param sdids: List of SDID strings
    :return: Dictionary with SDID keys and empty clinseq barcode information dictionaries as values.
    """

    scaffold_dict = {}
    for sdid in sdids:
        curr_clinseq_barcode_info = {"sdid": sdid, "N": [], "T": [], "CFDNA": []}
        scaffold_dict[sdid] = curr_clinseq_barcode_info

    return scaffold_dict


def populate_clinseq_barcode_info(clinseq_barcode_info, clinseq_barcode):
    """
    Populates the specified clinseq barcode information item with the information
    in the specified clinseq barcode.

    :param clinseq_barcode_info: A dictionary containing clinseq barcode information
    for a single SDID.
    :param clinseq_barcode: A validated clinseq barcode string.
    """

    if not clinseq_barcode_is_valid(clinseq_barcode):
        raise ValueError("Invalid clinseq barcode: ", clinseq_barcode)

    # Append this clinseq barcode to the relevant field in the specified clinseq barcode
    # info dictionary:
    sample_type = parse_sample_type(clinseq_barcode)
    clinseq_barcode_info[sample_type].append(clinseq_barcode)


def convert_barcodes_to_sampledict(clinseq_barcodes):
    """
    Coverts the specified clinseq barcode strings into a dictionary linking
    from SDID to clinseq barcode information for that individual.

    :param clinseq_barcodes: A list of valid clinseq barcode strings
    :return: A dictionary with the required structure.
    """

    for curr_clinseq_barcode in clinseq_barcodes:
        if not clinseq_barcode_is_valid(curr_clinseq_barcode):
            raise ValueError("Invalid clinseq barcode: ", curr_clinseq_barcode)

    # Extract set of unique SDIDs from the specified clinseq barcodes:
    sdids = set([parse_sdid(clinseq_barcode) for clinseq_barcode in clinseq_barcodes])

    # Create a scaffold dictionary from those SDIDs:
    sdid_to_clinseq_barcode_info = create_scaffold_sampledict(sdids)

    for clinseq_barcode in clinseq_barcodes:
        populate_clinseq_barcode_info(sdid_to_clinseq_barcode_info[parse_sdid(clinseq_barcode)],
                                      clinseq_barcode)

    return sdid_to_clinseq_barcode_info
