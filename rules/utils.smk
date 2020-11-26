

def get_fastq(all_clinseq_barcodes, libdir):
    fq1_files = []
    fq2_files = []
    for clinseq_barcode in all_clinseq_barcodes:
        fq1, fq2 = find_fastqs(clinseq_barcode, libdir)
        fq1_files.extend(fq1)
        fq2_files.extend(fq2)
    
    return fq1_files, fq2_files


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
    regex_fq1 = '(.+)(_1\.fastq.gz|_1\.fq.gz|R1_\d{3}.fastq.gz)'
    regex_fq2 = '(.+)(_2\.fastq.gz|_2\.fq.gz|R2_\d{3}.fastq.gz)'

    d = normpath(os.path.join(libdir, library))
    print("Looking for fastq files for library {library} in {libdir}".format(library=library, libdir=libdir))

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

    print("Found {}".format((fq1s, fq2s)))
    return fq1s, fq2s


def get_readgroup(wildcards):
    library_id = parse_prep_id(wildcards.sample)
    sample_string = compose_sample_str(extract_unique_capture(wildcards.sample))

    readgroup = "\"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:ILLUMINA\"".format(\
        rg_id=wildcards.sample, rg_sm=sample_string, rg_lb=library_id)
    
    return readgroup


def get_targets(wildcards, reference):
    """
    return capture id
    """
    unique_capture = extract_unique_capture(wildcards.sample)
    targets = get_capture_name(unique_capture.capture_kit_id)
    
    return reference['targets'][targets]['targets-bed-slopped20']


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
                            "PN": "pancancer2"
                            }

    if capture_kit_code == 'WG':
        return 'lowpass_wgs'

    else:
        return capture_kit_loopkup[capture_kit_code]


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
 