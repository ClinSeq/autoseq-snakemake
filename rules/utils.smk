

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
    logger.debug("Looking for fastq files for library {library} in {libdir}".format(library=library, libdir=libdir))

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

    logging.debug("Found {}".format((fq1s, fq2s)))
    return fq1s, fq2s
