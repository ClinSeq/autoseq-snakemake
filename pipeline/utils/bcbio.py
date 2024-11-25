"""
Functions taken from bcbio. These are methods that can be run on the command line like so:

python -c 'import sys; from autoseq.util.bcbio import depth_freq_filter; depth_freq_filter(sys.stdin.read(), 0, \"bwa\")'

and

python -c 'import sys; from autoseq.util.bcbio import call_somatic; print call_somatic(sys.stdin.read())'

"""


def depth_freq_filter_input_stream(input_stream, tumor_index, aligner):
    processed_output = ""
    for line in input_stream.readlines():
        processed_output = processed_output + depth_freq_filter(line, tumor_index, aligner)
    return processed_output


def depth_freq_filter(line, tumor_index, aligner):
    """Command line to filter VarDict calls based on depth, frequency and quality.

    Looks at regions with low depth for allele frequency (AF * DP < 6, the equivalent
    of < 13bp for heterogygote calls, but generalized. Within these calls filters if a
    calls has:

    - Low mapping quality and multiple mismatches in a read (NM)
        For bwa only: MQ < 55.0 and NM > 1.0 or MQ < 60.0 and NM > 2.0
    - Low depth (DP < 10)
    - Low QUAL (QUAL < 45)

    Also filters in low allele frequency regions with poor quality, if all of these are
    true:
    - Allele frequency < 0.2
    - Quality < 55
    - P-value (SSF) > 0.06
    """
    if line.startswith("#CHROM"):
        headers = [('##FILTER=<ID=LowAlleleDepth,Description="Low depth per allele frequency '
                    'along with poor depth, quality, mapping quality and read mismatches.">'),
                   ('##FILTER=<ID=LowFreqQuality,Description="Low frequency read with '
                    'poor quality and p-value (SSF).">')]
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        sample_ft = {a: v for (a, v) in zip(parts[8].split(":"), parts[9 + tumor_index].split(":"))}
        qual = _safe_to_float(parts[5])
        dp = _safe_to_float(sample_ft.get("DP"))
        af = _safe_to_float(sample_ft.get("AF"))
        nm = _safe_to_float(sample_ft.get("NM"))
        mq = _safe_to_float(sample_ft.get("MQ"))
        ssfs = [x for x in parts[7].split(";") if x.startswith("SSF=")]
        pval = _safe_to_float(ssfs[0].split("=")[-1] if ssfs else None)
        fname = None
        if dp is not None and af is not None:
            if dp * af < 6:
                if aligner == "bwa" and nm is not None and mq is not None:
                    if (mq < 55.0 and nm > 1.0) or (mq < 60.0 and nm > 2.0):
                        fname = "LowAlleleDepth"
                if dp < 10:
                    fname = "LowAlleleDepth"
                if qual is not None and qual < 45:
                    fname = "LowAlleleDepth"
        if af is not None and qual is not None and pval is not None:
            if af < 0.2 and qual < 55 and pval > 0.06:
                fname = "LowFreqQuality"
        if fname:
            if parts[6] in set([".", "PASS"]):
                parts[6] = fname
            else:
                parts[6] += ";%s" % fname
        line = "\t".join(parts)
        return line


def call_somatic(line):
    """Call SOMATIC variants from tumor/normal calls, adding REJECT filters and SOMATIC flag.

    Assumes tumor/normal called with tumor first and normal second, as done in bcbio
    implementation.

    Uses MuTect like somatic filter based on implementation in speedseq:
    https://github.com/cc2qe/speedseq/blob/e6729aa2589eca4e3a946f398c1a2bdc15a7300d/bin/speedseq#L62

    Extracts the genotype likelihoods (GLs) from FreeBayes, which are like phred scores
    except not multiplied by 10.0 (https://en.wikipedia.org/wiki/Phred_quality_score).
    For tumors, we retrieve the best likelihood to not be reference (the first GL) and
    for normal, the best likelhood to be reference.

    After calculating the likelihoods, we compare these to thresholds to pass variants
    at tuned sensitivity/precision. Tuning done on DREAM synthetic 3 dataset evaluations.

    We also check that the frequency of the tumor exceeds the frequency of the normal by
    a threshold to avoid calls that are low frequency in both tumor and normal. This supports
    both FreeBayes and VarDict output frequencies.
    """
    # Thresholds are like phred scores, so 3.5 = phred35
    tumor_thresh, normal_thresh = 3.5, 3.5
    if line.startswith("#CHROM"):
        headers = ['##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic event">',
                   ('##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency '
                    'or phred likelihoods: tumor: %s, normal %s.">')
                   % (int(tumor_thresh * 10), int(normal_thresh * 10))]
        return "\n".join(headers) + "\n" + line
    elif line.startswith("#"):
        return line
    else:
        parts = line.split("\t")
        if _check_lods(parts, tumor_thresh, normal_thresh) and _check_freqs(parts):
            parts[7] += ";SOMATIC"
        else:
            if parts[6] in set([".", "PASS"]):
                parts[6] = "REJECT"
            else:
                parts[6] += ";REJECT"
        line = "\t".join(parts)
        return line


def _safe_to_float(x):
    if x is None:
        return None
    else:
        try:
            return float(x)
        except ValueError:
            return None


def _check_lods(parts, tumor_thresh, normal_thresh):
    """Ensure likelihoods for tumor and normal pass thresholds.

    Skipped if no FreeBayes GL annotations available.
    """
    try:
        gl_index = parts[8].split(":").index("GL")
    except ValueError:
        return True
    try:
        tumor_gls = [float(x) for x in parts[9].split(":")[gl_index].split(",") if x != "."]
        if tumor_gls:
            tumor_lod = max(tumor_gls[i] - tumor_gls[0] for i in range(1, len(tumor_gls)))
        else:
            tumor_lod = -1.0
    # No GL information, no tumor call (so fail it)
    except IndexError:
        tumor_lod = -1.0
    try:
        normal_gls = [float(x) for x in parts[10].split(":")[gl_index].split(",") if x != "."]
        if normal_gls:
            normal_lod = min(normal_gls[0] - normal_gls[i] for i in range(1, len(normal_gls)))
        else:
            normal_lod = normal_thresh
    # No GL inofmration, no normal call (so pass it)
    except IndexError:
        normal_lod = normal_thresh
    return normal_lod >= normal_thresh and tumor_lod >= tumor_thresh


def _check_freqs(parts):
    """Ensure frequency of tumor to normal passes a reasonable threshold.

    Avoids calling low frequency tumors also present at low frequency in normals,
    which indicates a contamination or persistent error.
    """
    thresh_ratio = 2.7
    try:  # FreeBayes
        ao_index = parts[8].split(":").index("AO")
        ro_index = parts[8].split(":").index("RO")
    except ValueError:
        ao_index, ro_index = None, None
    try:  # VarDict
        af_index = parts[8].split(":").index("AF")
    except ValueError:
        af_index = None
    if af_index is None and ao_index is None:
        raise NotImplementedError("Unexpected format annotations: %s" % parts[0])

    def _calc_freq(item):
        try:
            if ao_index is not None and ro_index is not None:
                ao = sum([int(x) for x in item.split(":")[ao_index].split(",")])
                ro = int(item.split(":")[ro_index])
                freq = ao / float(ao + ro)
            elif af_index is not None:
                freq = float(item.split(":")[af_index])
        except (IndexError, ValueError, ZeroDivisionError):
            freq = 0.0
        return freq

    tumor_freq, normal_freq = _calc_freq(parts[9]), _calc_freq(parts[10])
    return normal_freq <= 0.001 or normal_freq <= tumor_freq / thresh_ratio
