import gzip


def _is_gzipped(path):
    """
    check if a file is gzipped. query the "magic numbers" and see if 
    they match those of gzip. 
    See https://stackoverflow.com/a/47080739/5666087 
    """
    with open(path, "rb") as f:
        return f.read(2) == b'\x1f\x8b'


def load_fasta(path):
    """
    Loads a zipped or unzipped fasta file.
    Returns a dict file, with the headers as keys and another dict as 
    value. The nested dict has the string "sequence" as key and the 
    actual DNA sequence as value.
    Removes the first character of the header (should be >)
    """
    templates = {}
    f = gzip.open(path, 'rb') if _is_gzipped(path) else open(path, 'rb')
    lines = []
    for line in f:
        lines.append(line.rstrip())
        if len(lines) == 2:
            templates[lines[0][1:]] = {'sequence': lines[1]}
            lines = []
    f.close()
    return templates


def _process_fastq_entry(lines):
    """ Returns a dict of the fastq entry """
    keys = ['name', 'sequence', 'optional', 'quality']
    read = {key: value for key, value in zip(keys, lines)}
    return read


def load_fastq(path):
    """
    Loads a zipped or unzipped fastq file.
    Returns a list of dicts. Each dict the following keys, one for each 
    line of a fastq entry: name, sequence, optional, quality. 
    """
    reads = []
    f = gzip.open(path, 'rb') if _is_gzipped(path) else open(path, 'rb')
    lines = []
    for line in f:
        lines.append(line.rstrip())
        if len(lines) == 4:
            # Create a dict of the fastq entry
            reads.append(_process_fastq_entry(lines))
            lines = []
    f.close()
    return reads


def _clean_up_fastq_header(header, seperator):
    """
    Cleans up the header by removes additions that is added to the
    header of original template by the trimming programs:
    - Removes the "@M_"/"@" at the start of the header line.
    - Removes characters from the end of the line up to the seperator.
    """
     # remove the addition at the beginning of the header
    if header.startswith(b'@M_'):
        # AdapterRemoval, ClipAndMerge
        clean_header = header.split(b'@M_', 1)[1]
    elif header.startswith(b'@'):
        # leeHom, seqtk/adna, bbmerge, fastp, SeqPrep
        clean_header = header.split(b'@', 1)[1]
    # remove the addition at the end of the header
    clean_header = clean_header.rsplit(seperator, 1)[0]
    return clean_header


def clean_merged_reads(reads, templates, seperator):
    """
    - Removes unmerged reads (header starting with @F_ or @R_.) 
    - Cleans up the header by removes additions that is added to the
    header of original template by the trimming programs, such as 
    "@M_"/"@" at the start and something at the end of the header.
    - Removes reads that can not be assigned to a template 
    """
    merged_reads = []
    for read in reads:
        # discard unmerged reads
        if not read["name"].startswith((b'@F_', b'@R_')):
            # Clean up header
            read["name"] = _clean_up_fastq_header(read["name"], seperator)
            # discard reads that can not be assigned to a template
            if templates.get(read["name"]) is not None:
                merged_reads.append(read)
    return merged_reads