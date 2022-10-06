import gzip


def _is_gzipped(path):
    """
    check if a file is gzipped. query the "magic numbers" and see if 
    they match those of gzip. 
    See https://stackoverflow.com/a/47080739/5666087 
    """
    with open(path, "rb") as f:
        return f.read(2) == b'\x1f\x8b'
    

def _clean_up_fastq_header(header, seperator):
    """
    Removes the "@M_"/"@" at the start of the header line.
    Removes characters from the end of the line up to the 
    seperator.
    """
    # make definition name uniform for dict lookup
    if header.startswith(b'@M_'):
        # AdapterRemoval, ClipAndMerge
        defname = header.split(b'@M_', 1)[1].rsplit(seperator, 1)[0]
    else:
        # leeHom, seqtk/adna, bbmerge, fastp, SeqPrep
        defname = header.split(b'@', 1)[1].rsplit(seperator, 1)[0]
    return defname


def _process_fastq_entry(lines, seperator):
    """ Returns a dict of the fastq entry """
    keys = ['name', 'sequence', 'optional', 'quality']
    reading = {k: v for k, v in zip(keys, lines)}
    # make definition name uniform for dict lookup
    reading['name'] = _clean_up_fastq_header(lines[0], seperator)
    return reading


def load_merged_fastq(path, templates, seperator):
    """
    Loads a zipped fastq file. Ignores unmerged reads.
    Returns a list of dicts. Each dict has the header of a sequence as 
    key. Header is without the "@M_"/"@" at the start and without the 
    "-" addition at the end of the header line
    """
    merged_reads = []
    f = gzip.open(path, 'rb') if _is_gzipped(path) else open(path, 'rb')
    lines = []
    for line in f:
        lines.append(line.rstrip())
        if len(lines) == 4:
            # skip unmerged reads
            if not lines[0].startswith((b'@F_', b'@R_')):
                # Create a dict of the fastq entry
                reading = _process_fastq_entry(lines, seperator)
                if templates.get(reading['name']) is not None:
                    merged_reads.append(reading)
            lines = []
    f.close()
    return merged_reads