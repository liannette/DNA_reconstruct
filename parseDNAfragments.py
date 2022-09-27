#!/usr/bin/python3

""" 
This script was initially written by Leonardo and then rewritten by 
Annette Lien (a.lien@posteo.de) (12.09.2022)

usage: parseDNAfragments.py [-h] [-o EXPORT_PATH] [-nf NFRAGS] [-fl FRAGLEN]
                            [-t PROGRAM_NAME]
                            templates_path readm_path

Performs analysis of the reconstructed reads.If the results should be
exported, all the optional arguments must be supplied.

positional arguments:
  templates_path        fasta file of the simulated DNA templates
  readm_path            fastq file of the trimmed and merged reads

optional arguments:
  -h, --help            show this help message and exit
  -o EXPORT_PATH, --out EXPORT_PATH
                        Path for the output csv file
  -nf NFRAGS, --nfrags NFRAGS
                        nfrags
  -fl FRAGLEN, --fraglen FRAGLEN
                        fraglen
  -t PROGRAM_NAME, --tool PROGRAM_NAME
                        Name of the program used for trimming

-----------------------------------------------------------------------
leonardos notes

We report:
- Are there missing reads? Dropped due to low quality or something?
- how many times we got the length right?
- the divergence between original and predicted sequence?

Probably have to adjust some parameters for the trimming algorithms if
a program has a high rate of dropping reads.

Now the difficulty is that the order of b the reads might not be the
same. What I would suggest is a dictionary using the definition lines
where the value would be the original sequence. The other difficulty
is that certain programs modify the def line by adding for example M_
for merged reads.

ok, also, it might happen that you have a collision i.e. we sampled
the same sequence twice. It is unlikely but I would not worry about
it.

Ok great, a potential pitfall is that it is possible that certain
programs might not output the same number of sequences that they take
in. If that is the case will have to modify the commands in order to
output all of the sequences. 
"""


import sys, os
import gzip
import edlib
import argparse


def parse_arguments():
    """
    """
    parser = argparse.ArgumentParser(
        description="Performs analysis of the reconstructed reads."
                    "If the results should be exported, all the optional "
                    "arguments must be supplied.")
    # required arguments
    parser.add_argument(
        "templates_path", type=str, 
        help='fasta file of the simulated DNA templates')
    parser.add_argument(
        "readm_path", type=str, 
        help='fastq file of the trimmed and merged reads')
    # optional arguments, only needed if exporting the results
    parser.add_argument(
        "-o", "--out", action="store", type=str, required=False,
        dest="export_path", help="Path for the output csv file")
    parser.add_argument(
        "-nf", "--nfrags", action="store", type=int, required=False,
        help="nfrags")  
    parser.add_argument(
        "-fl", "--fraglen", action="store", type=int, required=False,
        help="fraglen") 
    parser.add_argument(
        "-t", "--tool", action="store", type=str, required=False,
        dest="program_name", help="Name of the program used for trimming")

    args = parser.parse_args()

    if len(sys.argv) == 3:
        return args.templates_path, args.readm_path
    elif len(sys.argv) == 7:
        return (args.templates_path, 
                args.readm_path, 
                args.export_path, 
                args.nfrags, 
                args.fraglen, 
                args.program_name
        )
    else:
        parser.print_help()
        sys.exit(2)


def load_initial(templates_path):
    """
    Loads a zipped fasta file.
    Returns a dict file, with the headers as keys and another dict as 
    value. The nested dict has the string "sequence" as key and the 
    actual DNA sequence as value.
    Removes the first character of the header (should be >)
    """
    templates = {}
    templates_cnt = 0
    with gzip.open(templates_path, 'rb') as f:
        lines = []
        for line in f:
            lines.append(line.rstrip())
            if len(lines) == 2:
                templates[lines[0][1:]] = {'sequence': lines[1]}
                templates_cnt += 1
                lines = []
    # Note from annette:
    # I think templates_cnt is not nesseccary, because 
    # templates_cnt should be equal len(templates)
    return templates, templates_cnt


def _clean_up_fastq_header(header):
    """
    Removes the "@M_"/"@" at the start of the header line and removes 
    the "-" addition at the end of the header line
    """
    # make definition name uniform for dict lookup
    if header.startswith(b'@M_'):
        # AdapterRemoval, ClipAndMerge
        defname = header.split(b'@M_',1)[1].rsplit(b'-',1)[0]
    else:
        # leeHom, seqtk/adna, bbmerge, fastp
        defname = header.split(b'@',1)[1].rsplit(b'-',1)[0]
    return defname


def _process_fastq_entry(lines):
    """ Returns a dict of the fastq entry """
    keys = ['name', 'sequence', 'optional', 'quality']
    reading = {k: v for k, v in zip(keys, lines)}
    # make definition name uniform for dict lookup
    reading['name'] = _clean_up_fastq_header(lines[0])
    return reading


def load_readm(readm_path, templates):
    """
    Loads a zipped fastq file. Ignores unmerged reads.
    Returns a list of dicts. Each dict has the header of a sequence as 
    key. Header is without the "@M_"/"@" at the start and without the 
    "-" addition at the end of the header line
    """
    merged_reads = []
    with gzip.open(readm_path, 'rb') as f:
        lines = []
        for line in f:
            lines.append(line.rstrip())
            if len(lines) == 4:
                # don't consider forward (@F_) or reverse (@R_) reads
                # reads processed with AdapterRemoval have a @M_
                #
                if not lines[0].startswith((b'@F_', b'@R_')):
                    # Create a dict of the fastq entry
                    reading = _process_fastq_entry(lines)
                    if templates.get(reading['name']) is not None:
                        merged_reads.append(reading)
                lines = []
    return merged_reads


def _levenshtein_distance(template_seq, read_seq):
    """
    count how many nucleotides are different, use edlib's levenshtein
    algo for edit distance
    """
    return edlib.align(template_seq, read_seq)['editDistance']


def analyze_merged_reads(merged_reads, templates):
    """
    Returns:
        list of
            number of reads that have been perfectly reconstructed reads
            list of divergences of the not perfectly reconstructed reads
            number of reads that have the incorrect length

    """
    perfectly_reconstructed_cnt = 0
    incorrect_length_cnt = 0
    divergences = []

    for read in merged_reads:
        read_seq = read['sequence']
        template_seq = templates[read['name']]['sequence']

        if read_seq == template_seq:
            perfectly_reconstructed_cnt += 1
        else:
            # edit distance
            divergences.append(_levenshtein_distance(template_seq, read_seq))
            # length change
            if len(template_seq) != len(read_seq):
                incorrect_length_cnt += 1

    return perfectly_reconstructed_cnt, incorrect_length_cnt, divergences


def _get_template_sequence_length(templates):
    """
    All sequences in templates have the same length
    """
    return len(next(iter(templates.values()))['sequence'])


def _in_percentage(cnt, total_cnt):
    return round(cnt/total_cnt*100, 3)


def main(args):

    #################### Load data ####################

    template_path = args[0]
    readm_path = args[1]
    templates, templates_cnt = load_initial(template_path)
    reads = load_readm(readm_path, templates)


    #################### Analysis and Results ####################

    # all templates are assumed to have the same length
    templates_length = _get_template_sequence_length(templates)

    # Dropped reads
    dropped_reads_cnt = templates_cnt - len(reads)
    dropped_reads_percent = _in_percentage(dropped_reads_cnt, templates_cnt)

    results = analyze_merged_reads(reads, templates)
    # Perfectly reconstructed
    perfectly_reconstructed_cnt = results[0]
    perfectly_reconstructed_percent = _in_percentage(
        perfectly_reconstructed_cnt, templates_cnt)
    # Incorrect length
    incorrect_length_cnt = results[1]
    incorrect_length_percent = _in_percentage(incorrect_length_cnt, 
                                               templates_cnt)
    # Divergence
    divergences = results[2]
    divergences_cnt = len(divergences)
    if len(reads) != 0:
        divergent_reads_percent = _in_percentage(divergences_cnt/len(reads))
        avg_divergence = sum(divergences)/len(reads) # NT changes per read
        avg_divergence_percent = _in_percentage(
            avg_divergence, templates_length) # NT changes per NT (%)
    else:
        avg_divergence = 'NA'
        divergent_reads_percent = 'NA'
        avg_divergence_percent = 'NA'


    #################### Print results ####################

    print(f"{os.path.basename(readm_path)}:")
    print(f"Dropped reads: {dropped_reads_cnt} of {templates_cnt} total "
          f"sequences ({dropped_reads_percent}%)")
    print(f"Incorrect length reads: {incorrect_length_cnt} of {templates_cnt} "
          f"total sequences ({incorrect_length_percent}%)")
    print(f"Perfectly reconstructed fragments (edit distance = 0 and correct "
          f"length): {perfectly_reconstructed_cnt} of {templates_cnt} total "
          f"sequences ({perfectly_reconstructed_percent}%)")
    if len(reads) != 0:
        print(f"Divergent reads (edit distance > 0): {divergences_cnt} of "
              f"{len(reads)}  merged (non-dropped) reads "
              f"({divergent_reads_percent}%)")
        print(f"Average divergence (edit distance): {round(avg_divergence, 3)}"
              f" of {templates_length} nucleotides ({avg_divergence_percent}%)\n")

    #################### Export results to CSV ####################

    if len(args) == 7:
        export_path = args[3]
        nfrags = args[4]
        fraglen = args[5]
        program_name = args[6]
        with open(export_path, 'w') as f:
            f.write(
                "program,filename,nfrags,fraglen,total_sequences,total_reads,"
                "dropped_reads,incorrect_length_reads,divergent_reads,"
                "average_divergence,dropped_reads_percentage,"
                "incorrect_length_percentage,divergent_reads_percentage,"
                "average_divergence_percentage,perfectly_reconstructed,"
                "perfectly_reconstructed_percentage\n")
            f.write(f"{program_name},"
                    f"{os.path.basename(readm_path)},"
                    f"{nfrags},"
                    f"{fraglen},"
                    f"{templates_cnt},"
                    f"{len(reads)},"
                    f"{dropped_reads_cnt},"
                    f"{incorrect_length_cnt},"
                    f"{divergences_cnt},"
                    f"{avg_divergence},"
                    f"{dropped_reads_percent},"
                    f"{incorrect_length_percent},"
                    f"{divergent_reads_percent},"
                    f"{avg_divergence_percent},"
                    f"{perfectly_reconstructed_cnt},"
                    f"{perfectly_reconstructed_percent}")


if __name__ == "__main__":

    args = parse_arguments()
    main(args)