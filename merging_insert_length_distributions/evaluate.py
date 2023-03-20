#!/usr/bin/python3

""" 
This script was initially written by Leonardo and then improved by 
Annette Lien (a.lien@posteo.de) (12.09.2022)

------------------------------------------------------------------------


------------------------------------------------------------------------
Old comments from Leonardo: 

We report:
- Are there missing reads? Dropped due to low quality or something?
- how many times we got the length right?
- the divergence between original and predicted sequence?

If a program has a high rate of dropping reads, is might be neccessary
to adjust some parameters of the trimming algorithms.

It might happen that you have a collision i.e. the same sequence got
sampled twice. It is unlikely but I would not worry about it.

A potential pitfall is that it is possible that certain programs might
not output the same number of sequences that they take in. If that is 
the case will have to modify the commands in order to output all of the 
sequences. 
"""


import sys
import os
import edlib
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import common 



def parse_arguments():
    """
    """
    parser = argparse.ArgumentParser(
        description="Performs analysis of the reconstructed reads. "
                    "If the results should be exported, all the optional "
                    "arguments must be supplied.")
    
    # required arguments
    parser.add_argument(
        "-in1", "--templates", action="store", type=str, required=True, 
        dest="templates_path", help='gzipped or unzipped fasta file of '
                                    'the simulated DNA templates. ')
    parser.add_argument(
        "-in2", "--mreads", action="store", type=str,  required=True, 
        dest="readm_path", help='gzipped or unzipped fastq file of the '
                                'trimmed and merged reads')
    parser.add_argument(
        "-d", "--fraglendist", action="store", type=str, required=True,
        help="fragment length distribution name") 
    parser.add_argument(
        "-n", "--nfrags", action="store", type=int, required=True,
        help="nfrags")  
    parser.add_argument(
        "-qs", "--qualityshift", action="store", type=int, required=True,
        help="the quality shift used when simulating the reads (when shifting "
             "scores by x, the error rate is 1/(10^(x/10)) of the default "
             "profile.") 
    parser.add_argument(
        "-o", "--out", action="store", type=str, required=True,
        dest="export_path", help="Path for the output csv file")
    parser.add_argument(
        "-t", "--tool", action="store", type=str, required=True,
        dest="tool_name", help="Name of the tool used for trimming")

    args = parser.parse_args()
    arguments = [
        args.templates_path, 
        args.readm_path,
        args.nfrags, 
        args.fraglendist,
        args.qualityshift,
        args.export_path, 
        args.tool_name
        ]
    
    return arguments


def _levenshtein_distance(template_seq, read_seq):
    """
    count how many nucleotides are different, use edlib's levenshtein
    algo for edit distance
    """
    return edlib.align(template_seq, read_seq)['editDistance']


def get_edit_distances(merged_reads, templates):
    edit_distances = []
    for read in merged_reads:
        read_seq = read['sequence']
        template_seq = templates[read['name']]['sequence']
        edit_distances.append(_levenshtein_distance(template_seq, read_seq))  
    return edit_distances


def make_edit_distances_string(edit_dist_list):
    """Counts the occurence of each edit distance and returns a 
    dict-like string""" 
    edit_dist_string = ""
    for edit_dist in sorted(set(edit_dist_list)):
        edit_dist_string += f"{edit_dist}:"
        edit_dist_string += f"{edit_dist_list.count(edit_dist)} "
    return edit_dist_string


def main(template_path, readm_path, nfrags, distname, qs, export_path, 
         tool_name):

    # Load files --------------------------------------------------------------

    templates = common.load_fasta(template_path)

    # Check for duplicate fragments
    # if len(templates) != nfrags:
    #     print(f"ATTENTION: number of total_sequences is {len(templates)}, " 
    #           f"but the nfrags is {nfrags}. Possible reason: duplicate "
    #           "fragments")

    reads = common.load_fastq(readm_path)
    # seperator: this character and all charaters to the right of it
    # will be removed from the fastq header
    seperator = b'-'
    reads = common.clean_merged_reads(reads, templates, seperator)

    # Analysis and Results ----------------------------------------------------

    # edit distances of all reads
    edit_dist_list = get_edit_distances(reads, templates)
    edit_dist_string = make_edit_distances_string(edit_dist_list)
    # Number of dropped reads
    # Comment: why take the length of templates and not nfrags?
    dropped_reads_cnt = nfrags - len(reads)


    #################### export results ####################
            
    with open(export_path, 'w') as f:
        f.write(
            "program,"
            "filename,"
            "nfrags,"
            "fraglen_distribution,"
            "quality_shift,"
            "total_sequences,"
            "total_reads,"
            "dropped_reads,"
            "edit_distances"
            "\n")
        f.write(f"{tool_name},"
                f"{os.path.basename(readm_path)},"
                f"{nfrags},"
                f"{distname},"
                f"{qs},"
                f"{len(templates)},"
                f"{len(reads)},"
                f"{dropped_reads_cnt},"
                f"{edit_dist_string.rstrip()}"
                )


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)