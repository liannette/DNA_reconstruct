#!/usr/bin/python3

""" 
This script was initially written by Leonardo and then improved by 
Annette Lien (a.lien@posteo.de) (12.09.2022)
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
        "-l", "--fraglen", action="store", type=int, required=True,
        help="fraglen") 
    parser.add_argument(
        "-n", "--nfrags", action="store", type=int, required=True,
        help="nfrags")  
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
        args.fraglen,
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


def main(template_path, readm_path, nfrags, fraglen, export_path, tool_name):

    # Load files --------------------------------------------------------------

    templates = common.load_fasta(template_path)

    # Check for duplicate fragments
    if len(templates) != nfrags:
        print(f"ATTENTION: number of total_sequences is {len(templates)}, " 
              f"but the nfrags is {nfrags}. Possible reason: duplicate "
              "fragments")

    reads = common.load_fastq(readm_path)
    # seperator: this character and all charaters to the right of it
    # will be removed from the fastq header
    seperator = b'-'
    reads = common.clean_merged_reads(reads, templates, seperator)

    # Analysis and Results ----------------------------------------------------

    # edit distances of all reads
    edit_dist_list = get_edit_distances(reads, templates)
    edit_dist_string = ""
    for edit_dist in sorted(set(edit_dist_list)):
        edit_dist_string += f"{edit_dist}:"
        edit_dist_string += f"{edit_dist_list.count(edit_dist)} "
    # Number of dropped reads
    dropped_reads_cnt = nfrags - len(reads)
    # NT change per NT (%)
    if len(reads) > 0:
        avg_divergence_per_nt = sum(edit_dist_list) / len(reads) / fraglen * 100
        avg_divergence_per_nt = round(avg_divergence_per_nt, 3)
    else:
        avg_divergence_per_nt = "NA"


    #################### export results ####################
            
    with open(export_path, 'w') as f:
        f.write(
            "program,"
            "filename,"
            "nfrags,"
            "fraglen,"
            "total_sequences,"
            "total_reads,"
            "dropped_reads,"
            "avg_divergence_per_nt,"
            "edit_distances"
            "\n")
        f.write(f"{tool_name},"
                f"{os.path.basename(readm_path)},"
                f"{nfrags},"
                f"{fraglen},"
                f"{len(templates)},"
                f"{len(reads)},"
                f"{dropped_reads_cnt},"
                f"{avg_divergence_per_nt},"
                f"{edit_dist_string.rstrip()}"
                )


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)