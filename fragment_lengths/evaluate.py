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
        "-l", "--fraglen", action="store", type=int, required=True,
        help="fraglen") 
    parser.add_argument(
        "-n", "--nfrags", action="store", type=int, required=True,
        help="nfrags")  

    # optional arguments, only needed if exporting the results
    parser.add_argument(
        "-o", "--out", action="store", type=str, required=False,
        dest="export_path", help="Path for the output csv file")
    parser.add_argument(
        "-t", "--tool", action="store", type=str, required=False,
        dest="tool_name", help="Name of the tool used for trimming")

    args = parser.parse_args()
    required_arguments = [
        args.templates_path, 
        args.readm_path,
        args.nfrags, 
        args.fraglen,
        ]
    optional_arguments = [args.export_path, args.tool_name]
    
    if set(optional_arguments) == {None}:
        # no optional arguments have been passed
        return required_arguments
    elif None not in optional_arguments:
        # all optional arguments have been passed
        return required_arguments + optional_arguments
    else:
        # only some optional arguments have been passed
        print("parseDNAfragments.py: error: only some optional arguments have "
              "been passed")
        parser.print_help()
        sys.exit(2)


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
    dropped_reads_cnt = len(templates) - len(merged_reads)
    incorrect_length_cnt = 0
    divergent_with_correct_length_cnt = 0
    perfectly_reconstructed_cnt = 0
    divergences = []

    for read in merged_reads:
        read_seq = read['sequence']
        template_seq = templates[read['name']]['sequence']

        if read_seq == template_seq:
            perfectly_reconstructed_cnt += 1
        else:
            # length change
            if len(template_seq) == len(read_seq):
                divergent_with_correct_length_cnt += 1
            else:
                incorrect_length_cnt += 1
            # edit distance
            divergences.append(_levenshtein_distance(template_seq, read_seq))

    results = [
        dropped_reads_cnt, 
        incorrect_length_cnt, 
        divergent_with_correct_length_cnt, 
        perfectly_reconstructed_cnt,
        divergences]     
    
    return results


def _percentage(cnt, total_cnt):
    return round(cnt/total_cnt*100, 3)


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


    results = analyze_merged_reads(reads, templates)
    # Dropped reads
    dropped_reads_cnt = results[0]
    # Incorrect length
    incorrect_length_cnt = results[1]
    # Right lenght but other sequence
    correct_len_incorrect_seq_cnt = results[2]
    # Perfectly reconstructed
    perfectly_reconstructed_cnt = results[3]
    # list of all divergences
    divergences = results[4]
    # Number of difergent reads
    divergent_cnt = len(divergences)
    assert divergent_cnt == (incorrect_length_cnt 
                               + correct_len_incorrect_seq_cnt)
    
    if len(reads) != 0:
        # Average number NT changes per merged read
        avg_divergence = sum(divergences) / len(reads)
    else:
        avg_divergence = 'NA'
    
    # Percentages
    dropped_reads_percent = _percentage(dropped_reads_cnt, nfrags)
    incorrect_length_percent = _percentage(incorrect_length_cnt, nfrags)
    correct_len_incorrect_seq_percent = _percentage(
        correct_len_incorrect_seq_cnt, nfrags)
    perfectly_reconstructed_percent = _percentage(perfectly_reconstructed_cnt, 
                                                  nfrags)
    if len(reads) != 0:
        # Percent of non-dropped reads that are not perfectly reconstucted
        divergent_reads_percent = _percentage(divergent_cnt, len(reads))
        # NT changes per NT (%)
        avg_divergence_percent = _percentage(avg_divergence, fraglen) 
    else:
        divergent_reads_percent = 'NA'
        avg_divergence_percent = 'NA'


    #################### Print or export results ####################

    if export_path is None:
        print(f"{os.path.basename(readm_path)}:")
        print(f"Dropped reads: {dropped_reads_cnt} of {nfrags} total "
              f"sequences ({dropped_reads_percent}%)")
        print(f"Incorrect length reads: {incorrect_length_cnt} of "
              f"{nfrags} total sequences ({incorrect_length_percent}%)")
        print(f"Perfectly reconstructed fragments (edit distance = 0 and "
              f"correct length): {perfectly_reconstructed_cnt} of "
              f"{nfrags} total sequences "
              f"({perfectly_reconstructed_percent}%)")
        # if len(reads) != 0:
        #     print(f"Divergent reads (edit distance > 0): {divergent_cnt} of "
        #           f"{len(reads)} merged (non-dropped) reads "
        #           f"({divergent_reads_percent}%)")
        #     print(f"Average divergence (edit distance): "
        #           f"{round(avg_divergence, 3)} of {fraglen} "
        #           f"nucleotides ({avg_divergence_percent}%)\n")

    else:
            
        with open(export_path, 'w') as f:
            f.write(
                "program,"
                "filename,"
                "nfrags,"
                "fraglen,"
                "total_sequences,"
                "total_reads,"
                "dropped_reads,"
                "incorrect_length,"
                "correct_len_incorrect_seq,"
                "perfectly_reconstructed,"
                "dropped_reads_percentage,"
                "incorrect_length_percentage,"
                "correct_len_incorrect_seq_percentage,"
                "perfectly_reconstructed_percentage,"
                "divergent_reads,"
                #"average_divergence,"
                #"divergent_reads_percentage,"
                #"average_divergence_percentage,"
                "\n")
            f.write(f"{tool_name},"
                    f"{os.path.basename(readm_path)},"
                    f"{nfrags},"
                    f"{fraglen},"
                    f"{len(templates)},"
                    f"{len(reads)},"
                    f"{dropped_reads_cnt},"
                    f"{incorrect_length_cnt},"
                    f"{correct_len_incorrect_seq_cnt},"
                    f"{perfectly_reconstructed_cnt},"
                    f"{dropped_reads_percent},"
                    f"{incorrect_length_percent},"
                    f"{correct_len_incorrect_seq_percent},"
                    f"{perfectly_reconstructed_percent},"
                    f"{divergent_cnt},"
                    #f"{avg_divergence},"
                    #f"{divergent_reads_percent},"
                    #f"{avg_divergence_percent}"
                    )


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)