#!/usr/bin/python3

import argparse
import sys
import gzip
import math
import pandas as pd
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import common
import numpy as np
import scipy.stats as st


# running this script for with files that have 100 Mio fragments, uses 
# an insane amount of RAM, up to 120 gb. I could consider rewriting this
# script to work with less RAM


def parse_arguments():
    """
    """
    parser = argparse.ArgumentParser(
        description="Performs analysis of the reconstructed reads. "
                    "If the results should be exported, all the optional "
                    "arguments must be supplied.")
    
    # required arguments
    parser.add_argument(
        "--templates", action="store", type=str, required=True, 
        dest="templates_path", help='gzipped or unzipped fasta file of '
                                    'the simulated DNA templates. ')
    parser.add_argument(
        "--s1", action="store", type=str,  required=True, 
        dest="s1_path", help='gzipped or unzipped fastq file of the simulated '
                             's1 reads')
    parser.add_argument(
        "--s2", action="store", type=str,  required=True, 
        dest="s2_path", help='gzipped or unzipped fastq file of the simulated '
                             's2 reads')
    parser.add_argument(
        "--mergecsv", action="store", type=str,  required=True, 
        dest="merge_info_path", help='csv containing merge information')
    parser.add_argument(
        "--fraglen", action="store", type=int, required=True,
        help="fraglen") 
    parser.add_argument(
        "--nfrags", action="store", type=int, required=True,
        help="nfrags")  
    parser.add_argument(
        "--tool", action="store", type=str, required=True,
        dest="tool_name", help="Name of the program used for trimming")
    parser.add_argument(
        "--qualityshift", action="store", type=int, required=True,
        help="the quality shift used when simulating the reads (when shifting "
             "scores by x, the error rate is 1/(10^(x/10)) of the default "
             "profile.") 

    # optional arguments, only needed if exporting the results
    parser.add_argument(
        "-o", "--out", action="store", type=str, required=False,
        dest="export_path", help="Path for the output csv file")

    args = parser.parse_args()
    required_arguments = [
        args.templates_path, 
        args.s1_path,
        args.s2_path,
        args.merge_info_path, 
        args.nfrags, 
        args.fraglen,
        args.qualityshift,
        args.tool_name
        ]
    optional_arguments = [args.export_path]
    
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


def _next_fastq_entry(fastq_file, fraglen):
    entry = dict()
    entry["name"] = fastq_file.readline().rstrip()[1:-4]
    entry["sequence"] = fastq_file.readline().rstrip()[:fraglen]
    entry["optional"] = fastq_file.readline().rstrip()
    entry["quality"] = fastq_file.readline().rstrip()[:fraglen]
    return entry


def _next_fasta_entry(fasta_file):
    entry = dict()
    entry["name"] = fasta_file.readline().rstrip()[1:]
    entry["sequence"] = fasta_file.readline().rstrip()
    return entry


def get_next_entries(fragments_file, s1_file,s2_file, fraglen, dna_complement):
    
    orig = _next_fasta_entry(fragments_file)
    read1 = _next_fastq_entry(s1_file, fraglen)
    read2 = _next_fastq_entry(s2_file, fraglen)
    read2["sequence"] = read2["sequence"].translate(dna_complement)[::-1]
    read2["quality"] = read2["quality"][::-1]
    assert orig["name"] == read1["name"]
    assert orig["name"] == read2["name"]
    return orig, read1, read2


def merged_base_and_qual(merge_df, base1, base2, qual1, qual2):
    # find out if the base matches in read1 and read2
    if base1 == base2:
        match_type = "match"
    else:
        match_type = "mismatch"
    # get the right row of the data frame
    mask = (merge_df["type"] == match_type) \
        & (merge_df["qs1"] == qual1) \
        & (merge_df["qs2"] == qual2)
    row = merge_df[mask]
    # get the merged base
    if row["new_nt"].iloc[0] == "s1":
        merged_base = base1
    else: 
        merged_base = base2
    # get merged quality score
    merged_qual = row["new_qs"].iloc[0]
    return merged_base, merged_qual


def add_to_phred_counter(phred_counter, merged_qual, merged_base, orig_base):
    # make sure the merged phred is in the dict
    if merged_qual not in phred_counter:
        phred_counter[merged_qual] = {
            "total_cnt": 0, 
            "error_cnt": 0
            }
    # add to counter
    phred_counter[merged_qual]["total_cnt"] += 1
    if merged_base != orig_base:
        phred_counter[merged_qual]["error_cnt"] += 1


def process_reads(orig, read1, read2, merge_df, phred_counter):
    for pos in range(len(orig["sequence"])):
        orig_base = orig["sequence"][pos]
        base1 = read1["sequence"][pos]
        base2 = read2["sequence"][pos]
        # get the quality score of read1 and read2
        qual1 = read1["quality"][pos] - 33
        qual2 = read2["quality"][pos] - 33
        # get merged base and quality
        merged_base, merged_qual = merged_base_and_qual(merge_df, 
                                                        base1, base2, 
                                                        qual1, qual2)
        # add count for the quality score
        add_to_phred_counter(phred_counter, merged_qual, merged_base, orig_base)


def phred_2_p_error(q_value):
    p_error = 10**(-int(q_value)/10)
    return p_error


def p_mismatch(n_total, n_mismatch):
    if n_total == 0:
        p = np.NAN
    else:
        p = n_mismatch/n_total
    return p


def binomial_ci(n, k, alpha):
    """ 
    Exact Confidence Interval
    https://sigmazone.com/binomial-confidence-intervals/ 
    """
    if k == 0:
        p_lower = 0
    else: 
        p_lower = 1 - st.beta.ppf(1-(alpha/2), n-k+1 , k)
    if k == n:
        p_upper = 1
    else:
        p_upper = 1 - st.beta.ppf(alpha/2, n-k , k+1)
    return p_lower, p_upper


def p_error_2_phred(p_err, max_phred=100):
    # probability of incorrect base call
    if p_err <= 0:
        phred_observed = max_phred
    elif p_err == 1:
        phred_observed = 0
    else: 
        phred_observed = -10 * math.log10(p_err)
    return phred_observed


def get_results(phred_counter, alpha):
    # set a maximum value for the observed phred, this is needed if the
    # observed error rate is 0
    max_phred = 100
    
    results = {'predicted_phred': list(),
        'predicted_error': list(),
        'n_matches': list(),
        'n_mismatches': list(),
        'n_total': list(),
        'p_mismatch': list(),
        'p_mismatch_lower': list(),
        'p_mismatch_upper': list(),
        'observed_phred': list(),
        'observed_phred_lower': list(),
        'observed_phred_upper': list(),
        }
    
    for q_score in sorted(phred_counter.keys()):
        
        predicted_error = phred_2_p_error(q_score)
        n_total = phred_counter[q_score]["total_cnt"]
        n_mismatches = phred_counter[q_score]["error_cnt"]
        p_mm = p_mismatch(n_total, n_mismatches)
        p_mm_lower, p_mm_upper = binomial_ci(n_total, n_mismatches, alpha)
        observed_phred = p_error_2_phred(p_mm, max_phred)
        observed_phred_lower = p_error_2_phred(p_mm_upper, max_phred)
        observed_phred_upper = p_error_2_phred(p_mm_lower, max_phred)
    
        results["predicted_phred"].append(q_score)
        results["predicted_error"].append(predicted_error)
        results["n_matches"].append(n_total - n_mismatches)
        results["n_mismatches"].append(n_mismatches)
        results["n_total"].append(n_total)
        results["p_mismatch"].append(p_mm)
        results["p_mismatch_lower"].append(p_mm_lower)
        results["p_mismatch_upper"].append(p_mm_upper)
        results["observed_phred"].append(observed_phred)
        results["observed_phred_lower"].append(observed_phred_lower)
        results["observed_phred_upper"].append(observed_phred_upper)
    
    return results


def main(template_path, s1_path, s2_path, merge_info_path, nfrags, fraglen, qualityshift, tool_name,
         export_path=None):

    alpha = 0.01

    # Load files
    templates_file = gzip.open(template_path, "rb") if common._is_gzipped(template_path) else open(template_path, "rb")
    s1_file = gzip.open(s1_path, "rb") if common._is_gzipped(s1_path) else open(s1_path, "rb")
    s2_file = gzip.open(s2_path, "rb") if common._is_gzipped(s2_path) else open(s2_path, "rb")
    merge_df = pd.read_csv(merge_info_path)
    merge_df = merge_df[merge_df["program"] == tool_name]

    # Analysis
    phred_counter = dict()
    dna_complement = bytes.maketrans(b"ACTG", b"TGAC")
    for _ in range(nfrags):
        orig, read1, read2 = get_next_entries(templates_file, s1_file, s2_file, fraglen, dna_complement)
        process_reads(orig, read1, read2, merge_df, phred_counter)
    
    results = get_results(phred_counter, alpha)


    # Export results 
    if export_path is not None:
        df = pd.DataFrame.from_dict(results)
        df.insert(0, 'program', tool_name)
        df.insert(1, 'nfrags', nfrags)
        df.insert(2, 'fraglen', fraglen)
        df.insert(3, 'qual_shift', qualityshift)
        df.insert(4, 'alpha', alpha)
        df.to_csv(export_path, na_rep="NA")


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)
    