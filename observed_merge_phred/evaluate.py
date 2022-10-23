#!/usr/bin/python3

import argparse
import sys
import math
import pandas as pd
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import common
import numpy as np
import scipy.stats as st


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
        "-qs", "--qualityshift", action="store", type=int, required=True,
        help="the quality shift used when simulating the reads (when shifting "
             "scores by x, the error rate is 1/(10^(x/10)) of the default "
             "profile.") 
    
    # optional arguments, only needed if exporting the results
    parser.add_argument(
        "-o", "--out", action="store", type=str, required=False,
        dest="export_path", help="Path for the output csv file")
    parser.add_argument(
        "-t", "--tool", action="store", type=str, required=False,
        dest="tool_name", help="Name of the program used for trimming")

    args = parser.parse_args()
    required_arguments = [
        args.templates_path, 
        args.readm_path,
        args.nfrags, 
        args.fraglen,
        args.qualityshift,
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
        

def process_merged_reads(merged_reads, templates):
    phred_counter = dict()

    for read in merged_reads:
        name = read['name']
        orig_seq = templates[name]['sequence']
        read_seq = read['sequence']
        read_quality = read['quality']
        
        # try:
        #     if len(orig_seq) != len(read_seq):
        #         raise Exception("Error: Read with different length")
        # except Exception as e:
        #     print(e)
        #     print(name)
        #     sys.exit(1)
        
        # Make sure that the merged read can be compared with the orig 
        if len(orig_seq) == len(read_seq):
            for i in range(len(read_seq)):
                nt = read_seq[i]
                qual = read_quality[i]-33
                
                if qual not in phred_counter:
                    phred_counter[qual] = {
                        "match_cnt": 0, 
                        "mismatch_cnt": 0
                        }

                if nt == orig_seq[i]:
                    phred_counter[qual]["match_cnt"] += 1
                else: 
                    phred_counter[qual]["mismatch_cnt"] += 1

    return phred_counter


def p_mismatch(n_total, n_mismatch):
    if n_total == 0:
        p = np.NAN
    else:
        p = n_mismatch/n_total
    return p


def margin_of_error(n_total, p, alpha=0.01):
    """ value of maximum error given alpha for normal distribution """
    # z score of normal distribution
    z_value = st.norm.ppf(1 - alpha/2)
    # standard error of an estimate of a sample proportion
    stderr = math.sqrt( p*(1-p) / n_total )
    return z_value * stderr


def p_error_2_phred(p_err, max_observed_phred=100):
    # probability of incorrect base call
    if p_err <= 0:
        phred_observed = max_observed_phred
    elif p_err == 1:
        phred_observed = 0
    else: 
        phred_observed = -10 * math.log10(p_err)
    return phred_observed


def get_results(phred_counter, alpha):
    # set a maximum value for the observed phred, this is needed if the
    # observed error rate is 0
    max_phred = 100
    
    results = []
    for phred_value in sorted(phred_counter.keys()):
        n_matches = phred_counter[phred_value]["match_cnt"]
        n_mismatches = phred_counter[phred_value]["mismatch_cnt"]
        n_total = n_mismatches + n_matches
        p_mm = p_mismatch(n_total, n_mismatches)
        err_margin = margin_of_error(n_total, p_mm)
        observed_phred_mean = p_error_2_phred(p_mm, max_phred)
        observed_phred_lower = p_error_2_phred(p_mm+err_margin, max_phred)
        observed_phred_higher = p_error_2_phred(p_mm-err_margin, max_phred)
            
        results.append(
            [
                phred_value, 
                n_matches, 
                n_mismatches, 
                n_total, 
                observed_phred_mean, 
                observed_phred_lower,
                observed_phred_higher,
                alpha
                ]
            )
    return results


def main(template_path, readm_path, nfrags, fraglen, qualityshift, 
         export_path=None, tool_name=None):

    # Load files --------------------------------------------------------------

    templates = common.load_fasta(template_path)

    # Check for duplicate fragments
    if len(templates) != nfrags:
        print(f"ATTENTION: number of total_sequences is {len(templates)}," 
              f"but the nfrags is {nfrags}. Possible reason: duplicate "
              "fragments")

    reads = common.load_fastq(readm_path)
    # seperator: this character and all charaters to the right of it
    # will be removed from the fastq header
    seperator = b'-'
    merged_reads = common.clean_merged_reads(reads, templates, seperator)


    # Analysis ----------------------------------------------------------------

    phred_counter = process_merged_reads(merged_reads, templates)
    results = get_results(phred_counter, alpha=0.01)


    # Export results ----------------------------------------------------------
    
    if export_path is not None:

        df = pd.DataFrame(
            results, 
            columns=[
                'predicted_phred', 
                'n_matches', 
                'n_mismatches', 
                'n_total', 
                'observed_phred_mean'
                'observed_phred_lower'
                'observed_phred_higher'
                'alpha'
                ]
            )
        df.set_index('predicted_phred', inplace=True)
        df.insert(0, 'program', tool_name)
        df.insert(1, 'nfrags', nfrags)
        df.insert(2, 'fraglen', fraglen)
        df.insert(3, 'qual_shift', qualityshift)
        df.to_csv(export_path, na_rep="NA")


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)