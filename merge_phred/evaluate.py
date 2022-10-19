#!/usr/bin/python3

import sys, os
import gzip
import pandas as pd
import argparse
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import common 


def parse_arguments():
    """
    """
    parser = argparse.ArgumentParser(
        description="Performs analysis of quality score merging.")
    
    # required arguments
    parser.add_argument(
        "-s1", action="store", type=str, required=True, dest="s1_path", 
        help='gzipped or unzipped fastq file of the initial forward reads')
    parser.add_argument(
        "-s2", action="store", type=str,  required=True, dest="s2_path", 
        help='gzipped or unzipped fastq file of the initial reverse reads')
    parser.add_argument(
        "-m", action="store", type=str, required=True, dest="merged_path", 
        help="gzipped or unzipped fastq file of the merged reads")
    parser.add_argument(
        "-o", action="store", type=str, required=True, dest="out_path",
        help="path for the csv result file")  
    parser.add_argument(
        "-t", "--tool", action="store", type=str, required=True,
        dest="program_name", help="Name of the program used for merging")

    args = parser.parse_args()
    arguments = [
        args.s1_path, 
        args.s2_path, 
        args.merged_path,
        args.out_path, 
        args.program_name,
        ]
    return arguments


def reverse_complement_bytes(seq):
    tab_b = bytes.maketrans(b"ACTG", b"TGAC")
    return seq.translate(tab_b)[::-1]


def load_initial_fastq(path, rev_complement = False):
    """
    Loads a zipped or unzipped fastq file.
    Returns a dict file, with the cleaned up headers as keys and another
    dict as value. The nested dict has the keys "sequence" and 
    "quality", the real sequence and quality scores (values) are limited
    to 31 characters. This is the size of the fragment.
    """
    seqs = {}
    f = gzip.open(path, 'rb') if common._is_gzipped(path) else open(path, 'rb')
    lines = []
    for line in f:
        lines.append(line.rstrip())
        if len(lines) == 4:
            header = common._clean_up_fastq_header(lines[0], b'/')
            if rev_complement:
                sequence =  reverse_complement_bytes(lines[1][0:31])
                quality = lines[3][0:31][::-1]
            else:
                sequence =  lines[1][0:31]
                quality = lines[3][0:31]
            seqs[header] = {
                'sequence': sequence, 
                'quality': quality
            }
            lines = []
    f.close()
    return seqs


def analyze_merged_reads(merged_reads, s1_seqs, s2_seqs):

    incorrect_length_cnt = 0
    matching_nt = []
    mismatching_nt = []

    for read in merged_reads:
        name = read['name']
        
        # indexing a byte string retrieves the integer form of the byte      
        s1_nt = chr(s1_seqs[name]['sequence'][15])
        s2_nt = chr(s2_seqs[name]['sequence'][15])
        merged_nt = chr(read['sequence'][15]) 
        s1_qual = s1_seqs[name]['quality'][15] - 33
        s2_qual = s2_seqs[name]['quality'][15] - 33
        merged_qual = read['quality'][15]-33 
            
        if len(read['sequence']) == 31:
            nt_info = [
                name.decode("utf-8"), 
                s1_nt, s2_nt, merged_nt, 
                s1_qual, s2_qual, merged_qual
            ]
            if s1_nt == s2_nt:
                matching_nt.append(nt_info)
            else:
                mismatching_nt.append(nt_info)
        else:
            incorrect_length_cnt += 1
    return matching_nt, mismatching_nt, incorrect_length_cnt



def main(s1_path, s2_path, merged_path, out_path, program_name):
    
    # Load files --------------------------------------------------------------

    # initial fastq files
    s1_seqs = load_initial_fastq(s1_path)
    s2_seqs = load_initial_fastq(s2_path, rev_complement = True)

    # Merged reads
    merged_reads = common.load_fastq(merged_path)
    # seperator: this character and all charaters to the right of it
    # will be removed from the fastq header
    seperator = b'/'
    merged_reads = common.clean_merged_reads(merged_reads, s1_seqs, seperator)
    
    # Analysis and Results ----------------------------------------------------

    result = analyze_merged_reads(merged_reads, s1_seqs, s2_seqs)
    matching_nt = result[0]
    mismatching_nt = result[1]
    incorrect_length_cnt = result[2]
    
    if incorrect_length_cnt > 0:
        # This should not happen
        try:
            raise Exception("Merged read of incorrect length")
        except Exception as e:
            print(e)
            

    # Export results -------------------------------------------------
        
    print(f"{os.path.basename(merged_path)}:")
    print(f"total seqs: {len(s1_seqs)}")
    print(f"total merged: {len(merged_reads)}")
    print(f"matching count: {len(matching_nt)}")
    print(f"mismatching count: {len(mismatching_nt)}")
    print(f"incorrect length count: {incorrect_length_cnt}")

    df1 = pd.DataFrame(
        matching_nt, 
        columns=['name', 'nt1', 'nt2', 'new_nt', 'qs1', 'qs2', 'new_qs']
        )
    df1.set_index('name', inplace=True)
    df1.insert(0, 'program', program_name)
    df1.insert(1, 'type', 'match')

    df2 = pd.DataFrame(
        mismatching_nt, 
        columns=['name', 'nt1', 'nt2', 'new_nt', 'qs1', 'qs2', 'new_qs']
        )
    df2.set_index('name', inplace=True)
    df2.insert(0, 'program', program_name)
    df2.insert(1, 'type', 'mismatch')

    df3 = pd.concat([df1, df2], axis=0)
    if df3.duplicated().any():
        try:
            raise Exception("Duplicates in the results")
        except Exception as e:
            print(e)
    df3.to_csv(out_path, na_rep="NA")

    print("Data exported sucessfully\n")


if __name__ == "__main__":

    args = parse_arguments()
    main(*args)
