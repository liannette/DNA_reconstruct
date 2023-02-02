import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots the pca")
    
    # required arguments
    parser.add_argument(
        "-evel", action="store", type=str, required=True, dest="evec_file", 
        help='path of the evec file')
    parser.add_argument(
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args
    
    
    
df = read_csv()


if __name__ == "__main__":
    
    