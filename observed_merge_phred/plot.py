import argparse
import pandas as pd
from evaluate import p_mismatch, margin_of_error, p_error_2_phred
from matplotlib import pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots the heatmaps")
    
    # required arguments
    parser.add_argument(
        "-i", action="store", type=str, required=True, dest="infile", 
        help='path for the csv result file')
    parser.add_argument(
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args.infile, args.outdir


def plot_confidence_interval(x, mean, lower, upper, color='#2187bb', 
                             horizontal_line_width=0.6):
    left = x - horizontal_line_width / 2
    right = x + horizontal_line_width / 2
    plt.plot([x, x], [upper, lower], color=color, lw=0.6)
    plt.plot([left, right], [upper, upper], color=color, lw=0.6)
    plt.plot([left, right], [lower, lower], color=color, lw=0.6)
    plt.plot(x, mean, '.', color=color)


def plot_observed_vs_predicted_phred(df_program, program, alpha, outdir):

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot the expected line
    ax.plot(df_program.index, df_program.index, "--r", lw=0.8)
    # Plot each conf interval
    for x in df_program.index:
        mean = df_program["observed_phred_mean"][x]
        lower = df_program["observed_phred_lower"][x]
        upper = df_program["observed_phred_upper"][x]
        plot_confidence_interval(x, mean, lower, upper, color="k")

    # Add alpha information (used to calculate conf interval)
    ax.text(df_program.index[-1]-5, 2, f"alpha = {alpha}",)
    # Show all values on x axes
    plt.xticks(df_program.index, size="x-small")
    # Set x axes limits
    ax.set_xlim(left=df_program.index[0]-1, right=df_program.index[-1]+1) 
    # Add grid
    ax.set_axisbelow(True)
    plt.grid(alpha=0.5)
    # Add labels
    ax.set_title(f"{program}")
    ax.set_xlabel(f"Merged Phred")
    ax.set_ylabel(f"Observed Phred")

    # Save plot
    fig.tight_layout()
    plt.savefig(f"{outdir}/{program}.png", 
                dpi='figure', 
                format="png")
    

def main(infile, outdir):
    
    alpha = 0.01
    
    df = pd.read_csv(infile)
    for program in list(df["program"].unique()):
        
        # Only take data for one program 
        df_program = df[df["program"] == program]
        # add the counts of the datasets with the different qual_shift
        df_program = df_program[["predicted_phred", 
                                "n_matches", 
                                "n_mismatches", 
                                "n_total"]
                                ].groupby("predicted_phred").sum()
        
        # Calculate the new observed phred with confidence intervall
        df_program["p_mismatch"] = df_program.apply(
            lambda x: p_mismatch(x["n_total"], x["n_mismatches"]), 
            axis=1
            )
        df_program["error_margin"] = df_program.apply(
            lambda x: margin_of_error(x["n_total"], x["p_mismatch"], 
                                      alpha=alpha), 
            axis=1
            )
        df_program["observed_phred_mean"] = df_program.apply(
            lambda x: p_error_2_phred(x["p_mismatch"]), 
            axis=1
            )
        df_program["observed_phred_lower"] = df_program.apply(
            lambda x: p_error_2_phred(x["p_mismatch"]+x["error_margin"]),
            axis=1
            )
        df_program["observed_phred_upper"] = df_program.apply(
            lambda x: p_error_2_phred(x["p_mismatch"]-x["error_margin"]),
            axis=1
            )   
        plot_observed_vs_predicted_phred(df_program, program, alpha, outdir)
        
        
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)