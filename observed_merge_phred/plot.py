import argparse
import pandas as pd
from evaluate import p_mismatch, p_error_2_phred, phred_2_p_error
from matplotlib import pyplot as plt
import scipy.stats as st


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


def combine_quality_scores(df, alpha):
    df2 = df[["predicted_phred", "n_matches", "n_mismatches", "n_total"]] \
           .groupby("predicted_phred") \
           .sum()
    df2["predicted_phred"] = df2.index
    df2["predicted_p_mismatch"] = df2.apply(
        lambda x: phred_2_p_error(x["predicted_phred"]), axis=1)

    df2["p_mismatch"] = df2.apply(
        lambda x: p_mismatch(x["n_total"], x["n_mismatches"]), 
        axis=1
        )
    df2["p_mismatch_lower"] = df2.apply(
        lambda x: binomial_ci(x["n_total"], x["n_mismatches"], alpha)[0], 
        axis=1)
    df2["p_mismatch_upper"] = df2.apply(
        lambda x: binomial_ci(x["n_total"], x["n_mismatches"], alpha)[1], 
        axis=1
        )
    df2["observed_phred"] = df2.apply(
        lambda x: p_error_2_phred(x["p_mismatch"]), 
        axis=1
        )
    df2["observed_phred_lower"] = df2.apply(
        lambda x: p_error_2_phred(x["p_mismatch_upper"]),
        axis=1
        )
    df2["observed_phred_upper"] = df2.apply(
        lambda x: p_error_2_phred(x["p_mismatch_lower"]),
        axis=1
        )   
    return df2


def weighted_r_squared(df):
    """
    weighted mean: https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/weighted-mean/
    weighted r squared: https://stats.stackexchange.com/questions/83826/is-a-weighted-r2-in-robust-linear-model-meaningful-for-goodness-of-fit-analys/375752#375752 
    """
    y = df["observed_phred"]
    # Weights are the confidence interval
    conf_intervals = df["observed_phred_upper"] - df["observed_phred_lower"]
    weights = [1/w for w in conf_intervals]
    # difference between merged and observed phred
    residuals = y - df["predicted_phred"]
    
    # weighted residual sum of squared (SSe)
    sse = sum([w * e**2 for w, e in zip(weights, residuals)])
    # weighted total sum of squared (SSt)
    weighted_mean = sum([w * f for w, f in zip(weights, y)]) / sum(weights)
    sst = sum([w * (f-weighted_mean)**2 for w, f in zip(weights, y)])
    # weighted r squared
    r2 = 1 - (sse/sst)

    return r2


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


def plot_confidence_interval(x, y, lower, upper, color='#2187bb', 
                             horizontal_line_width=0.6):
    left = x - horizontal_line_width / 2
    right = x + horizontal_line_width / 2
    plt.plot([x, x], [upper, lower], color=color, lw=0.6)
    plt.plot([left, right], [upper, upper], color=color, lw=0.6)
    plt.plot([left, right], [lower, lower], color=color, lw=0.6)
    plt.plot(x, y, '.', color=color)
    
    
def plot_phred_calibration(df_qs, program, r2, alpha, outdir, qs=None):

    fig, ax = plt.subplots(figsize=(15, 8))

    # Plot the expected line
    ax.plot(df_qs["predicted_phred"], df_qs["predicted_phred"], "--r", lw=0.8)
    # Plot each conf interval
    idx = df_qs.index
    for i in idx:
        x = df_qs["predicted_phred"][i]
        y = df_qs["observed_phred"][i]
        lower = df_qs["observed_phred_lower"][i]
        upper = df_qs["observed_phred_upper"][i]
        plot_confidence_interval(x, y, lower, upper, color="k")
    # Add grid
    ax.set_axisbelow(True)
    plt.grid(alpha=0.5)
    # twin object for two different y-axis on the sample plot
    ax2=ax.twinx()
    # Plot the phred distribution
    ax2.bar(df_qs["predicted_phred"], df_qs["n_total"], alpha=0.2)

    # Add alpha information (used to calculate conf interval)
    ax.text(.01, .95, f"weighted R2 = {round(r2, 5)}", fontsize=14, 
            style='normal', ha='left', va='top', transform=ax.transAxes)
    ax.text(.01, .9, f"α={alpha}", fontsize=14, color='grey', 
            style='italic', ha='left', va='top', transform=ax.transAxes)
    
    # Show all values on x axes
    plt.xticks(range(df_qs["predicted_phred"][idx[0]], 
                     df_qs["predicted_phred"][idx[-1]] + 1), 
               size="xx-small")
    # Set x axes limits
    ax.set_xlim(left=df_qs["predicted_phred"][idx[0]]-0.5, 
                right=df_qs["predicted_phred"][idx[-1]]+0.5) 
    # Add labels
    if qs is None:
        plt.title(f"{program}", loc='left', fontsize=16)
    else:
        plt.title(f"{program}, quality shift: {qs}", loc='left', fontsize=16)
    ax.set_xlabel(f"Merged Phred")
    ax.set_ylabel(f"Observed error rate (phred scale)")
    ax2.set_ylabel("Total count")

    # Save plot
    fig.tight_layout()
    if outdir is not None:
        if qs is None:
            plt.savefig(f"{outdir}/{program}_combined.png", 
                dpi='figure', 
                format="png")
        else:
            plt.savefig(f"{outdir}/{program}_{qs}.png", 
                        dpi='figure', 
                        format="png")
    else:
        plt.show()
    plt.close(fig)
    
    
def main(infile, outdir):
    
    alpha = 0.01
    
    df = pd.read_csv(infile)
    
    # Only take data for one program 
    for program in list(df["program"].unique()):
        df_program = df[df["program"] == program]
        
        # For each quality shift
        for qs in list(df_program["qual_shift"].unique()):
            df_qs = df_program[df_program["qual_shift"] == qs]
            r2 = weighted_r_squared(df_qs)
            plot_phred_calibration(df_qs, program, r2, alpha, outdir, qs=qs)
        
        # combine the counts from the different quality shifts
        df_combined = combine_quality_scores(df_program, alpha)
        r2 = weighted_r_squared(df_combined)
        plot_phred_calibration(df_combined, program, r2, alpha, outdir, 
                               qs=None)



# def plot_observed_vs_predicted_phred(df_program, program, alpha, outdir):

#     fig, ax = plt.subplots(figsize=(12, 8))

#     # Plot the expected line
#     ax.plot(df_program.index, df_program.index, "--r", lw=0.8)
#     # Plot each conf interval
#     for x in df_program.index:
#         mean = df_program["observed_phred_mean"][x]
#         lower = df_program["observed_phred_lower"][x]
#         upper = df_program["observed_phred_upper"][x]
#         plot_confidence_interval(x, mean, lower, upper, color="k")

#     # Add alpha information (used to calculate conf interval)
#     ax.text(df_program.index[-1]-5, 2, f"alpha = {alpha}",)
#     # Show all values on x axes
#     plt.xticks(df_program.index, size="x-small")
#     # Set x axes limits
#     ax.set_xlim(left=df_program.index[0]-1, right=df_program.index[-1]+1) 
#     # Add grid
#     ax.set_axisbelow(True)
#     plt.grid(alpha=0.5)
#     # Add labels
#     ax.set_title(f"{program}")
#     ax.set_xlabel(f"Merged Phred")
#     ax.set_ylabel(f"Observed Phred")

#     # Save plot
#     fig.tight_layout()
#     plt.savefig(f"{outdir}/{program}.png", 
#                 dpi='figure', 
#                 format="png")
    

# def main(infile, outdir):
    
#     alpha = 0.01
    
#     df = pd.read_csv(infile)
#     for program in list(df["program"].unique()):
        
#         # Only take data for one program 
#         df_program = df[df["program"] == program]
#         # add the counts of the datasets with the different qual_shift
#         df_program = df_program[["predicted_phred", 
#                                 "n_matches", 
#                                 "n_mismatches", 
#                                 "n_total"]
#                                 ].groupby("predicted_phred").sum()
        
#         # Calculate the new observed phred with confidence intervall
#         df_program["p_mismatch"] = df_program.apply(
#             lambda x: p_mismatch(x["n_total"], x["n_mismatches"]), 
#             axis=1
#             )
#         df_program["error_margin"] = df_program.apply(
#             lambda x: margin_of_error(x["n_total"], x["p_mismatch"], 
#                                       alpha=alpha), 
#             axis=1
#             )
#         df_program["observed_phred_mean"] = df_program.apply(
#             lambda x: p_error_2_phred(x["p_mismatch"]), 
#             axis=1
#             )
#         df_program["observed_phred_lower"] = df_program.apply(
#             lambda x: p_error_2_phred(x["p_mismatch"]+x["error_margin"]),
#             axis=1
#             )
#         df_program["observed_phred_upper"] = df_program.apply(
#             lambda x: p_error_2_phred(x["p_mismatch"]-x["error_margin"]),
#             axis=1
#             )   
#         plot_observed_vs_predicted_phred(df_program, program, alpha, outdir)
        
        
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)