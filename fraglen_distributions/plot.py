import argparse
import pandas as pd
from matplotlib import pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots the heatmaps")
    
    # required arguments
    parser.add_argument(
        "-i", action="store", type=str, required=True, dest="infile", 
        help='path for the csv result file')
    parser.add_argument(
        "-l", action="store", type=str, required=True, dest="dist_dir", 
        help='directory with the fragment length files')
    parser.add_argument(
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args.infile, args.dist_dir, args.outdir


def plot_histogram(dist_dir, distnames, outdir):

    fig, axes = plt.subplots(1, 4, figsize=[15, 5])

    for i in range(4):
        ax = axes[i]
        distname = distnames[i]
        infile = dist_dir + f"/{distname}.gz"
        df = pd.read_csv(infile, header = None) 
        
        x_min = 0
        if distname == "A9180":
            x_max = 100
        elif distname == "cfDNA":
            x_max = 600
        elif distname == "chagyrskaya8":
            x_max = 300
        elif distname == "Vi33.19":
            x_max = 175

        binwidth = int(x_max/60)

        ax.hist(df, density=1, bins=range(x_min, x_max + 10, binwidth))
        ax.set_title(f"{distname}")
        ax.set_xlabel('Fragment length')
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel('Probability')
        
    fig.suptitle(f"Fragment length distributions")
    fig.tight_layout()
    plt.savefig(f"{outdir}/fraglen_histograms.png", 
                dpi='figure', 
                format="png")


def get_edit_distance_matrix(df):
    # convert the edit_distance strings to a matrix for bar plotting
    nfrags = list(df["nfrags"])[0]
    edit_distance_matrix = []
    for row in df["edit_distances"]:
        # create the row of the matrx
        occurences_per_edit_distance = 8 * [0]
        
        # Skip rows that dont contain a string -> no merged reads
        if isinstance(row, str):
            for element in row.split():
                edit_dist, cnt = [int(x) for x in element.split(":")]
                percent = round(cnt/nfrags*100, 3)
                # edit distances bigger than the matrix will be put into the last row
                if edit_dist < 5:
                    occurences_per_edit_distance[int(edit_dist)] = percent
                elif edit_dist <= 10:
                    occurences_per_edit_distance[-3] += percent
                elif edit_dist <= 25:
                    occurences_per_edit_distance[-2] += percent
                else:
                    occurences_per_edit_distance[-1] += percent
                    
        edit_distance_matrix.append(occurences_per_edit_distance)
    return list(enumerate(zip(*edit_distance_matrix)))


def plot_edit_distance(program, df_program, distnames, outdir):
    
    fig, axes = plt.subplots(1, 4, figsize=[15, 7]) #, gridspec_kw={'width_ratios': [0.98, 0.02]})
    for ax, distname in zip(axes, distnames):
        
        # get data
        df_dist = df_program[df_program["fraglen_distribution"] == distname]
        quality_shift = df_dist["quality_shift"]
        edit_distances = get_edit_distance_matrix(df_dist)
        dropped_percent = df_dist["dropped_reads"] / df_dist["nfrags"] * 100
        
        labels = [
            "prefectly reconstructed", 
            "edit distance: 1",
            "edit distance: 2",
            "edit distance: 3",
            "edit distance: 4",
            "edit distance: 5-10",
            "edit distance: 11-25",
            "edit distance: >25",
            "not merged",
            ]
        colors = [
            "#fcffa4",
            "#fac228",
            "#f57d15",
            "#d44842",
            "#9f2a63",
            "#65156e",
            "#280b53",
            "#000004",
            "tab:grey",
            ]
        width = 1
        
        bottom = len(df_dist) * [0]
        # Plot divergent and perfectly reconstructed reads
        for i, percent in list(reversed(edit_distances)):
            ax.bar(quality_shift, percent, width, bottom=bottom, 
                   label=labels[i], color=colors[i])
            bottom = [sum(x) for x in zip(bottom, percent)]
        # Plot unmerge reads
        ax.bar(quality_shift, dropped_percent, width, bottom=bottom,
               label=labels[-1], color=colors[-1])
        bottom += [sum(x) for x in zip(bottom, dropped_percent)]
        

        # Set y axis limit
        ax.set_xlim(-30, 0)
        ax.set_ylim(0, max(bottom))
        ax.set_yticks(range(0,101,10))
        # Add grid
        ax.grid(axis="y", alpha=0.3)
        # Add title and xlabel
        ax.set_title(f"{distname}")
        ax.set_xlabel('quality shift')
        ax.invert_xaxis()

    # Add suptitle and ylabel
    plt.suptitle(f"{program}, merged reads")
    axes[0].set_ylabel('Percentage of simulated fragments')
    # Add legend and change order
    handles, labels = plt.gca().get_legend_handles_labels()
    order = list(reversed(range(len(edit_distances)+1)))
    axes[-1].legend([handles[idx] for idx in order], 
                    [labels[idx] for idx in order],
                    loc='center left', 
                    bbox_to_anchor=(1, 0.5))

    # And save it
    fig.tight_layout() 
    plt.savefig(f"{outdir}/fraglen_distributions_{program}.png", 
                dpi='figure', 
                format="png")


def main(infile, dist_dir, outdir):
    
    distnames = ["A9180", "Vi33.19", "chagyrskaya8", "cfDNA"]
    
    plot_histogram(dist_dir, distnames, outdir)
    
    df = pd.read_csv(infile)
    for program in list(df["program"].unique()):
        df_program = df[df["program"] == program]
        plot_edit_distance(program, df_program, distnames, outdir)
    
            
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)