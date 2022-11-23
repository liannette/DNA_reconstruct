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
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args.infile, args.outdir


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


def plot_edit_distance(df, outdir):
    fig, axes = plt.subplots(1, 7, figsize=[20, 7]) #, gridspec_kw={'width_ratios': [0.98, 0.02]})
    
    programs = list(df["program"].unique())
    for ax, program in zip(axes, programs):
        
        # get data
        df_program = df[df["program"] == program]
        quality_shift = df_program["quality_shift"]
        edit_distances = get_edit_distance_matrix(df_program)
        dropped_percent = df_program["dropped_reads"] / df_program["nfrags"] * 100
        
        labels = [
            "prefectly reconstructed", 
            "edit distance: 1",
            "edit distance: 2",
            "edit distance: 3",
            "edit distance: 4",
            "edit distance: 5-10",
            "edit distance: 11-20",
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
        width = 10
        
        bottom = len(df_program) * [0]
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
        ax.set_xlim(-33, 5)
        ax.set_ylim(0, max(bottom))
        ax.set_yticks(range(0,101,10))
        # Add grid
        ax.grid(axis="y", alpha=0.3)
        # Add title and xlabel
        ax.set_title(f"{program}")
        ax.set_xlabel('quality shift')
        ax.invert_xaxis()

    # Add suptitle and ylabel
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
    plt.savefig(f"{outdir}/droprate_qualityscores.png", 
                dpi='figure', 
                format="png")


def main(infile, outdir):
    
    df = pd.read_csv(infile)
    plot_edit_distance(df, outdir)
    
            
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)