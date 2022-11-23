import argparse
import pandas as pd
from matplotlib import pyplot as plt
import copy
import math


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


def get_edit_distance_matrix(df_program):
    # convert the edit_distance strings to a matrix for bar plotting
    nfrags = list(df_program["nfrags"])[0]
    edit_distance_matrix = []
    for row in df_program["edit_distances"]:
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


def break_xaxis(ax, ax2, xlim_ax1, xlim_ax2):
    ax.set_xlim(*xlim_ax1)
    ax2.set_xlim(*xlim_ax2)

    # hide the spines between ax and ax2
    ax.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    ax.yaxis.tick_left()
    ax2.yaxis.tick_right()

    # Draw the diagonal lines
    d = 0.7  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax.plot([1, 1], [1, 0], transform=ax.transAxes, **kwargs)
    ax2.plot([0, 0], [1, 0], transform=ax2.transAxes, **kwargs)

    # we can vary the distance between
    # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
    # the diagonal lines will move accordingly, and stay right at the tips
    # of the spines they are 'breaking'
    

def plot_edit_distance_plot(program, frag_len, dropped_percent, edit_distances, 
                            outdir):

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', 
                                   figsize=[12, 7], 
                                   gridspec_kw={'width_ratios': [0.98, 0.02]})
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
    
    for ax in (ax1, ax2):
        bottom = len(dropped_percent) * [0]
        # Plot divergent and perfectly reconstructed reads
        for i, percent in list(reversed(edit_distances)):
            ax.bar(frag_len, percent, width, bottom=bottom, 
                   label=labels[i], color=colors[i])
            bottom = [sum(x) for x in zip(bottom, percent)]
        # Plot unmerge reads
        ax.bar(frag_len, dropped_percent, width, bottom=bottom,
               label=labels[-1], color=colors[-1])
        bottom += [sum(x) for x in zip(bottom, dropped_percent)]


    # break axis
    fig.subplots_adjust(wspace=0.03)
    break_xaxis(ax1, ax2, (0, 253), (997, 1001))
    # Set y axis limit
    ax1.set_ylim(0, max(bottom))
    # Add ticks
    ax2.set_xticks([1000])
    ax1.set_xticks(range(0,251,10))
    ax1.set_yticks(range(0,101,10))
    # Add a line at the read length
    ax1.plot([125.5, 125.5], [0, 100], color='green', linestyle='--', lw=1)
    ax1.text(124, 80, f"read length", color='green', fontsize=7,
             rotation=90, rotation_mode='anchor')
    # Add grid
    ax1.grid(alpha=0.5)
    ax2.grid(alpha=0.5)
    # Add labels, title
    ax1.set_xlabel('Fragment lengths')
    ax1.set_ylabel('Percentage of simulated fragments')
    ax1.set_title(f"{program}, merged reads")
    # Add legend and change order
    handles, labels = plt.gca().get_legend_handles_labels()
    order = list(reversed(range(len(edit_distances)+1)))
    ax2.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='center left', bbox_to_anchor=(1, 0.5))

    fig.tight_layout()
    plt.savefig(f"{outdir}/edit_distances_{program}.png", 
                dpi='figure', 
                format="png")


def main(infile, outdir):
    
    df = pd.read_csv(infile)
    
    for program in list(df["program"].unique()):
        df_program = df[df["program"] == program]
        
        frag_len = df_program["fraglen"]
        nfrags = df_program["nfrags"]
        dropped_percent = df_program["dropped_reads"] / nfrags * 100
        edit_distances = get_edit_distance_matrix(df_program)
        
        plot_edit_distance_plot(
            program, 
            frag_len, 
            dropped_percent, 
            edit_distances, 
            outdir
            )
    
            
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)