import argparse
import pandas as pd
from matplotlib import pyplot as plt
import copy


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


def get_edit_distance_matrix(df_program, num_categories):
    # convert the edit_distance strings to a matrix for bar plotting
    nfrags = list(df_program["nfrags"])[0]
    edit_distance_matrix = []
    for row in df_program["edit_distances"]:
        # create the row of the matrx
        occurences_per_edit_distance = num_categories * [0]
        for element in row.split():
            edit_dist, cnt = [int(x) for x in element.split(":")]
            percent = round(cnt/nfrags*100, 3)
            # edit distances bigger than the matrix will be put into the last row
            if edit_dist >= num_categories:
                occurences_per_edit_distance[-1] += percent
            else:
                occurences_per_edit_distance[int(edit_dist)] = percent
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
        

def plot_read_lengths(program, frag_len, dropped, incorrect_len, 
                      correct_len_incorrect_seq, perfect, outdir):
      
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w', 
                                   figsize=[12, 7], 
                                   gridspec_kw={'width_ratios': [0.98, 0.02]})
    
    for ax in (ax1, ax2):
        width = 1
        ax.bar(frag_len, dropped, width, label='not merged', color='#57106e')
        bottom = copy.deepcopy(dropped)
        ax.bar(frag_len, incorrect_len, width, bottom=bottom, 
               label='incorrect length', color='#bc3754')
        bottom += incorrect_len
        ax.bar(frag_len, correct_len_incorrect_seq, width, bottom=bottom, 
              label='correct len, incorrect seq', color='#f98e09')
        bottom += correct_len_incorrect_seq
        ax.bar(frag_len, perfect, width, bottom=bottom, 
              label='perfectly reconstructed', color='#fcffa4')
        bottom += perfect

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
    ax1.text(124, 80, f"raw read length", color='green', fontsize=7,
             rotation=90, rotation_mode='anchor')
    # Add grid
    ax1.grid(alpha=0.5)
    ax2.grid(alpha=0.5)
    # Add labels, title
    ax1.set_xlabel('Fragment lengths')
    ax1.set_ylabel('Percentage')
    ax1.set_title(f"{program}, merged reads")
    # Add legend and change order
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [3, 2, 1, 0]
    ax2.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='center left', bbox_to_anchor=(1, 0.5))

    fig.tight_layout()
    plt.savefig(f"{outdir}/fraglen_{program}.png", 
                dpi='figure', 
                format="png")
    #plt.show()
    

def main(infile, outdir):
    
    df = pd.read_csv(infile)
    df = df.sort_values(by=['fraglen'])
    for program in list(df["program"].unique()):
        
        # Only take data for one program 
        df_program = df[df["program"] == program]
        
        frag_len = df_program["fraglen"]
        dropped = df_program["dropped_reads_percentage"]
        incorrect_len = df_program["incorrect_length_percentage"]
        correct_len_incorrect_seq = df_program["correct_len_incorrect_seq_percentage"]
        perfect = df_program["perfectly_reconstructed_percentage"]
        
        plot_read_lengths(program, frag_len, dropped, incorrect_len,
                          correct_len_incorrect_seq, perfect, outdir)
    
            
if __name__ == "__main__":

    args = parse_arguments()
    main(*args)