from matplotlib import pyplot as plt
import pandas as pd
import os


def plot_simulated_phred_occurence(df, outdir):
    
    
    fig, axes = plt.subplots(4, 2, figsize=(20, 12))
    
    reads = sorted(df["read"].unique())
    quality_shifts = sorted(df["quality_shift"].unique(), reverse=True)
    
    
    for i in range(len(reads)):
        read = reads[i]
        df_read = df[df["read"] == read]
        for j in range(len(quality_shifts)):
            q_shift = quality_shifts[j]
            df_qs = df_read[df_read["quality_shift"] == q_shift]
            
            ax = axes[j][i]
            ax.bar(df_qs["quality_score"], df_qs["count"])
    
            ax.set_xlim(-1, 42)
            ax.set_xticks(range(0, 42))
            #ax.xaxis.set_tick_params(labelsize=5)
            ax.set_xlabel(f"Simulated phred quality score")
            ax.set_ylabel(f"Total count")
            ax.set_title(f"Quality shift {q_shift}, {read}")
            ax.grid(axis="y", alpha=0.5)
            
    fig.suptitle("Occurence of quality scores in simulated data", 
                 fontsize=16)
    
    # Save plot
    fig.tight_layout()
    if outdir is not None:
        plt.savefig(f"{outdir}/simulated_phred_occurence.png", 
                    dpi='figure', 
                    format="png")
    else:
        plt.show()
    plt.close(fig)
    
    
infile = os.path.join("output", "evaluation", "simulated_reads.csv")
outdir = os.path.join("output", "plots")
df = pd.read_csv(infile)
plot_simulated_phred_occurence(df, outdir)