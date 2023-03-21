from matplotlib import pyplot as plt
import pandas as pd
import os


def plot_simulated_phred_occurence(df, outdir):
    
    
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    
    reads = sorted(df["read"].unique())
    quality_shifts = sorted(df["quality_shift"].unique(), reverse=True)
    
    
    for i in range(len(reads)):
        read = reads[i]
        df_read = df[df["read"] == read]
        read_verbos = "forward" if read == "s1" else "reverse"
        for j in range(len(quality_shifts)):
            q_shift = quality_shifts[j]
            df_qs = df_read[df_read["quality_shift"] == q_shift]
            
            ax = axes[j][i]
            ax.bar(df_qs["quality_score"], df_qs["count"])
    
            ax.set_xlim(-1, 42)
            ax.set_xticks(range(0, 42, 5))
            ax.set_xticks(range(0, 42), minor=True)
            #ax.xaxis.set_tick_params(labelsize=5)
            ax.set_xlabel(f"Phred quality score")
            ax.set_ylabel(f"Total count")
            ax.set_title(f"{read_verbos} read, quality shift {q_shift},")
            ax.grid(axis="y", alpha=0.5)
            
    
    # Save plot
    fig.tight_layout()
    if outdir is not None:
        plt.savefig(f"{outdir}/phred_count_simulated_reads.png", 
                    dpi='figure', 
                    format="png")
    else:
        plt.show()
    plt.close(fig)
    
    
infile = os.path.join("output", "phred_count_simulated_reads.csv")
outdir = os.path.join("output", "plots", "final")
df = pd.read_csv(infile)
plot_simulated_phred_occurence(df, outdir)