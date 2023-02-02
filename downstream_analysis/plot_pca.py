import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plots the pca")
    
    # required arguments
    parser.add_argument(
        "-eval", action="store", type=str, required=True, dest="evec_file", 
        help='path to the evec file from the smartPCA output')
    parser.add_argument(
        "-p1", action="store", type=str, required=True, dest="pop_pca_file", 
        help='path to the population list file, that was used to generate the '
        'principal components')
    parser.add_argument(
        "-p2", action="store", type=str, required=True, 
        dest="pop_projected_file", help='path of the evec file')
    parser.add_argument(
        "-o", action="store", type=str, required=False, dest="outdir",
        help="path of the output directory")  

    args = parser.parse_args()
    return args
    
    
def load_pca_data(evec_file):
    col_names = ["Name", "PC1", "PC2", "PC3", "PC4", "Population"]
    df_pca = pd.read_csv(evec_file, sep="\t", header=0, names=col_names)
    return df_pca


def load_population_list(pop_pca_file):
    """ Reads the csv file of the populations that have been used to generate
    the principal components. Adds color and symbol indeces to each 
    population for plotting """
    df_population = pd.read_csv(
        pop_pca_file, names=["Population"]).sort_values(by="Population")
    n_populations= len(pop_pca_list)
    n_colors = 8
    n_symbols = int(n_populations / n_colors)
    color_indices = [int(i / n_symbols) for i in range(n_populations)]
    symbol_indices = [i % n_symbols for i in range(n_populations)]
    df_population = df_population.assign(color_index=color_indices, 
                                         symbol_index=symbol_indices)
    return df_population


def load_population_list(pop_pca_file):
    """ Reads the csv file of the populations that have been used to generate
    the principal components. Adds color and symbol indeces to each 
    population for plotting """
    df_population = pd.read_csv(
        pop_pca_file, names=["Population"]).sort_values(by="Population")
    n_populations= len(pop_pca_list)
    n_colors = 8
    n_symbols = int(n_populations / n_colors)
    color_indices = [int(i / n_symbols) for i in range(n_populations)]
    symbol_indices = [i % n_symbols for i in range(n_populations)]
    df_population = df_population.assign(color_index=color_indices, 
                                         symbol_index=symbol_indices)
    return df_population


def plot_pca(df_pca, df_population):
    """ flipping the x axis to make the correlation to Geography more obvious
    """
    
    symbols = ["8", "s", "p", "P", "*", "h", "H", "+", "x", "X", "D", "d"]
    colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd',
                u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    
    plt.figure(figsize=(10, 10))
    
    # The populations
    for i, row in df_population.iterrows():
        d = df_pca[df_pca.Population == row["Population"]]
        plt.scatter(x=-d["PC1"], y=d["PC2"], 
                    c=colors[row["color_index"]],
                    marker=symbols[row["symbol_index"]], 
                    label=row["Population"])
    # The samples HG002 and HG005
    for i, row in df_samples.iterrows():
        d = df_pca[df_pca.Samples == row["Population"]]




    plt.xlabel("PC1")
    plt.ylabel("PC2")
    
    plt.legend(loc=(1.1, 0), ncol=3)


if __name__ == "__main__":
    args = parse_arguments()
    df_pca = load_pca_data(args.evec_file)
    df_population = load_population_list(args.pop_pca_file)
    
    
    