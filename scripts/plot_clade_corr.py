import os
import argparse
from itertools import combinations
from libraries_plot import *


def main(tsv_input_list: list, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    clades_dico = {}
    datasets = set()
    for tsv_input in tsv_input_list:
        clade = os.path.basename(os.path.dirname(tsv_input)).split("_")[1]
        clades_dico[clade] = pd.read_csv(tsv_input, sep='\t')
        datasets.update(clades_dico[clade]["dataset"].unique())

    # First giving priority to mammals, then to noprimates, then to primates
    clades = sorted(clades_dico.keys(),
                    key=lambda x: sum((("mammal" in x) * 0.5, "NoPrimates" in x, ("primate" in x) * 2)))
    list_pairs = list(combinations(clades, 2))
    datasets = sorted(datasets)
    fig, axs = plt.subplots(nrows=len(list_pairs), ncols=len(datasets),
                            figsize=(5 * len(datasets), 4 * len(list_pairs)))
    for x_1, (clade_1, clade_2) in enumerate(list_pairs):
        df_out = pd.merge(clades_dico[clade_1], clades_dico[clade_2], on=["trait", "dataset"],
                          suffixes=(f"_{clade_1}", f"_{clade_2}"), how="inner")
        gpby = {dataset: df_gr for dataset, df_gr in df_out.groupby("dataset")}

        for x_2, dataset in enumerate(datasets):
            if dataset not in gpby:
                continue
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs[x_1]
            df_gr = gpby[dataset]
            plot_bins_box(ax, x_label=f"ratio_{clade_1}", y_label=f"ratio_{clade_2}", df=df_gr, q=4, scale="log")
            ax.set_xlabel(f"Ratio {clade_1}")
            ax.set_ylabel(f"Ratio {clade_2}")
            ax.set_title(dataset)
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv files", nargs="+", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
