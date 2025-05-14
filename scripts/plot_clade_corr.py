import os
import argparse
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input_list: list, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    clades = {}
    datasets = set()
    for tsv_input in tsv_input_list:
        clade = os.path.basename(os.path.dirname(tsv_input)).split("_")[1]
        clades[clade] = pd.read_csv(tsv_input, sep='\t')
        datasets.update(clades[clade]["dataset"].unique())

    list_pairs = list(combinations(clades.keys(), 2))
    datasets = sorted(datasets)
    fig, axs = plt.subplots(nrows=len(list_pairs), ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4 * len(list_pairs)))
    for x_1, (clade_1, clade_2) in enumerate(list_pairs):
        df_out = pd.merge(clades[clade_1], clades[clade_2], on=["trait", "dataset"],
                          suffixes=(f"_{clade_1}", f"_{clade_2}"), how="inner")
        gpby = {dataset: df_gr for dataset, df_gr in df_out.groupby("dataset")}

        for x_2, dataset in enumerate(datasets):
            if dataset not in gpby:
                continue
            df_gr = gpby[dataset]
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs[x_1]
            ax.set_xlabel(f"Ratio {clade_1}")
            ax.set_ylabel(f"Ratio {clade_2}")
            ax.axvline(1, linestyle="--", color="black", alpha=0.5)
            ax.axhline(1, linestyle="--", color="black", alpha=0.5)
            ax.set_xscale("log")
            ax.set_yscale("log")
            if len(df_gr) < 2:
                continue
            corr = df_gr[f"ratio_{clade_1}"].corr(df_gr[f"ratio_{clade_2}"])
            ax.plot(df_gr[f"ratio_{clade_1}"], df_gr[f"ratio_{clade_2}"], "o", alpha=0.5, label=f"r={corr:.2g}")
            ax.set_title(dataset)
            ax.legend()
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
