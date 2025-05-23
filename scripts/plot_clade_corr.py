import os
import argparse
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
    clades = sorted(clades_dico.keys(), key=lambda x: sum((("mammal" in x) * 0.5, "NoPrimates" in x, ("primate" in x) * 2)))
    list_pairs = list(combinations(clades, 2))
    datasets = sorted(datasets)
    fig, axs = plt.subplots(nrows=len(list_pairs), ncols=len(datasets), sharex='col',
                            figsize=(5 * len(datasets), 4 * len(list_pairs)))
    for x_1, (clade_1, clade_2) in enumerate(list_pairs):
        df_out = pd.merge(clades_dico[clade_1], clades_dico[clade_2], on=["trait", "dataset"],
                          suffixes=(f"_{clade_1}", f"_{clade_2}"), how="inner")
        gpby = {dataset: df_gr for dataset, df_gr in df_out.groupby("dataset")}

        for x_2, dataset in enumerate(datasets):
            if dataset not in gpby:
                continue
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs[x_1]
            ax.set_xlabel(f"Ratio {clade_1}")
            ax.set_ylabel(f"Ratio {clade_2}")
            ax.axvline(0, linestyle="--", color="black", alpha=0.5)
            ax.axhline(0, linestyle="--", color="black", alpha=0.5)
            #ax.set_xscale("log")
            #ax.set_yscale("log")
            df_gr = gpby[dataset]
            # Show as 2d histogram
            x = df_gr[f"ratio_{clade_1}"]
            y = df_gr[f"ratio_{clade_2}"]
            f = (~(np.isnan(x) | np.isnan(y)) & (x > 0) & (y > 0))
            # Clip to -5 <-> 2
            x = np.clip(np.log(x[f]), -5, 4)
            y = np.clip(np.log(y[f]), -5, 4)
            ax.hist2d(x, y, bins=30, cmap="Blues", cmin=1)
            # Bins
            # df = df_gr.groupby("ratioqcut", observed=False).agg({f"ratio_{clade_1}": "mean", f"ratio_{clade_2}": "mean"}).reset_index()
            # df[f"ln_ratio_{clade_1}"] = df[f"ratio_{clade_1}"].apply(lambda x: np.log(x) if x > 0 else np.nan)
            # df[f"ln_ratio_{clade_2}"] = df[f"ratio_{clade_2}"].apply(lambda x: np.log(x) if x > 0 else np.nan)
            # corr = df[f"ln_ratio_{clade_1}"].corr(df[f"ln_ratio_{clade_2}"])
            # ax.plot(df[f"ratio_{clade_1}"], df[f"ratio_{clade_2}"], "o", alpha=0.5, label=f"r={corr:.2g}")
            ax.set_title(dataset)
            # ax.legend()
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
