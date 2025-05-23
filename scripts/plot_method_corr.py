import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input: list, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    assert len(tsv_input) == 2, f"Two input files are required, got {len(tsv_input)}"
    method_1 = os.path.basename(tsv_input[0]).split("_")[1].replace(".tsv", "").replace(".gz", "")
    method_2 = os.path.basename(tsv_input[1]).split("_")[1].replace(".tsv", "").replace(".gz", "")
    assert method_1 != method_2, f"Methods should be different, got {method_1} and {method_2}"
    df_1 = pd.read_csv(tsv_input[0], sep='\t')
    df_2 = pd.read_csv(tsv_input[1], sep='\t')
    df_out = pd.merge(df_1, df_2, on=["trait", "dataset"], suffixes=(f"_{method_1}", f"_{method_2}"), how="inner")
    df_out = df_out[(df_out[f"ratio_{method_1}"] > 0) & (df_out[f"ratio_{method_2}"] > 0)]
    min_trait = min(df_out[f"ratio_{method_1}"].min(), df_out[f"ratio_{method_2}"].min())
    max_trait = max(df_out[f"ratio_{method_1}"].max(), df_out[f"ratio_{method_2}"].max())
    datasets = sorted(set(df_out["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4))
    for x_1, (dataset, df_gr) in enumerate(df_out.groupby("dataset")):
        ax = axs[x_1] if len(datasets) > 1 else axs
        ax.set_xlabel(f"Ratio {method_1}")
        ax.set_ylabel(f"Ratio {method_2}")
        ax.plot([1, 1], [min_trait, max_trait], "--", color="black", alpha=0.5)
        ax.plot([min_trait, max_trait], [1, 1], "--", color="black", alpha=0.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim((min_trait, max_trait))
        ax.set_ylim((min_trait, max_trait))
        if len(df_gr) < 2:
            continue
        corr = df_gr[f"ratio_{method_1}"].corr(df_gr[f"ratio_{method_2}"])
        ax.plot(df_gr[f"ratio_{method_1}"], df_gr[f"ratio_{method_2}"], "o", alpha=0.5, label=f"r={corr:.2g}")
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
