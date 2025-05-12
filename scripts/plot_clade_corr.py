import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input: list, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    assert len(tsv_input) == 2, f"Two input files are required, got {len(tsv_input)}"
    clade_1 = os.path.basename(os.path.dirname(tsv_input[0])).split("_")[1]
    clade_2 = os.path.basename(os.path.dirname(tsv_input[1])).split("_")[1]
    assert clade_1 != clade_2, f"Methods should be different, got {clade_1} and {clade_2}"
    df_1 = pd.read_csv(tsv_input[0], sep='\t')
    df_2 = pd.read_csv(tsv_input[1], sep='\t')
    df_out = pd.merge(df_1, df_2, on=["trait", "dataset"], suffixes=(f"_{clade_1}", f"_{clade_2}"), how="inner")
    datasets = sorted(set(df_out["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4))
    for x_1, (dataset, df_gr) in enumerate(df_out.groupby("dataset")):
        ax = axs[x_1] if len(datasets) > 1 else axs
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
