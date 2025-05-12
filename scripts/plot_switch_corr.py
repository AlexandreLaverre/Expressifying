import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input_ratio: str, tsv_input_switch: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input_ratio, sep='\t')
    
    df_switch = pd.read_csv(tsv_input_switch, sep='\t')
    # Renamed the column "simu" to "trait" in df_switch
    df_switch.rename(columns={"simu": "trait", "dataset_switch": "dataset"}, inplace=True)
    df_switch["dataset"] = df_switch["dataset"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_join = pd.merge(df_trait, df_switch, on=["trait", "dataset"], suffixes=("_trait", "_switch"), how="inner")
    assert len(df_join) > 0, f"No overlap between {tsv_input_ratio} and {tsv_input_switch}"
    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4))
    for x_1, (data, df_gr) in enumerate(df_join.groupby("dataset")):
        df_gr["ratioqcut"] = pd.qcut(df_gr["ratio"], q=50)

        ax = axs[x_1] if len(datasets) > 1 else axs
        ax.set_xlabel("ratio")
        ax.set_ylabel("switch")
        ax.set_title(data)
        ax.axvline(1, linestyle="--", color="black", alpha=0.5)
        ax.set_xscale("log")
        if len(df_gr) < 2:
            continue

        df = df_gr.groupby("ratioqcut", observed=False).agg({"ratio": "mean", "is_nuc": "mean"}).reset_index()
        corr = df[f"ratio"].corr(df["is_nuc"])
        ax.plot(df[f"ratio"], df["is_nuc"], "o", alpha=0.5, label=f"r={corr:.2g}")
        ax.legend()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_ratio", help="Input ratio tsv file", required=False)
    parser.add_argument("--input_switch", help="Input switch tsv file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input_ratio, args.input_switch, args.output)
