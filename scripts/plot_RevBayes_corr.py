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
    df_switch["dataset"] = df_switch["dataset"].apply(lambda x: x.replace("Switchnodes_", ""))
    df_join = pd.merge(df_trait, df_switch, on=["trait", "dataset"], suffixes=("_trait", "_switch"), how="inner")
    assert len(df_join) > 0, f"No overlap between {tsv_input_ratio} and {tsv_input_switch}"
    
    min_trait, max_trait = df_trait["ratio"].min(), df_trait["ratio"].max()
    groups = df_join.groupby("dataset")
    group_by = {key1: group for key1, group in groups}
    datasets = sorted(group_by.keys())
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets),
                            figsize=(5 * len(datasets), 4), sharex='all', sharey='row')
    for x_1, data_1 in enumerate(datasets):
        df_gr = group_by[data_1]
        ax = axs[x_1] if len(datasets) > 1 else axs
        ax.set_xlabel(f"{data_1} ratio")
        ax.set_ylabel("switch")
        # vertical line
        ax.plot([1, 1], [0.0, 1.0], "--", color="black", alpha=0.5)
        ax.set_xscale("log")
        ax.plot([min_trait, max_trait], [0.5, 0.5], "--", color="black", alpha=0.5)
        if len(df_gr) < 2:
            continue
        corr = df_gr[f"ratio"].corr(df_gr["is_nuc"])
        ax.plot(df_gr[f"ratio"], df_gr["is_nuc"], "o", alpha=0.5, label=f"r={corr:.2g}")
        ax.set_ylim((0.0, 1.0))
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
