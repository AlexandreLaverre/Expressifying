import os
import argparse
from libraries_plot import *

def main(tsv_input_ratio: str, tsv_input_switch: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input_ratio, sep='\t')
    
    df_switch = pd.read_csv(tsv_input_switch, sep='\t')
    # Renamed the column "simu" to "trait" in df_switch
    df_switch.rename(columns={"simu": "trait", "dataset_switch": "dataset"}, inplace=True)
    df_switch["dataset"] = df_switch["dataset"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_join = pd.merge(df_trait, df_switch, on=["trait", "dataset"], suffixes=("_trait", "_switch"), how="inner")
    if len(df_join) == 0:
        print(f"No overlap between {tsv_input_ratio} and {tsv_input_switch}")
        plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
        plt.savefig(output_pdf)
        return
    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), sharex='col', sharey='row',
                            figsize=(5 * len(datasets), 4))

    for x_1, (data, df_gr) in enumerate(df_join.groupby("dataset")):
        ax = axs[x_1] if len(datasets) > 1 else axs
        ax.set_xlabel("Ratio (log10)")
        ax.set_ylabel("switch")
        ax.set_title(f"{data} (n={len(df_gr)})")
        ax.axvline(0.0, linestyle="--", color="black", alpha=0.5)
        df_gr["log_ratio"] = df_gr["ratio"].apply(lambda x: np.log10(x) if x > 0 else np.nan)
        plot_scatter_bins(ax, x_label="log_ratio", y_label="is_nuc", df=df_gr, q=200)
        # plot_2d_histogram(ax, df_gr["log_ratio"], df_gr["is_nuc"], bins=30, cmap="Blues")

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
