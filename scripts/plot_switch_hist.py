import os
import argparse
from libraries_plot import *

def main(tsv_input_switch: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    
    df_switch = pd.read_csv(tsv_input_switch, sep='\t')
    # Renamed the column "simu" to "trait" in df_switch
    df_switch.rename(columns={"simu": "trait", "dataset_switch": "dataset"}, inplace=True)
    df_switch["dataset"] = df_switch["dataset"].apply(lambda x: "_".join(x.split("_")[1:]))

    datasets = sorted(set(df_switch["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), sharex='col', sharey='row',
                            figsize=(5 * len(datasets), 4))

    for x_1, (data, df_gr) in enumerate(df_switch.groupby("dataset")):
        ax = axs[x_1] if len(datasets) > 1 else axs
        ax.set_xlabel("Switch (phylogram favored)")
        ax.set_ylabel("frequency")
        ax.set_title(f"{data} (n={len(df_gr)})")
        # Plot the histogram of switch frequencies
        ax.hist(df_gr["is_nuc"], bins=30, color=BLUE, alpha=1.0, edgecolor='black')
        ax.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_switch", help="Input switch tsv file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input_switch, args.output)
