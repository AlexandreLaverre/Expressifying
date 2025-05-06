import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input: str, gene_table: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    omega_list = ["ωA_phy", "ω", "ω0"]
    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input, sep='\t')
    df_sequence = pd.read_csv(gene_table, sep='\t')
    df_sequence["trait"] = df_sequence["ENSG"].apply(lambda x: x.split("_")[0])
    omega_list = [col for col in omega_list if col in df_sequence.columns]
    assert len(omega_list) > 0, f"No omega found in {gene_table}"
    df_join = pd.merge(df_trait, df_sequence, on="trait", suffixes=("_trait", "_sequence"), how="inner")
    assert len(df_join) > 0, f"No overlap between {tsv_input} and {gene_table}"
    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=len(omega_list), ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4 * len(omega_list)))
    for x_1, (data, df_gr) in enumerate(df_join.groupby("dataset")):
        df_gr["ratioqcut"] = pd.qcut(df_gr["ratio"], q=50)
        for x_2, omega in enumerate(omega_list):
            ax = axs[x_2, x_1] if len(datasets) > 1 else axs[x_2]
            ax.set_xlabel("ratio")
            ax.set_ylabel(omega)
            ax.set_title(data)
            ax.axvline(1, linestyle="--", color="black", alpha=0.5)
            ax.set_xscale("log")
            if omega == "ωA_phy":
                ax.axhline(0, linestyle="--", color="black", alpha=0.5)
            else:
                ax.set_yscale("log")

            df = df_gr.groupby("ratioqcut", observed=False).agg({"ratio": "mean", omega: "mean"}).reset_index()
            corr = df[f"ratio"].corr(df[omega])
            ax.plot(df[f"ratio"], df[omega], "o", alpha=0.5, label=f"r={corr:.2g}")
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv file", required=False)
    parser.add_argument("--gene_table", help="Input gene table file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.gene_table, args.output)
