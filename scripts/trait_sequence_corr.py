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
    min_trait, max_trait = df_trait["ratio"].min(), df_trait["ratio"].max()
    minmax_omega_dict = {omega: (df_join[omega].min(), df_join[omega].max()) for omega in omega_list}
    groups = df_join.groupby("dataset")
    group_by = {key1: group for key1, group in groups}
    datasets = sorted(group_by.keys())
    fig, axs = plt.subplots(nrows=len(omega_list), ncols=len(datasets),
                            figsize=(5 * len(datasets), 4 * len(omega_list)), sharex='all', sharey='row')
    for x_1, data_1 in enumerate(datasets):
        df_gr = group_by[data_1]
        for x_2, omega in enumerate(omega_list):
            ax = axs[x_2, x_1] if len(datasets) > 1 else axs[x_2]
            ax.set_xlabel(f"{data_1} ratio")
            ax.set_ylabel(omega)
            ax.plot([1, 1], [minmax_omega_dict[omega][0], minmax_omega_dict[omega]][1], "--", color="black", alpha=0.5)
            ax.set_xscale("log")
            if omega == "ωA_phy":
                ax.plot([min_trait, max_trait], [0, 0], "--", color="black", alpha=0.5)
            else:
                ax.set_yscale("log")
            if len(df_gr) < 2:
                continue
            corr = df_gr[f"ratio"].corr(df_gr[omega])
            ax.plot(df_gr[f"ratio"], df_gr[omega], "o", alpha=0.5, label=f"r={corr:.2g}")
            ax.set_ylim(minmax_omega_dict[omega])
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
