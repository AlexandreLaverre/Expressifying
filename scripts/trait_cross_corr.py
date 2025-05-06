import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def logit(x):
    # Bound to -4, 4
    x = np.clip(x, 1e-3, 1 - 1e-3)
    return np.log(x / (1 - x))


def main(tsv_input: str, gene_table: str, cross_table: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    value_list = ["logit p-value", "parental difference (B6-DBA)"]
    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input, sep='\t')

    df_gene = pd.read_csv(gene_table, sep='\t')
    dico_gene = {row["ENSG"].split('_')[1].upper(): row["ENSG"].split('_')[0] for _, row in df_gene.iterrows()}

    df_cross = pd.read_csv(cross_table, sep=',')
    df_cross["trait"] = df_cross["gene"].apply(lambda x: str(x).upper())
    df_cross["trait"] = df_cross["trait"].apply(lambda x: dico_gene[x] if x in dico_gene else x)
    df_cross["logit p-value"] = logit(df_cross["nucleus accumbens p-value"])

    value_list = [col for col in value_list if col in df_cross.columns]
    assert len(value_list) > 0, f"No value found in {gene_table}"
    df_join = pd.merge(df_trait, df_cross, on="trait", suffixes=("_trait", "_cross"), how="inner")
    assert len(df_join) > 0, f"No overlap between {tsv_input} and {gene_table}"

    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=len(value_list), ncols=len(datasets), sharex='all', sharey='row',
                            figsize=(5 * len(datasets), 4 * len(value_list)))
    for x_1, (dataset, df_gr) in enumerate(df_join.groupby("dataset")):
        df_gr["ratioqcut"] = pd.qcut(df_gr["ratio"], q=50)
        for x_2, value in enumerate(value_list):
            ax = axs[x_2, x_1] if len(datasets) > 1 else axs[x_2]
            ax.set_title(dataset)
            ax.set_xlabel("ratio")
            ax.set_ylabel(value)

            df = df_gr.groupby("ratioqcut", observed=False).agg({"ratio": "mean", value: "mean"}).reset_index()
            ax.axvline(1, linestyle="--", color="black", alpha=0.5)
            ax.set_xscale("log")

            corr = df[f"ratio"].corr(df[value])
            ax.plot(df[f"ratio"], df[value], "o", alpha=0.5, label=f"r={corr:.2g}")
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv file", required=False)
    parser.add_argument("--gene_table", help="Input gene table file", required=False)
    parser.add_argument("--cross_table", help="Input cross table file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.gene_table, args.cross_table, args.output)
