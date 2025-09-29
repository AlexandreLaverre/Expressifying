import os
import argparse
from libraries_plot import *


def main(tsv_input: str, gene_table: str, cross_table: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input, sep='\t')

    df_gene = pd.read_csv(gene_table, sep='\t')
    dico_gene = {row["ENSG"].split('_')[1].upper(): row["ENSG"].split('_')[0] for _, row in df_gene.iterrows()}

    df_cross = pd.read_csv(cross_table, sep=',')
    df_cross["trait"] = df_cross["Gene.Symbol"].apply(lambda x: str(x).upper())
    df_cross["trait"] = df_cross["trait"].apply(lambda x: dico_gene[x] if x in dico_gene else x)
    # groupby trait and take the mean

    df_join = pd.merge(df_trait, df_cross, on="trait", suffixes=("_trait", "_cross"), how="inner")
    # Filter only to M3 and M4 datasets
    df_join = df_join[df_join["Table"].str.contains("M3|M4")]
    # Replace M3 by M3: Dir and M4 by M4: Stab
    df_join["Table"] = df_join["Table"].replace({"M3": "M2 (Directional)", "M4": "M1: (Stabilizing)"})
    assert len(df_join) > 0, f"No overlap between {tsv_input} and {gene_table}"

    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=1, ncols=len(datasets), figsize=(5 * len(datasets), 4))
    gpby = df_join.groupby("dataset")
    for x_1, (dataset, df_gr) in enumerate(gpby):
        df_gr["log_ratio"] = df_gr["ratio"].apply(lambda x: np.log10(x) if x > 0 else np.nan)
        for x_2, value in enumerate(["Table"]):
            ax = axs[x_1] if len(datasets) > 1 else axs
            ax.set_title(f"{str(dataset).split('_')[0].capitalize()} (n={len(df_gr)})")
            ax.set_ylabel("Ratio")
            plot_cat_box(ax, x_label=value, y_label="ratio", df=df_gr, scale="log")

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
