import os
import argparse
from libraries_plot import *
from trait_cross_corr import logit, categorize_p_value


def main(tsv_input_list: str, gene_table: str, cross_table: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)

    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    clades_dico = {}
    for tsv_input in tsv_input_list:
        clade = os.path.basename(os.path.dirname(tsv_input)).split("_")[1]
        clades_dico[clade] = pd.read_csv(tsv_input, sep='\t')
    clade_1 = "primates"
    clade_2 = "NoPrimates"
    df_out = pd.merge(clades_dico[clade_1], clades_dico[clade_2], on=["trait", "dataset"],
                      suffixes=(f"_{clade_1}", f"_{clade_2}"), how="inner")

    df_gene = pd.read_csv(gene_table, sep='\t')
    dico_gene = {row["ENSG"].split('_')[1].upper(): row["ENSG"].split('_')[0] for _, row in df_gene.iterrows()}

    df_cross = pd.read_csv(cross_table, sep=',')
    df_cross["trait"] = df_cross["gene"].apply(lambda x: str(x).upper())
    df_cross["trait"] = df_cross["trait"].apply(lambda x: dico_gene[x] if x in dico_gene else x)
    # groupby trait and take the mean
    df_cross = df_cross.groupby("trait").agg({"p-value": "mean", "parental difference (B6-DBA)": "mean"}).reset_index()
    df_cross["logit p-value"] = logit(df_cross["p-value"])
    df_cross["abs difference"] = np.abs(df_cross["parental difference (B6-DBA)"])

    df_cross["category"] = df_cross["p-value"].apply(categorize_p_value)
    value_list = ["logit p-value", "abs difference"]

    gpby = {dataset: df_gr for dataset, df_gr in df_out.groupby("dataset")}
    gb_join = {d: pd.merge(gpby[d], df_cross, on="trait", suffixes=("_trait", "_pbs"), how="inner") for d in gpby}
    datasets = sorted(gb_join.keys())
    fig, axs = plt.subplots(nrows=len(value_list), ncols=len(datasets),
                            figsize=(5 * len(datasets), 4 * len(value_list)))

    for x_2, dataset in enumerate(datasets):
        df_join = gb_join[dataset]
        # Median for which half the genes are above and half are below
        df_join["category"] = annotate_category(df_join, clade_1, clade_2, "ratio")
        for x_1, value in enumerate(value_list):
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs[x_1]
            plot_cat_box(ax, x_label="category", y_label=value, df=df_join)
            ax.set_title(dataset)
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv files", nargs="+", required=True)
    parser.add_argument("--gene_table", help="Input gene table file", required=False)
    parser.add_argument("--cross_table", help="Input cross table file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.gene_table, args.cross_table, args.output)
