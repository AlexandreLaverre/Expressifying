import os
import argparse
from libraries_plot import *


def main(tsv_input_list: str, fst_table: str, pst_table: str, output_pdf: str):
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

    df_fst = pd.read_csv(fst_table, sep=',')
    df_pst = pd.read_csv(pst_table, sep=',')
    df = pd.merge(df_fst, df_pst, on="Ensembl ID", suffixes=("_fst", "_pst"), how="inner")
    df["trait"] = df["Ensembl ID"].apply(lambda x: str(x).upper())
    pops = ["TSI", "GBR", "FIN", "YRI"]
    v = "pst"
    value_list = [f"{pop}_{v}" for pop in pops]
    df[f"max_{v}"] = df[value_list].max(axis=1)
    df[f"mean_{v}"] = df[value_list].mean(axis=1)
    df["r"] = df[f"mean_{v}"] / df[[f"{pop}_fst" for pop in pops]].mean(axis=1)
    value_list += [f"max_{v}", f"mean_{v}", "r"]
    assert len(value_list) > 0, f"No value found in {fst_table}"

    gpby = {dataset: df_gr for dataset, df_gr in df_out.groupby("dataset")}
    gb_join = {d: pd.merge(gpby[d], df, on="trait", suffixes=("_trait", "_pbs"), how="inner") for d in gpby}
    datasets = sorted(gb_join.keys())
    fig, axs = plt.subplots(nrows=len(value_list), ncols=len(datasets),
                            figsize=(5 * len(datasets), 4 * len(value_list)))

    for x_2, dataset in enumerate(datasets):
        df_join = gb_join[dataset]
        df_join["category"] = annotate_category(df_join, clade_1, clade_2, "ratio")
        for x_1, value in enumerate(value_list):
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs[x_1]
            assert len(df_join) > 0, f"No overlap between dataset and {fst_table}"
            plot_cat_box(ax, x_label="category", y_label=value, df=df_join)
            ax.set_title(dataset)
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv files", nargs="+", required=True)
    parser.add_argument("--fst_table", help="Input fst table file", required=False)
    parser.add_argument("--pst_table", help="Input pbs table file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.fst_table, args.pst_table, args.output)
