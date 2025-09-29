import os
import argparse
from libraries_plot import *


def main(tsv_input: str, fst_table: str, pst_table: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)

    os.makedirs(output_dir, exist_ok=True)
    df_trait = pd.read_csv(tsv_input, sep='\t')

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
    df_join = pd.merge(df_trait, df, on="trait", suffixes=("_trait", "_pbs"), how="inner")
    assert len(df_join) > 0, f"No overlap between {tsv_input} and {fst_table}"

    datasets = sorted(set(df_join["dataset"]))
    fig, axs = plt.subplots(nrows=len(value_list), ncols=len(datasets), sharex='col', sharey='row',
                            figsize=(5 * len(datasets), 4 * len(value_list)))
    gpby = df_join.groupby("dataset")
    n = min([len(df) for _, df in gpby])
    for x_1, (dataset, df_gr) in enumerate(df_join.groupby("dataset")):
        # subsample the dataframe with 1000 rows
        #if len(df_gr) > n:
        #    df_gr = df_gr.sample(n, random_state=42)
        df_gr["log_ratio"] = df_gr["ratio"].apply(lambda x: np.log10(x) if x > 0 else np.nan)
        df_gr = df_gr.dropna(subset="log_ratio")
        df_gr["qcut"] = pd.qcut(df_gr["log_ratio"], q=4, duplicates="drop", labels=["Q1", "Q2", "Q3", "Q4"])
        for x_2, value in enumerate(value_list):
            ax = axs[x_2, x_1] if len(datasets) > 1 else axs[x_2]
            ax.set_title(f"{dataset} (n={len(df_gr)})")
            ax.set_xlabel("Ratio (log10)")
            ax.set_ylabel(value)
            # ax.axvline(0.0, linestyle="--", color="black", alpha=0.5)
            plot_cat_box(ax, x_label="qcut", y_label=value, df=df_gr)
            # plot_2d_histogram(ax, df_gr["log_ratio"], df_gr[value], bins=30, cmap="Blues", cmin=1)

    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv file", required=False)
    parser.add_argument("--fst_table", help="Input fst table file", required=False)
    parser.add_argument("--pst_table", help="Input pbs table file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.fst_table, args.pst_table, args.output)
