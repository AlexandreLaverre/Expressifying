import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    df_out = pd.read_csv(tsv_input, sep='\t')
    datasets = sorted(set(df_out["dataset"]))
    gp = {group: df_group for group, df_group in df_out.groupby("dataset")}
    fig, axs = plt.subplots(nrows=len(datasets), ncols=len(datasets), sharex=True, sharey=True,
                            figsize=(len(datasets) * 5, len(datasets) * 4))
    for x_1, data_1 in enumerate(datasets):
        for x_2, data_2 in enumerate(datasets):
            df_1 = gp[data_1]
            df_2 = gp[data_2]
            joint_df = pd.merge(df_1, df_2, on="trait", suffixes=("_1", "_2"), how="inner")
            joint_df["qcut"] = pd.qcut(joint_df["ratio_1"], q=50)
            cols = ["ratio_1", "ratio_2"]
            if "pp_ratio_greater_1" in df_1.columns and "pp_ratio_greater_1" in df_2.columns:
                cols += ["pp_ratio_greater_1_1", "pp_ratio_greater_1_2"]
            joint_df = joint_df.groupby("qcut", observed=False).agg({col: "mean" for col in cols}).reset_index()
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs
            ax.set_xlabel(f"{data_1} ratio")
            ax.set_ylabel(f"{data_2} ratio")
            # vertical line and horizontal line for 1
            ax.axvline(1, linestyle="--", color="black", alpha=0.5)
            ax.axhline(1, linestyle="--", color="black", alpha=0.5)
            ax.set_xscale("log")
            ax.set_yscale("log")
            if len(joint_df) < 2:
                continue
            corr = joint_df["ratio_1"].corr(joint_df["ratio_2"])
            pp, pp_val, pp_label = ("ratio", 1.0, "ρ")
            if "pp_ratio_greater_1" in df_1.columns and "pp_ratio_greater_1" in df_2.columns:
                pp, pp_val, pp_label = ("pp_ratio_greater_1", 0.95, "pp")

            for i in [1, 2]:
                joint_df[f"pos_{i}"] = (joint_df[f"{pp}_{i}"] > pp_val)
                joint_df[f"neg_{i}"] = ~joint_df[f"pos_{i}"]
            for p1, p2, c, l in [("pos_1", "pos_2", "red", f"Both with {pp_label} > {pp_val}"),
                                 ("pos_1", "neg_2", "blue", f"{data_1} with {pp_label} > {pp_val}"),
                                 ("neg_1", "pos_2", "green", f"{data_2} with {pp_label} > {pp_val}"),
                                 ("neg_1", "neg_2", "orange", f"Both with {pp_label} < {pp_val}")]:
                pp1_pp2 = joint_df[p1] & joint_df[p2]
                if sum(pp1_pp2) == 0:
                    continue
                label = f"{l} (n={sum(pp1_pp2)})".capitalize()
                ax.scatter(joint_df["ratio_1"][pp1_pp2], joint_df["ratio_2"][pp1_pp2], c=c, label=label)
                ax.set_title(f"{len(joint_df)} bins ($R^2$={corr * corr:.2f})")
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv file", required=False)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
