import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main(tsv_input: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    df_out = pd.read_csv(tsv_input, sep='\t')
    min_trait = df_out["ratio"].min()
    max_trait = df_out["ratio"].max()
    datasets = sorted(set(df_out["dataset"]))
    fig, axs = plt.subplots(nrows=len(datasets), ncols=len(datasets), figsize=(16, 16), sharex=True, sharey=True)
    for x_1, data_1 in enumerate(datasets):
        for x_2, data_2 in enumerate(datasets):
            df_1 = df_out[df_out["dataset"] == data_1]
            df_2 = df_out[df_out["dataset"] == data_2]
            joint_df = pd.merge(df_1, df_2, on="trait", suffixes=("_1", "_2"), how="inner")
            ax = axs[x_1, x_2]
            ax.set_xlabel(f"{data_1} ratio")
            ax.set_ylabel(f"{data_2} ratio")
            # vertical line and horizontal line for 1
            ax.plot((1, 1), (min_trait, max_trait), linestyle="--", color="black")
            ax.plot((min_trait, max_trait), (1, 1), linestyle="--", color="black")
            ax.set_xscale("log")
            ax.set_yscale("log")
            if len(joint_df) < 2:
                continue
            corr = joint_df["ratio_1"].corr(joint_df["ratio_2"])
            pp, pp_val, pp_label = ("ratio", 1.0, "ρ")
            if "pp_ratio_greater_1" in joint_df.columns:
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
                ax.set_title(f"{len(joint_df)} genes ($R^2$={corr * corr:.2f})")
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
