import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main(tsv_input: str, output_file: str):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df_out = pd.read_csv(tsv_input, sep='\t')
    datasets = sorted(set(df_out["dataset"]))
    gp = {group: df_group for group, df_group in df_out.groupby("dataset")}
    fig, axs = plt.subplots(nrows=len(datasets), ncols=len(datasets),
                            figsize=(len(datasets) * 6, len(datasets) * 4))
    for x_1, data_1 in enumerate(datasets):
        for x_2, data_2 in enumerate(datasets):
            df_1 = gp[data_1]
            df_2 = gp[data_2]
            joint_df = pd.merge(df_1, df_2, on="trait", suffixes=("_1", "_2"), how="inner")
            joint_df["qcut"] = pd.qcut(joint_df["ratio_1"], q=min(35, len(joint_df) // 2), duplicates="drop")
            cols = ["ratio_1", "ratio_2"]
            if "pp_ratio_greater_1" in df_1.columns and "pp_ratio_greater_1" in df_2.columns:
                cols += ["pp_ratio_greater_1_1", "pp_ratio_greater_1_2"]
            # joint_df = joint_df.groupby("qcut", observed=False).agg({col: "mean" for col in cols}).reset_index()
            ax = axs[x_1, x_2] if len(datasets) > 1 else axs
            n_1 = data_1.split("_")[0].title()
            n_2 = data_2.split("_")[0].title()
            ax.set_xlabel(f"{n_1} rate")
            ax.set_ylabel(f"{n_2} rate")
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
                pp, pp_val, pp_label = ("pp_ratio_greater_1", 0.5, "pp")

            for i in [1, 2]:
                joint_df[f"pos_{i}"] = (joint_df[f"{pp}_{i}"] > pp_val)
                joint_df[f"neg_{i}"] = ~joint_df[f"pos_{i}"]
            for p1, p2, c, l in [("pos_1", "pos_2", "red", f"Positive selection in both organs"),
                                 ("pos_1", "neg_2", "blue", f"Positive selection in {n_1}"),
                                 ("neg_1", "pos_2", "green", f"Positive selection in {n_2}"),
                                 ("neg_1", "neg_2", "orange", f"No evidence of positive selection")]:
                pp1_pp2 = joint_df[p1] & joint_df[p2]
                if sum(pp1_pp2) == 0:
                    continue
                label = f"{l} (n={sum(pp1_pp2)})".capitalize()
                ax.scatter(joint_df["ratio_1"][pp1_pp2], joint_df["ratio_2"][pp1_pp2], c=c, label=label, alpha=0.25)
                ax.set_title(f"$R^2$={corr * corr:.2f} (n={len(joint_df)})")
            ax.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close("all")
    plt.clf()

    violin_output = ".".join(output_file.split(".")[:-1]) + "_violin.pdf"
    fig = plt.figure(figsize=(2 + len(datasets), 5))
    ax = fig.add_subplot(111)
    # Prepare data for violinplot
    df = pd.DataFrame({"y": np.concatenate([gp[ds]["ratio"] for ds in datasets]),
                       "x": np.concatenate([[ds] * len(gp[ds]["ratio"]) for ds in datasets])})
    df = df[np.isfinite(df["y"]) & df["y"] > 0.0]
    sns.violinplot(data=df, x="x", y="y", ax=ax, cut=0, log_scale=True, legend=False, linewidth=0.1, linecolor="auto")
    # Add the number of ratio above 1 for each violin
    legend_labels = []
    for i, ds in enumerate(datasets):
        mean_txt = f"rate > 1\nn=${sum(gp[ds]["ratio"] > 1.0)}$"
        ax.text(i, 7.0, mean_txt,
                bbox=dict(facecolor="white", edgecolor="#EC6231", boxstyle="round", pad=0.15),
                fontsize=11, ha="center", va="center", zorder=10)
        # Add the mean of the ratio for each violin
        mean_v = gp[ds]["ratio"].mean()
        ax.text(i, mean_v, f"mean\n${mean_v:.2f}$",
                bbox=dict(facecolor="white", edgecolor="black", boxstyle="round", pad=0.15),
                fontsize=11, ha="center", va="center", zorder=10)
        legend_labels.append(f"{ds.split("_")[0].title()}\nn={len(gp[ds]["ratio"])}")
    ax.axhspan(1.0, np.max(df["y"]), color="#EC6231", alpha=0.25, zorder=-1)
    ax.axhspan(np.min(df["y"]), 1.0, color="#8FB03E", alpha=0.25, zorder=-1)
    ax.set_xlabel("")
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels(legend_labels, rotation=0, ha='center')
    ax.set_ylabel("Rate (log scale)")
    ax.axhline(1.0, c="black")
    ax.set_ylim(np.min(df["y"]), np.max(df["y"]))
    plt.tight_layout()
    plt.savefig(violin_output)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tsv file", required=False)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
