import os
import argparse
import pandas as pd


def main(tsv_list: str, output: str):
    output_dir = os.path.dirname(output)
    os.makedirs(output_dir, exist_ok=True)
    list_df = []
    for path in tsv_list:
        if not os.path.exists(path):
            continue
        if open(path).read().strip() == "":
            continue
        df = pd.read_csv(path, sep='\t')
        df["dataset"] = "_".join(os.path.basename(path).split("_")[:-1])
        list_df.append(df)
    df_out = pd.concat(list_df)
    sort_key = [i for i in ["ratio_pv", "ratio", "is_nuc"] if i in df_out.columns][0]
    df_out = df_out.sort_values(by=[sort_key], ascending=False)
    df_out.to_csv(output, sep="\t", index=False, float_format="%.4f")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv_results", nargs="+", help="Input tsv trait files", required=False)
    parser.add_argument("--output", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.tsv_results, args.output)
