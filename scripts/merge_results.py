import os
import argparse
import pandas as pd


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(tsv_ML_list: str, tsv_Bayes_list: str, output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    list_df = []
    for method, tsv_list in [("ML", tsv_ML_list), ("Bayesian", tsv_Bayes_list)]:
        if tsv_list is None:
            continue
        for path in tsv_list:
            df = pd.read_csv(path, sep='\t')
            df["method"] = method
            df["dataset"] = os.path.basename(path).replace(".tsv", "").split("_")[0]
            list_df.append(df)
    df_out = pd.concat(list_df)
    # Sort by trait, dataset, sex, logT and then method
    df_out = df_out.sort_values(by=["ratio"], ascending=False)
    df_out.to_csv(output, sep="\t", index=False, float_format="%.3f")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv_ML", nargs="+", help="Input tsv ML file", required=False)
    parser.add_argument("--tsv_Bayes", nargs="+", help="Input tsv BayesCode file", required=False)
    parser.add_argument("--output", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.tsv_ML, args.tsv_Bayes, args.output)
