import os
import argparse
import pandas as pd


def main(input_path: str, output_path: str):
    out_list_df = []
    if os.path.exists(input_path):
        for file_path in os.listdir(input_path):
            if not file_path.endswith(".tsv"):
                continue
            df = pd.read_csv(f"{input_path}/{file_path}", sep="\t")
            out_list_df.append(df)
        df_out = pd.concat(out_list_df, axis=0)
        df_out.to_csv(output_path, sep="\t", index=False)
    else:
        print(f"Input path {input_path} does not exist.")
        open(output_path, "w").close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input traces directory", required=False)
    parser.add_argument("--output", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
