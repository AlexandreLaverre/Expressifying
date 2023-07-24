import os
import argparse
import pandas as pd
from neutrality_index import open_tree, prune_tree


def main(input_tree, input_traits, input_var_within, neutrality_index, output_tree, output_traits, output_var_within):
    for path in [input_tree, input_traits, input_var_within, neutrality_index]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_tree, output_traits, output_var_within]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    tree = open_tree(input_tree, format_ete3=1)
    df_traits = pd.read_csv(input_traits, sep="\t")
    df_var_within = pd.read_csv(input_var_within, sep="\t")
    df_neutrality_index = pd.read_csv(neutrality_index, sep="\t")
    greater_than_one = (df_neutrality_index["ratio"] > 1.0)
    if greater_than_one.sum() == 0:
        df_neutrality_index = df_neutrality_index.sort_values(by="ratio", ascending=False)
        # Keep the 3 first traits
        df_neutrality_index = df_neutrality_index.iloc[:3]
    else:
        df_neutrality_index = df_neutrality_index[df_neutrality_index["ratio"] > 1.0]
    ortho_list = df_neutrality_index["trait"].tolist()

    set_taxa_names = set(tree.get_leaf_names())
    col_traits = ["TaxonName"] + [f"{i}_mean" for i in ortho_list]
    df_traits = df_traits[col_traits]
    df_traits = df_traits.dropna(subset=[f"{i}_mean" for i in ortho_list], how='all')
    columns = ["TaxonName", "Nucleotide_diversity"] + [f"{i}_variance" for i in ortho_list]
    df_var_within = df_var_within[columns]
    df_var_within = df_var_within.dropna(subset=[f"{i}_variance" for i in ortho_list], how='all')
    set_taxa_names = set_taxa_names.intersection(set(df_traits["TaxonName"].tolist()))
    tree = prune_tree(tree, list(set_taxa_names))

    df_traits.to_csv(output_traits, sep="\t", index=False, na_rep="NaN")
    df_var_within.to_csv(output_var_within, sep="\t", index=False, na_rep="NaN")
    tree.write(outfile=output_tree, format=3)
    print(f"The tree has {len(set_taxa_names)} leaves")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tree", help="Input tree file", required=True)
    parser.add_argument("--input_traits", help="Input traits file", required=True)
    parser.add_argument('--input_var_within', help="Input var_within file", required=True)
    parser.add_argument('--neutrality_index', help="Input neutrality index file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    parser.add_argument('--output_var_within', help="Output var_within file", required=True)
    args = parser.parse_args()
    main(args.input_tree, args.input_traits, args.input_var_within, args.neutrality_index, args.output_tree,
         args.output_traits, args.output_var_within)
