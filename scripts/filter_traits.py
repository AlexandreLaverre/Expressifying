import os
import argparse
import pandas as pd
from neutrality_index import open_tree, prune_tree
from pre_processed_traits import name_internal_nodes


def main(input_tree, input_traits, input_var_within, neutrality_index, output_dir):
    for path in [input_tree, input_traits, input_var_within, neutrality_index]:
        assert os.path.exists(path), f"Path {path} does not exist"
    os.makedirs(output_dir, exist_ok=True)

    tree = open_tree(input_tree, format_ete3=1)
    df_traits = pd.read_csv(input_traits, sep="\t")
    df_var_within = pd.read_csv(input_var_within, sep="\t")
    df_neutrality_index = pd.read_csv(neutrality_index, sep="\t")
    df_neutrality_index = df_neutrality_index.sort_values(by="ratio", ascending=False)
    enough_species = ((df_neutrality_index["nbr_taxa_between"] >= 5) & (df_neutrality_index["nbr_taxa_within"] >= 5))
    print(f"Keeping {sum(enough_species)} traits with more than 5 species out of {len(df_neutrality_index)}")
    df_neutrality_index = df_neutrality_index[enough_species]

    ortho_list = df_neutrality_index["trait"].tolist()
    set_taxa_names = set(tree.get_leaf_names())
    for gene in ortho_list:
        col_traits = ["TaxonName", f"{gene}_mean"]
        df_gene = df_traits[col_traits].copy()
        df_gene = df_gene.dropna(subset=f"{gene}_mean", how='all')
        columns = ["TaxonName", "Nucleotide_diversity", f"{gene}_variance", f"{gene}_heritability_lower",
                   f"{gene}_heritability_upper"]
        df_gene_within = df_var_within[columns].copy()
        df_gene_within = df_gene_within.dropna(subset=columns, how='all')
        gene_taxa_names = set_taxa_names.intersection(set(df_gene["TaxonName"].tolist()))
        if len(gene_taxa_names) < 5:
            print(f"Skipping {gene} because it has less than 10 taxa")
            continue
        gene_tree = prune_tree(tree, list(gene_taxa_names))

        df_gene.to_csv(f"{output_dir}/{gene}.traits.tsv", sep="\t", index=False, na_rep="NaN")
        df_gene_within.to_csv(f"{output_dir}/{gene}.var_within.tsv", sep="\t", index=False, na_rep="NaN")
        gene_tree = name_internal_nodes(gene_tree)
        gene_tree.write(outfile=f"{output_dir}/{gene}.tree", format=3)
        print(f"The tree has {len(set_taxa_names)} leaves")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tree", help="Input tree file", required=True)
    parser.add_argument("--input_traits", help="Input traits file", required=True)
    parser.add_argument('--input_var_within', help="Input var_within file", required=True)
    parser.add_argument('--neutrality_index', help="Input neutrality index file", required=True)
    parser.add_argument("--output_dir", help="Output directory", required=True)
    args = parser.parse_args()
    main(args.input_tree, args.input_traits, args.input_var_within, args.neutrality_index, args.output_dir)
