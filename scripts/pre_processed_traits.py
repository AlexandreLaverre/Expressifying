import os
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from ete3 import Tree
from neutrality_index import open_tree, prune_tree


def rename_tree(tree: Tree) -> Tree:
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            if node.name.endswith("_E"):
                node.name = "_".join(node.name.replace("_E", "").split("_")[:-1])
            else:
                node.name = node.name.replace(" ", "_")
    return tree


def name_internal_nodes(tree: Tree) -> Tree:
    node_i = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    return tree


def convert_orthogroups_df_to_ensg(df_input: pd.DataFrame, path_ortho: str) -> pd.DataFrame:
    if path_ortho == "":
        return df_input
    df_ortho = pd.read_csv(path_ortho, sep=",", dtype=str, na_filter=False)
    # Rename the columns based on the df_ortho ENSG in human
    dico_rename = {i: (j if j.startswith("ENSG") else i) for i, j in zip(df_ortho["Orthogroups"], df_ortho["Homo_sapiens"])}
    return df_input.rename(columns=dico_rename)


def trait_heritability(trait_list: list, path_input_heritability: str) -> (list, dict):
    if path_input_heritability == "":
        return trait_list, {}
    df_h2 = pd.read_csv(path_input_heritability, sep="\t")
    set_ensg = set(trait_list)
    # Convert the df to a dictionary with the first column as key and row as value
    out_dico = {}
    for k, row in df_h2.iterrows():
        ensg = row["GENE"]
        if ensg not in set_ensg:
            continue
        h2 = row['GCTA_Sum']
        if not np.isfinite(h2):
            continue
        if h2 < 0.1:
            continue
        h2_se = row['GCTA_Sum_SE']
        out_dico[ensg] = (h2, h2_se)
    return list(out_dico.keys()), out_dico


def main(path_input_traits, path_input_orthogroups, path_input_heritability, path_input_pS, path_input_dS,
         path_output_tree, path_output_traits, path_output_var_within):
    for path in [path_input_traits, path_input_pS, path_input_dS]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [path_output_tree, path_output_traits, path_output_var_within]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    tree = open_tree(path_input_dS, format_ete3=1)
    tree = rename_tree(tree)
    set_taxa_names = set(tree.get_leaf_names())
    print(f"The tree has {len(set_taxa_names)} leaves")

    # Open the pS/hererozygosity file
    pS_df = pd.read_csv(path_input_pS, sep=",")
    # Remove all spaces in the specie name
    pS_df["species"] = pS_df["species"].apply(lambda x: x.replace(" ", "_"))
    print(f"The pS dataframe has {len(pS_df)} rows before filtering.")
    pS_df = pS_df[pS_df["species"].isin(set_taxa_names)]
    pS_col = [p for p in ["heterozygosity", "pS", 'Hzoo'] if p in pS_df.columns][0]
    print(f"Keeping only species with available '{pS_col}'.")
    pS_df = pS_df[np.isfinite(pS_df[pS_col]) & (pS_df[pS_col] > 0)]
    print(f"The pS dataframe has {len(pS_df)} rows after filtering taxon name and available pS.")
    assert len(pS_df) >= 5, "Not enough species with pS. Exiting."

    df_traits = pd.read_csv(path_input_traits)
    df_traits = convert_orthogroups_df_to_ensg(df_traits, path_input_orthogroups)
    assert "species" in df_traits.columns
    print(f"The trait dataframe has {len(df_traits)} rows before filtering taxa.")
    df_traits = df_traits[df_traits["species"].isin(set_taxa_names)]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering taxa also in the tree.")
    df_traits.to_csv(path_output_traits.replace("traits.tsv", "dataframe.tsv"), sep="\t", index=False)

    # Filter the tree and create dictionaries for variance and mean of traits
    set_taxa_names = set_taxa_names.intersection(set(df_traits["species"]))
    print(f"The intersection of the tree and trait dataframe has {len(set_taxa_names)} taxa.")
    tree = prune_tree(tree, list(set_taxa_names))
    assert len(tree.get_leaves()) == len(set_taxa_names)
    dico_var_within, dico_traits = defaultdict(list), defaultdict(list)
    for taxa_name in set_taxa_names:
        leaf_pS_df = pS_df[pS_df["species"] == taxa_name]
        if len(leaf_pS_df) == 0:
            pS = np.nan
        else:
            assert len(leaf_pS_df) == 1
            pS = float(leaf_pS_df[pS_col].iloc[0])
        dico_var_within["TaxonName"].append(taxa_name)
        dico_var_within[f"Nucleotide_diversity"].append(pS)
        dico_traits["TaxonName"].append(taxa_name)

    trait_list = df_traits.columns[2:]
    print(f"The trait dataframe has {len(trait_list)} traits before filtering heritability.")
    trait_list, heritability_dico = trait_heritability(trait_list, path_input_heritability)
    print(f"The trait dataframe has {len(trait_list)} traits after filtering heritability.")
    for trait in trait_list:
        trait_name = trait.replace(" ", "_").replace("(", "").replace(")", "")
        print(f"\nPhenotype considered is {trait}")
        trait_df = df_traits[["species", trait]].copy()
        trait_df = trait_df[np.isfinite(trait_df[trait])]
        print(f"The trait dataframe has {len(trait_df)} rows after filtering not finite values.")

        # Filter out the species with a unique row in the trait dataframe
        var_df = trait_df.copy()
        print(f"The trait dataframe has {len(var_df)} rows after filtering for Nsize != 1.")
        var_grouped = {k: v for k, v in var_df.groupby("species") if len(v) > 1}
        print(f"The trait dataframe has {len(var_grouped)} taxa after keeping taxon with more than 1 individuals.")
        for taxa_name in set_taxa_names:
            if taxa_name in var_grouped:
                phenotype_var = np.var(var_grouped[taxa_name][trait], ddof=1)
                if trait_name in heritability_dico:
                    h2_mean, h2_std = heritability_dico[trait_name]
                    h2_min, h2_max = max(0.0, h2_mean - h2_std), min(1.0, h2_mean + h2_std)
                else:
                    h2_min, h2_max = 1.0, 1.0
            else:
                phenotype_var = np.nan
                h2_min, h2_max = np.nan, np.nan
            dico_var_within[f"{trait_name}_variance"].append(phenotype_var)
            dico_var_within[f"{trait_name}_heritability_lower"].append(h2_min)
            dico_var_within[f"{trait_name}_heritability_upper"].append(h2_max)
        assert len(set_taxa_names) == len(dico_var_within[f"{trait_name}_variance"])
        print(f"{len(var_grouped)} species with variance computed.")

        mean_df = trait_df.copy()
        mean_grouped = {k: v for k, v in mean_df.groupby("species")}
        for taxa_name in set_taxa_names:
            if taxa_name in mean_grouped:
                filtered = mean_grouped[taxa_name]
                phenotype_mean = np.mean(filtered[trait])
            else:
                phenotype_mean = np.nan
            dico_traits[f"{trait_name}_mean"].append(phenotype_mean)
        print(f"{sum(np.isfinite(dico_traits[f'{trait_name}_mean']))} species with mean computed.")

    df_traits = pd.DataFrame(dico_traits)
    df_traits = df_traits[np.isfinite(df_traits.drop(["TaxonName"], axis=1)).any(axis=1)]
    df_traits.to_csv(path_output_traits, sep="\t", index=False, na_rep="NaN")
    set_taxa_names_traits = set(df_traits["TaxonName"])

    # Remove the species with only nan values across the traits
    df_var_within = pd.DataFrame(dico_var_within)
    df_var_within = df_var_within[df_var_within["TaxonName"].isin(set(pS_df["species"]))]
    df_var_within = df_var_within[df_var_within["TaxonName"].isin(set_taxa_names_traits)]
    df_var_within = df_var_within[np.isfinite(df_var_within.drop(["TaxonName", "Nucleotide_diversity"], axis=1)).any(axis=1)]
    # Write NaN for the species with no variance
    df_var_within.to_csv(path_output_var_within, sep="\t", index=False, na_rep="NaN")

    # Prune the tree and write it
    set_taxa_names = set(df_traits["TaxonName"])
    tree = prune_tree(tree, list(set_taxa_names))
    for taxa in df_var_within["TaxonName"]:
        assert taxa in set_taxa_names, f"{taxa} not in the tree"
    print(f"The final tree has {len(tree.get_leaves())} taxa.")
    tree_length = sum([node.dist for node in tree.traverse()])
    print(f"The tree length is {tree_length}.")
    tree = name_internal_nodes(tree)
    tree.write(outfile=path_output_tree, format=3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_traits", help="Input trait file", required=True)
    parser.add_argument("--input_orthogroups", help="Input orthogroups file", required=False, default="")
    parser.add_argument("--input_heritability", help="Input heritability file", required=False, default="")
    parser.add_argument("--input_pS", help="Input pS file", required=True)
    parser.add_argument("--input_dS", help="Input dS tree file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    parser.add_argument('--output_var_within', help="Output var_within file", required=True)
    args = parser.parse_args()
    main(args.input_traits, args.input_orthogroups, args.input_heritability, args.input_pS, args.input_dS,
         args.output_tree, args.output_traits, args.output_var_within)
