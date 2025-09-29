import os
import argparse
from ete3 import Tree
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from neutrality_index import replace_last

preamble = """
\\documentclass{article}

\\usepackage{hyperref}
\\usepackage[margin=40pt]{geometry}
\\usepackage{amssymb,amsfonts,amsmath,amsthm,mathtools}
\\usepackage{lmodern}
\\usepackage{bm,bbold}
\\usepackage{verbatim}
\\usepackage{float}
\\usepackage[export]{adjustbox}
\\usepackage{graphicx}
\\pdfinclusioncopyfonts=1
\\begin{document}

\\tableofcontents \n
"""

postamble = """
\\end{document}
"""


def open_orthogroups(input_orthogroups: str) -> dict:
    df_ortho = pd.read_csv(input_orthogroups, sep=",", dtype=str, na_filter=False)
    # Rename the columns based on the df_ortho ENSG in human
    dico_rename = {(j if j.startswith("ENSG") else i): i for i, j in
                   zip(df_ortho["Orthogroups"], df_ortho["Homo_sapiens"])}
    return dico_rename


def find_between(input_str: str, begin: str, last: str):
    assert begin in input_str
    assert last in input_str
    i = input_str.index(begin)
    j = input_str.index(last)
    return input_str[i + len(begin): j]


def main(tsv_traits_list: str, tree_file: str, tsv_Bayes: str, input_orthogroups: str, output_pdf: str):
    output_dir = os.path.dirname(output_pdf)
    t = Tree(tree_file)
    # Traverse the tree and get species names
    species_names = [n.name for n in t.traverse("postorder") if n.is_leaf()]
    species_order = {v: k for k, v in enumerate(species_names)}
    folder_plots = f"{output_dir}/boxplots"
    os.makedirs(folder_plots, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    df_out = pd.read_csv(tsv_Bayes, sep='\t')
    if "pp_ratio_greater_1" in df_out.columns:
        df_out = df_out[df_out["pp_ratio_greater_1"] > 0.95]
    output_tex = replace_last(output_pdf, '.pdf', '.tex')
    o = open(output_tex, 'w')
    o.write(preamble)
    dico_ortho = open_orthogroups(input_orthogroups)
    dico_df = {}
    for i in tsv_traits_list:
        dico_df[os.path.basename(i).replace(".csv", "")] = pd.read_csv(i, sep=',')

    for dataset, group_df in df_out.groupby("dataset"):
        print(dataset)
        group_df = group_df.sort_values(by="ratio", ascending=False)
        dataset_df = dico_df[dataset]
        o.write(f"\\section{{ {dataset.replace('_', ' ')} }} \n")
        for _, row in group_df.iterrows():
            print(row["trait"])
            title = f"{row['trait']} in {row['dataset']} with ρ={row['ratio']:.2f}"
            if 'pp_ratio_greater_1' in row:
                title += f" (P[ρ>1]={row['pp_ratio_greater_1']:.3f})"
            # Boxplot of the trait grouped by the species
            orthogroup = dico_ortho[row["trait"]]
            gp = dataset_df[["species", orthogroup]]
            # Remove the species with missing values
            gp = gp.dropna()
            # filter the species to keep only those in the tree
            gp = gp[gp["species"].isin(set(species_names))]
            # Sort the species by the order in the tree
            gp["species_index"] = gp["species"].map(lambda x: species_order[x])
            gp = gp.sort_values("species_index")
            # Remove "_" in the species names and capitalize the first letter of each word
            gp["species"] = gp["species"].str.replace("_", " ").str.title()
            # Rotate the labels on the x-axis
            sns.boxplot(data=gp, x="species", y=orthogroup)
            plt.xticks(rotation=45, ha='right')
            plt.title(title, fontsize=20)
            plt.suptitle("")
            plt.xlabel("")
            plt.tight_layout()
            output_box = f"{folder_plots}/{row['dataset']}_{row['trait']}.pdf"
            plt.savefig(output_box)
            plt.close("all")
            plt.clf()
            o.write(f"\\subsection{{ {row['trait'].replace('_', ' ')} }} \n")
            o.write(f"\\includegraphics[width=\\linewidth, page=1]{{ {output_box} }} \n\n")
        o.write("\\\\ \n")
    o.write(postamble)
    o.close()
    tex_to_pdf = f"pdflatex -synctex=1 -interaction=nonstopmode -output-directory={output_dir} {output_tex}"
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv_traits", nargs="+", help="Input tsv trait files", required=False)
    parser.add_argument("--tree_file", help="Input tree file", required=False)
    parser.add_argument("--tsv_Bayes", help="Input tsv BayesCode", required=False)
    parser.add_argument("--input_orthogroups", help="Input orthogroups file", required=False, default="")
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.tsv_traits, args.tree_file, args.tsv_Bayes, args.input_orthogroups, args.output)
