import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

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


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(tsv_traits_list: str, tsv_ML_list: str, tsv_Bayes_list: str, output: str):
    output_dir = os.path.dirname(output)
    folder_plots = f"{output_dir}/boxplots"
    os.makedirs(folder_plots, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
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
    sort_key = [i for i in ["ratio_pv", "ratio"] if i in df_out.columns][0]
    df_out = df_out.sort_values(by=[sort_key], ascending=False)
    df_out.to_csv(output, sep="\t", index=False, float_format="%.4f")
    if tsv_traits_list is None:
        return
    output_tex = f'{output_dir}/merged_boxplots.tex'
    o = open(output_tex, 'w')
    o.write(preamble)
    dico_df = {os.path.basename(i).split("_")[1]: pd.read_csv(i, sep=',') for i in tsv_traits_list}
    for dataset, group_df in df_out.groupby("dataset"):
        print(dataset)
        group_df = group_df.sort_values(by="ratio", ascending=False)
        dataset_df = dico_df[dataset]
        o.write(f"\\section{{ {dataset} }} \n")
        for _, row in group_df.iterrows():
            print(row["trait"])
            title = f"{row['trait']} in {row['dataset']} with ρ={row['ratio']:.2f}"
            if 'pp_ratio_greater_1' in row:
                title += f" (P[ρ>1]={row['pp_ratio_greater_1']:.3f})"
            # Boxplot of the trait grouped by the species
            gp = dataset_df[["species", row["trait"]]]
            # Remove the species with missing values
            print(gp)
            gp = gp.dropna()
            # Rotate the labels on the x-axis
            gp.boxplot(column=row["trait"], by="species", figsize=(12, 6))
            plt.xticks(rotation=45, ha='right')
            plt.title(title, fontsize=20)
            plt.suptitle("")
            plt.xlabel("")
            plt.tight_layout()
            output_pdf = f"{folder_plots}/{row['dataset']}_{row['trait']}.pdf"
            plt.savefig(output_pdf)
            plt.close("all")
            plt.clf()
            o.write(f"\\subsection{{ {row['trait'].replace('_', ' ')} }} \n\n")
            o.write(f"\\includegraphics[width=\\linewidth, page=1]{{ {output_pdf} }} \n\n")
        o.write("\\\\ \n")
    o.write(postamble)
    o.close()
    tex_to_pdf = f"pdflatex -synctex=1 -interaction=nonstopmode -output-directory={output_dir} {output_tex}"
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv_traits", nargs="+", help="Input tsv trait files", required=False)
    parser.add_argument("--tsv_ML", nargs="+", help="Input tsv ML files", required=False)
    parser.add_argument("--tsv_Bayes", nargs="+", help="Input tsv BayesCode files", required=False)
    parser.add_argument("--output", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.tsv_traits, args.tsv_ML, args.tsv_Bayes, args.output)
