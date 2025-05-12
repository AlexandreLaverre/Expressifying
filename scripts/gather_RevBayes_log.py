import os
from os.path import isdir
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd


def format_label(l):
    l = l.replace("b ", "").replace("m ", "♂").replace("f ", "♀ ")
    l = l.replace("bodyMass", "body mass").replace("brainMass", "brain mass")
    return l.replace("REML", "(REML)").replace("nodes", "")


def get_number_taxa(nexus_file):
    with open(nexus_file) as f:
        for line in f:
            if "ntax" in line:
                i = line.find("ntax=")
                return int(line[i + 5:].split(" ")[0])
    return -1


def main(folder, output_tsv):
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
    models = {"REML": "simple_BM_SwitchREML", "nodes": "simple_BM_Switchnodes"}
    x_input = {}
    output_dico = defaultdict(list)
    if not isdir(folder):
        df = pd.DataFrame(output_dico)
        df.to_csv(output_tsv, sep="\t", index=False)
        return
    for directory in sorted(os.listdir(folder)):
        if not isdir(f"{folder}/{directory}"):
            continue
        for rb_name, rb in models.items():
            path = f"{folder}/{directory}/{rb}.log.gz"
            if not os.path.exists(path):
                continue
            label = f"{directory}_{rb_name}"
            df_rb = pd.read_csv(path, sep="\t")

            for col in df_rb.columns:
                if ("[" in col) or ("]" in col) or ("traitroot" == col):
                    continue
                output_dico[col].append(np.nanmean(df_rb[col]))
            output_dico["rb"].append(rb_name)
            output_dico["simu"].append(directory)
            output_dico["ntax"].append(get_number_taxa(f"{folder}/{directory}/traits.nex"))
            # cut in 10 slices
            p_slice = np.array_split(df_rb["is_nuc"], 10)
            print(len(df_rb), "points in", label)
            # compute the mean across 100 points
            x_input[label] = [np.mean(p) for p in p_slice]

    df = pd.DataFrame(output_dico)
    df.to_csv(output_tsv, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-t', '--output_tsv', required=True, type=str, dest="output_tsv", help="Output tsv path")
    args = parser.parse_args()
    main(args.folder, args.output_tsv)
