#!/usr/bin/env python3
import os
import argparse
import pandas as pd
from Bio import Phylo


def main(input_nuctree: str, input_timetree: str, input_traits: str,
         output_nuctree: str, output_timetree: str, output_traits: str):
    for path in [input_traits, input_nuctree, input_timetree]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_traits, output_nuctree, output_timetree]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Convert newick to nexus including root
    Phylo.convert(input_nuctree, "newick", output_nuctree, "nexus")
    Phylo.convert(input_timetree, "newick", output_timetree, "nexus")
    traits = pd.read_csv(input_traits, sep="\t")
    cols = traits.columns.tolist()
    with open(output_traits, "w") as f:
        f.write('#NEXUS\n\n')
        f.write("Begin data;\n")
        f.write(f"Dimensions ntax={len(traits)} nchar=1;\n")
        f.write("Format datatype=Continuous missing=? gap=-;\n")
        f.write("Matrix\n")
        for row in traits.itertuples():
            f.write(f"{row.TaxonName}\t{getattr(row, cols[1])}\n")
        f.write(";\n")
        f.write("End;\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_nuctree", help="Input nuc tree file", required=True)
    parser.add_argument("--input_timetree", help="Input time tree file", required=True)
    parser.add_argument("--input_traits", help="Input traits file", required=True)
    parser.add_argument("--output_nuctree", help="Output nuc tree file", required=True)
    parser.add_argument("--output_timetree", help="Output time tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    args = parser.parse_args()
    main(args.input_nuctree, args.input_timetree, args.input_traits, args.output_nuctree, args.output_timetree,
         args.output_traits)
