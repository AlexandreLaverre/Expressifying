#!/usr/bin/env python3
# coding=utf-8

import glob
import pandas as pd
import networkx as nx
from collections import defaultdict
import csv

path = "/Users/alaverre/Documents/Expressifying/"
species_file = path + "data/Bgee_libraries/species_list.csv"
all_orthogroups = path + "data/gene_orthologies/all_orthogroups_2022_mammals.csv"
one2one_orthogroups = path + "data/gene_orthologies/one2one_orthogroups_2022_mammals.csv"

########################################################################################################################
# Read matching between NCBI ID and species names
species_index = {}
with open(species_file, "r") as match_file:
    reader = csv.reader(match_file)
    next(reader)  # Skip the header row
    for row in reader:
        if row[11] in ["Mammal", "Marsupial"]:
            NCBI_ID = row[0]
            #sp_name = row[1].replace(" ", "_")
            sp_name = row[9].replace(" ", "_")
            species_index[NCBI_ID] = sp_name

all_species = species_index.keys()
print(f"Number of species: {len(all_species)}")

########################################################################################################################
print("Finding all pairwise ortholog genes...")

file_pattern = path + "data/gene_orthologies/Bgee-OMA-orthologs-Nov2022/orthologs_*.csv"
files = glob.glob(file_pattern)

nb_pairwise = 1
all_pairs = []
for file in files:
    # Extract the species names from the file name
    file_name = file.split("/")[-1]  # Get the file name without the directory path
    species_names = file_name.split("_")[-1].split(".")[0].split("-")
    ref_sp = species_names[0]
    tg_sp = species_names[1]

    if ref_sp in all_species and tg_sp in all_species:
        df = pd.read_csv(file)

        # Edit gene names to integrate species ID
        df['gene1'] = ref_sp + '-' + df['gene1'].astype(str)
        df['gene2'] = tg_sp + '-' + df['gene2'].astype(str)

        # Combine columns to create the list of tuples IDs
        ids = list(zip(df['gene1'], df['gene2']))
        all_pairs.extend(ids)
        nb_pairwise += 1

print(f"Number of pairwise files analyzed: {nb_pairwise}")
print(f"Total number of gene pairs: {len(all_pairs)}")

########################################################################################################################
# Make connected graph
print("Making connected graph of orthologous pairs...")
G = nx.Graph()
G.add_edges_from(all_pairs)

# Attribute genes to species in each orthogroup
species_genes = {}
nb_group = 1
for orthogroup in nx.connected_components(G):
    orthogroup_dic = defaultdict(list)
    orthogroup_ID = "Orthogroup_" + str(nb_group)

    for gene in orthogroup:
        sp_ID, gene_ID = gene.split('-')
        sp_name = species_index[sp_ID]  # retrieve complete species name
        orthogroup_dic[sp_name].append(gene_ID)

    species_genes[orthogroup_ID] = orthogroup_dic
    nb_group += 1

print(f"Number of Orthogroups: {nb_group}")

########################################################################################################################
print("Writing output...")
# Create the output file and write the header
with open(all_orthogroups, "w", newline="") as file_all, open(one2one_orthogroups, "w", newline="") as file_one2one:
    writer_all = csv.writer(file_all)
    writer_one2one = csv.writer(file_one2one)

    # Write the header row
    header = ['Orthogroups'] + list(species_index.values())
    writer_all.writerow(header)
    writer_one2one.writerow(header)

    # Write the gene information for each orthogroup
    for orthogroup, genes in species_genes.items():
        # Initialize the row with NA for all species
        row_all = [orthogroup] + ['NA'] * len(all_species)
        row_one2one = [orthogroup] + ['NA'] * len(all_species)

        # Update the row with the actual gene names
        for species, species_genes in genes.items():
            col_index = header.index(species)
            row_all[col_index] = ','.join(species_genes)
            if len(species_genes) == 1:
                row_one2one[col_index] = species_genes[0]

        writer_all.writerow(row_all)
        writer_one2one.writerow(row_one2one)

########################################################################################################################
