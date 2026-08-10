""" Gene Group Symbol Capture

Annotating gene symbol aliases with "gene group symbol" if it contains a prefix that matches the group root symbol(HGNC).

HGNC groups genes with similar functions and sequences(https://www.genenames.org/data/genegroup/#!/).

Gene Group ID and abbreviation reference from HGNC:
hgnc_genefamily.csv = family.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)
hgnc_genefamily_hierarchy = hierarchy.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)
hgnc_id_symbol_genegroupid.txt = gene_has_family.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)

This module generates: 
- hgnc_gene_group_root_df.h5 (HGNC_ID, Approved symbol, Gene group ID, abbreviation)- DataFrame mapping genes to specific and ancestral gene groups and their respective gene group root symbols
- hgnc_gene_group_analysis_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Matching Abbreviation)- subset_genes_df with alias annotations 
- hgnc_gene_group_match_subset_genes_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Matching Abbreviation)- subset_genes_df with alias annotations, only the alias symbols that are annotated as "gene group" symbols
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import re
from pybiomart import Server
from tqdm import tqdm
import json
import os
tqdm.pandas()
from collections import deque
import requests
from io import StringIO
import nest_asyncio
nest_asyncio.apply()
import ast
import time

import nltk
import gene_group_symbol_capture_functions as ggscfx
from nltk.corpus import words
from requests.exceptions import RequestException

# Import files
hgnc_genefamilies_df = pd.read_csv(
    "../input/hgnc_genefamily.csv", sep=",")

hgnc_genefamily_hierarchy_df = pd.read_csv(
    "../input/hgnc_genefamily_hierarchy.csv")

hgnc_gene_groupid_df = pd.read_csv(
    "../input/hgnc_id_symbol_genegroupid.txt", sep="\t")

subset_genes_df = pd.read_csv(
    "../output/subset_genes_df.csv")

# Convert data types of the group IDs to strings
hgnc_genefamily_hierarchy_df['child_fam_id'] = hgnc_genefamily_hierarchy_df['child_fam_id'].astype(str)
hgnc_genefamily_hierarchy_df['parent_fam_id'] = hgnc_genefamily_hierarchy_df['parent_fam_id'].astype(str)

# Extract relevant columns (group ID and group root symbol)
hgnc_genefamilies_df = hgnc_genefamilies_df[["id","abbreviation"]]
hgnc_genefamilies_df.rename(columns={'id': 'Gene group ID'}, inplace=True)

# Remove the "HGNC:" prefix from the gene ID
hgnc_gene_groupid_df['HGNC ID'] = hgnc_gene_groupid_df['HGNC ID'].str.replace('^HGNC:', '', regex=True)

# Split the gene group IDs to one per row
hgnc_gene_groupid_df['Gene group ID'] = hgnc_gene_groupid_df['Gene group ID'].str.split('|')
hgnc_gene_groupid_df = hgnc_gene_groupid_df.explode('Gene group ID')

# Expand the gene group IDs to include ancestor gene group IDs
expanded_hgnc_gene_groupid_df = ggscfx.expand_gene_groups(hgnc_gene_groupid_df, hgnc_genefamily_hierarchy_df)

# Clean up
expanded_hgnc_gene_groupid_df = expanded_hgnc_gene_groupid_df.dropna(subset=['Gene group ID'])
expanded_hgnc_gene_groupid_df['Gene group ID'] = expanded_hgnc_gene_groupid_df['Gene group ID'].astype(int)

# Add ancestor group ID info into main
hgnc_gene_group_root_df = expanded_hgnc_gene_groupid_df.merge(hgnc_genefamilies_df, on='Gene group ID', how='left')
hgnc_gene_group_root_df = hgnc_gene_group_root_df.dropna(subset=['abbreviation'])

# Remove the trailing .0 from the HGNC ID
hgnc_gene_group_root_df["HGNC ID"] = hgnc_gene_group_root_df["HGNC ID"].apply(
    lambda x: f"HGNC:{int(x)}" if pd.notna(x) and x == int(x) else f"HGNC:{x}" if pd.notna(x) else x
)

# A file with all gene in gene groups and all of their specific and ancestral gene groups and representative gene group root symbols (abbreviations)
hgnc_gene_group_root_df.to_hdf(
    "../output/hgnc_gene_group_root_df.h5", key='df', mode='w'
)

# Create uppercase helper columns for abbreviation mapping
hgnc_gene_group_root_df["Approved symbol upper"] = (
    hgnc_gene_group_root_df["Approved symbol"].str.upper()
)
hgnc_gene_group_root_df["abbreviation upper"] = (
    hgnc_gene_group_root_df["abbreviation"].str.upper()
)

# Create dictionary mapping Approved symbol -> list of abbreviations
abbrev_map = (
    hgnc_gene_group_root_df
    .dropna(subset=["Approved symbol upper", "abbreviation upper"])
    .groupby("Approved symbol upper")["abbreviation upper"]
    .apply(list)
    .to_dict()
)

# Remove temporary helper columns
hgnc_gene_group_root_df = hgnc_gene_group_root_df.drop(
    ["Approved symbol upper", "abbreviation upper"],
    axis=1
)

# Check whether an gene symbol alias symbol starts with a gene group abbreviation
hgnc_gene_group_analysis_df = subset_genes_df.copy()

# Create uppercase helper columns once
hgnc_gene_group_analysis_df["alias_symbol_upper"] = (
    hgnc_gene_group_analysis_df["alias_symbol"].str.upper()
)
hgnc_gene_group_analysis_df["primary_gene_symbol_upper"] = (
    hgnc_gene_group_analysis_df["primary_gene_symbol"].str.upper()
)

# Calculate both result columns
hgnc_gene_group_analysis_df["Gene Group Symbol Match"] = (
    hgnc_gene_group_analysis_df.progress_apply(
        ggscfx.has_prefix_gene_group_match,
        axis=1,
        args=(abbrev_map,),
    )
)

hgnc_gene_group_analysis_df["Matching Abbreviations"] = (
    hgnc_gene_group_analysis_df.progress_apply(
        ggscfx.get_matching_abbreviations,
        axis=1,
        args=(abbrev_map,),
    )
)

# Remove temporary helper columns once
hgnc_gene_group_analysis_df = hgnc_gene_group_analysis_df.drop(
    ["alias_symbol_upper", "primary_gene_symbol_upper"],
    axis=1
)

hgnc_gene_group_analysis_df.to_hdf(
    "../output/hgnc_gene_group_analysis_df.h5", key='df', mode='w'
)

# Subset of rows with "gene group symbol" annotated genes
hgnc_gene_group_match_subset_genes_df = hgnc_gene_group_analysis_df[hgnc_gene_group_analysis_df["Gene Group Symbol Match"]]

hgnc_gene_group_match_subset_genes_df.to_hdf(
    "../output/hgnc_gene_group_match_subset_genes_df.h5", key='df', mode='w'
)
