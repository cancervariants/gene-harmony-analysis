""" Gene Group Symbol Capture

Anotating gene symbol aliases with "gene group symbol" if it contains a prefix or completely matches the group root symbol(HGNC) or related gene(NCBI Gene).

HGNC groups genes with similar functions and sequences(https://www.genenames.org/data/genegroup/#!/).
NCBI Gene groups genes with the following relationships to other genes: (https://ftp.ncbi.nlm.nih.gov/gene/DATA/README)
               - Potential readthrough sibling
               - Readthrough child
               - Readthrough parent
               - Readthrough sibling
               - Region member
               - Region parent
               - Related functional gene
               - Related pseudogene
Gene Group ID and abbreviation reference from HGNC:
hgnc_genefamily.csv = family.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)
hgnc_genefamily_hierarchy = hierarchy.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)
hgnc_id_symbol_genegroupid.txt = gene_has_family.csv (from https://www.genenames.org/download/gene-groups/#!/#tocAnchor-1-1)
Gene Group reference from NCBI Gene:
gene_group20251007 = gene_group.gz(from https://ftp.ncbi.nlm.nih.gov/gene/DATA/)

This module generates: 
- hgnc_gene_group_root_df.h5 (HGNC_ID, Approved symbol, Gene group ID, abbreviation)- DataFrame mapping genes to specific and ancestral gene groups and their respective gene group root symbols
- hgnc_gene_group_analysis_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Matching Abbreviation)- subset_genes_df with alias annotations 
- hgnc_gene_group_match_subset_genes_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Matching Abbreviation)- subset_genes_df with alias annotations, only the alias symbols that are annotated as "gene group" symbols
- ncbi_gene_group_analysis_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, relationship)- subset_genes_df with alias annotations 
- ncbi_gene_group_match_subset_genes_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, relationship)- subset_genes_df with alias annotations, only the alias symbols that are annotated as "gene group" symbols

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
import gene_group_symbol_capture_functions as pggscfx
from nltk.corpus import words
from requests.exceptions import RequestException

# HGNC

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
expanded_hgnc_gene_groupid_df = pggscfx.expand_gene_groups(hgnc_gene_groupid_df, hgnc_genefamily_hierarchy_df)

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
        pggscfx.has_prefix_gene_group_match,
        axis=1
    )
)

hgnc_gene_group_analysis_df["Matching Abbreviations"] = (
    hgnc_gene_group_analysis_df.progress_apply(
        pggscfx.get_matching_abbreviations,
        axis=1
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

# NCBI Gene

# Import
raw_ncbi_gene_group_df = pd.read_csv(
    "../input/gene_group20251007",
    sep="\t",
    usecols=["GeneID", "relationship","Other_GeneID"]
)

# Clean up
ncbi_gene_group_df = raw_ncbi_gene_group_df.replace('-', np.nan)
ncbi_gene_group_df = ncbi_gene_group_df.dropna(subset=['Other_GeneID'])
ncbi_gene_group_df = ncbi_gene_group_df.rename(columns={'GeneID': 'NCBI_ID'})

import ast

ncbi_gene_group_analysis_df = subset_genes_df.copy()

ncbi_gene_group_analysis_df["NCBI_ID"] = (
    ncbi_gene_group_analysis_df["NCBI_ID"]
    .apply(
        lambda x: list(ast.literal_eval(x))
        if isinstance(x, str)
        else []
    )
)

ncbi_gene_group_analysis_df = (
    ncbi_gene_group_analysis_df
    .explode("NCBI_ID")
    .reset_index(drop=True)
)

ncbi_gene_group_analysis_df = (
    ncbi_gene_group_analysis_df
    .dropna(subset=["alias_symbol", "NCBI_ID"])
)

ncbi_gene_group_analysis_df["NCBI_ID"] = (
    ncbi_gene_group_analysis_df["NCBI_ID"]
    .str.removeprefix("GENE ID:")
)

cols_to_clean = ['NCBI_ID', 'Other_GeneID']
ncbi_gene_group_df[cols_to_clean] = ncbi_gene_group_df[cols_to_clean].astype(str).apply(lambda col: col.str.strip())
ncbi_gene_group_analysis_df['NCBI_ID'] = ncbi_gene_group_analysis_df['NCBI_ID'].astype(str).str.strip()

# Add gene symbol aliases based on NCBI ID
ncbi_gene_group_with_aliases_df = ncbi_gene_group_df.merge(
    ncbi_gene_group_analysis_df[['NCBI_ID', 'alias_symbol']],
    left_on='NCBI_ID',
    right_on='NCBI_ID',
    how='left'
)
ncbi_gene_group_with_aliases_df = ncbi_gene_group_with_aliases_df.dropna(subset=['alias_symbol'])

# Add the primary gene symbol of the related gene
ncbi_gene_group_with_aliases_and_symbols_df = ncbi_gene_group_with_aliases_df.merge(
    ncbi_gene_group_analysis_df[['NCBI_ID', 'primary_gene_symbol']],
    left_on='Other_GeneID',
    right_on='NCBI_ID',
    how='left',
    suffixes=('', '_associated')
)

# Check whether an gene symbol alias symbol is a primary gene symbol of a related gene
ncbi_gene_group_with_aliases_and_symbols_df['Gene Group Symbol Match'] = (
    ncbi_gene_group_with_aliases_and_symbols_df['alias_symbol'] == ncbi_gene_group_with_aliases_and_symbols_df['primary_gene_symbol']
)

ncbi_gene_group_analysis_df = ncbi_gene_group_analysis_df.merge(
    ncbi_gene_group_with_aliases_and_symbols_df[['NCBI_ID', 'alias_symbol', 'Gene Group Gene Symbol Match', 'relationship']],
    on=['NCBI_ID', 'alias_symbol'],
    how='left'
)

ncbi_gene_group_analysis_df['NCBI_ID'] = ncbi_gene_group_analysis_df['NCBI_ID'].apply(
        lambda x: set(x.split(',')) if isinstance(x, str) and x else set() if x == '' else x
    )

ncbi_gene_group_analysis_df = ncbi_gene_group_analysis_df[ncbi_gene_group_analysis_df["Gene Group Gene Symbol Match"].fillna(False)]
ncbi_gene_group_analysis_df.loc[ncbi_gene_group_analysis_df['Gene Group Gene Symbol Match'] == False, 'relationship'] = np.nan
ncbi_gene_group_match_subset_genes_df = ncbi_gene_group_analysis_df.loc[ncbi_gene_group_analysis_df["Gene Group Gene Symbol Match"]]

ncbi_gene_group_analysis_df.to_hdf(
    "../output/ncbi_gene_group_analysis_df.h5", key='df', mode='w'
)

ncbi_gene_group_match_subset_genes_df.to_hdf(
    "../output/ncbi_gene_group_match_subset_genes_df.h5", key='df', mode='w'
)

# Merge the HGNC and NCBI gene_group_analysis_dfs
