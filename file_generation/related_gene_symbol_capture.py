""" Related Gene Symbol Capture

Anotating gene symbol aliases with "related gene symbol" if it  is identical to the related gene(NCBI Gene).

NCBI Gene groups genes with the following relationships to other genes: (https://ftp.ncbi.nlm.nih.gov/gene/DATA/README)
               - Potential readthrough sibling
               - Readthrough child
               - Readthrough parent
               - Readthrough sibling
               - Region member
               - Region parent
               - Related functional gene
               - Related pseudogene

Related gene reference from NCBI Gene:
gene_group20251007 = gene_group.gz(from https://ftp.ncbi.nlm.nih.gov/gene/DATA/)

This module generates: 
- ncbi_related_gene_analysis_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Relationship)- subset_genes_df with alias annotations 
- ncbi_related_gene_match_subset_genes_df (primary_gene_symbol, alias_symbol, HGNC_ID, ENSG_ID, NCBI_ID, Gene Group Symbol Match, Relationship)- subset_genes_df with alias annotations, only the alias symbols that are annotated as "related gene" symbols

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
from nltk.corpus import words
from requests.exceptions import RequestException

# Import
raw_ncbi_related_gene_df = pd.read_csv(
    "../input/gene_group20251007",
    sep="\t",
    usecols=["GeneID", "relationship","Other_GeneID"]
)

subset_genes_df = pd.read_csv(
    "../output/subset_genes_df.csv")

# Clean up
ncbi_related_gene_df = raw_ncbi_related_gene_df.replace('-', np.nan)
ncbi_related_gene_df = ncbi_related_gene_df.dropna(subset=['Other_GeneID'])
ncbi_related_gene_df = ncbi_related_gene_df.rename(
    columns={
        "GeneID": "NCBI_ID",
        "relationship": "Relationship"
    }
)

ncbi_related_gene_analysis_df = subset_genes_df.copy()

ncbi_related_gene_analysis_df["NCBI_ID"] = (
    ncbi_related_gene_analysis_df["NCBI_ID"]
    .apply(
        lambda x: list(ast.literal_eval(x))
        if isinstance(x, str)
        else []
    )
)

ncbi_related_gene_analysis_df = (
    ncbi_related_gene_analysis_df
    .explode("NCBI_ID")
    .reset_index(drop=True)
)

ncbi_related_gene_analysis_df = (
    ncbi_related_gene_analysis_df
    .dropna(subset=["alias_symbol", "NCBI_ID"])
)

ncbi_related_gene_analysis_df["NCBI_ID"] = (
    ncbi_related_gene_analysis_df["NCBI_ID"]
    .str.removeprefix("GENE ID:")
)

cols_to_clean = ['NCBI_ID', 'Other_GeneID']
ncbi_related_gene_df[cols_to_clean] = ncbi_related_gene_df[cols_to_clean].astype(str).apply(lambda col: col.str.strip())
ncbi_related_gene_analysis_df['NCBI_ID'] = ncbi_related_gene_analysis_df['NCBI_ID'].astype(str).str.strip()

# Add gene symbol aliases based on NCBI ID
ncbi_related_gene_with_aliases_df = ncbi_related_gene_df.merge(
    ncbi_related_gene_analysis_df[['NCBI_ID', 'alias_symbol']],
    left_on='NCBI_ID',
    right_on='NCBI_ID',
    how='left'
)
ncbi_related_gene_with_aliases_df = ncbi_related_gene_with_aliases_df.dropna(subset=['alias_symbol'])

# Add the primary gene symbol of the related gene
ncbi_related_gene_with_aliases_and_symbols_df = ncbi_related_gene_with_aliases_df.merge(
    ncbi_related_gene_analysis_df[['NCBI_ID', 'primary_gene_symbol']],
    left_on='Other_GeneID',
    right_on='NCBI_ID',
    how='left',
    suffixes=('', '_associated')
)

# Check whether an gene symbol alias symbol is a primary gene symbol of a related gene
ncbi_related_gene_with_aliases_and_symbols_df['Related Gene Symbol Match'] = (
    ncbi_related_gene_with_aliases_and_symbols_df['alias_symbol'] == ncbi_related_gene_with_aliases_and_symbols_df['primary_gene_symbol']
)

ncbi_related_gene_analysis_df = ncbi_related_gene_analysis_df.merge(
    ncbi_related_gene_with_aliases_and_symbols_df[['NCBI_ID', 'alias_symbol', 'Related Gene Symbol Match', 'Relationship']],
    on=['NCBI_ID', 'alias_symbol'],
    how='left'
)

ncbi_related_gene_analysis_df['NCBI_ID'] = ncbi_related_gene_analysis_df['NCBI_ID'].apply(
        lambda x: set(x.split(',')) if isinstance(x, str) and x else set() if x == '' else x
    )

ncbi_related_gene_analysis_df = ncbi_related_gene_analysis_df[ncbi_related_gene_analysis_df["Related Gene Symbol Match"].fillna(False)]
ncbi_related_gene_analysis_df.loc[ncbi_related_gene_analysis_df['Related Gene Symbol Match'] == False, 'Relationship'] = np.nan
ncbi_related_gene_match_subset_genes_df = ncbi_related_gene_analysis_df.loc[ncbi_related_gene_analysis_df["Related Gene Symbol Match"]]

ncbi_related_gene_analysis_df.to_hdf(
    "../output/ncbi_related_gene_analysis_df.h5", key='df', mode='w'
)

ncbi_related_gene_match_subset_genes_df.to_hdf(
    "../output/ncbi_related_gene_match_subset_genes_df.h5", key='df', mode='w'
)
