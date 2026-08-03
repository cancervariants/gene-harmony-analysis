import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pybiomart import Server
from tqdm import tqdm

tqdm.pandas()
from collections import deque
import requests
from io import StringIO
import nest_asyncio
nest_asyncio.apply()
import ast
import time

import nltk
import symbol_capture_generation_functions as scgfx
from automatic_alternate_abbreviation_annotation.alt_abbrev_llm_functions import *
from nltk.corpus import words
from requests.exceptions import RequestException

subset_genes_df = pd.read_csv("../output/subset_genes_df.csv")

## <a id='toc2_1_'></a>[Download an Ensembl Biomart export file with the Gene Name and the Ortholog Gene Name](#toc0_)
server = Server(host='http://www.ensembl.org')
dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

homologs = []
for key in dataset.attributes.keys():
    if '_homolog_associated_gene_name' in key:
        homologs.append(key)
print(f'{len(homologs)} found!')

for homolog in tqdm(homologs):
    attributes = ['ensembl_gene_id', 'external_gene_name', homolog]
    df = dataset.query(attributes=attributes)
    df.to_csv(f'../input/orthologs/{homolog}.csv')
## <a id='toc2_2_'></a>[Match aliases to orthologs!](#toc0_)
#can comment out after run
BATCH_SIZE = 5

files = sorted([file for file in os.listdir('../input/orthologs/') if file.endswith('.csv')])

for i in tqdm(range(0, len(files), BATCH_SIZE)):
    batch_files = files[i:i + BATCH_SIZE]
    
    batch_df = pd.concat(
        [pd.read_csv(f'../input/orthologs/{file}') for file in batch_files],
        ignore_index=True
    )
    
    analysis, counts = scgfx.match_alias_to_ortholog(subset_genes_df, batch_df)
    
    batch_number = i // BATCH_SIZE
    analysis.to_csv(f'../output/orthologs/analysis_{batch_number}.csv')

    counts_serialized = {k: int(v) for k, v in counts.items()} # Required to serialize numpy to int, otherwise cannot save to json
    with open(f'../output/orthologs/counts_{batch_number}.json', 'w') as f:
        json.dump(counts_serialized, f, indent=2)

## <a id='toc2_3_'></a>[SGD Ortholog Symbol](#toc0_)
Not all of the yeast orthologs are present in the Ensembl database. The Saccharomyces Genome Database was identified to potentially fill in the gap through manual annotation of the alias symbols LYS5 (correlates to HGNC:14235) and BFR2 (HGNC:19235)
In the Download tab on the SGD site: http://sgd-archive.yeastgenome.org / curation / chromosomal_feature / dbxref.tab
sgd_dbref_df = pd.read_csv(
    "../input/dbxref20251213.tab", 
    sep="\t", 
    header=None)
#from the dbxref.README
sgd_dbref_df.columns = [
    "DBXREF ID", 
    "DBXREF ID source", 
    "DBXREF ID type", 
    "S. cerevisiae feature name", 
    "SGDID",
    "gene_symbol"
]
isolate the rows that point the yeast gene to a HGNC xref
hgnc_sgd_dbref_df = sgd_dbref_df[sgd_dbref_df["DBXREF ID"].fillna('').str.startswith("HGNC:")]