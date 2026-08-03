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

# Download gene records from ENSG, HGNC, and NCBI
# These subset_genes_xxxx_dfs was created by the following modifications from the mini_xxxx_dfs
# - Primary gene symbol- alias symbol pairs where the alias was an exact match to the primary symbol were removed.
# - Primary gene symbol- alias symbol pairs that were duplicated were removed.

databases = ['ensg','hgnc','ncbi']

dfs = [
    scgfx.read_subset_genes_csv(f'../output/subset_genes_{db}_df.csv')
    for db in databases
]

subset_genes_df = pd.concat(dfs, axis=0, ignore_index=True)

subset_genes_df = subset_genes_df.groupby(
    ["primary_gene_symbol", "alias_symbol"],
    as_index=False
).agg({
    "HGNC_ID": lambda x: set(x.dropna()),
    "ENSG_ID": lambda x: set(x.dropna()),
    "NCBI_ID": lambda x: set(x.dropna()),
})

subset_genes_df['NCBI_ID'] = subset_genes_df['NCBI_ID'].apply(scgfx.remove_nan_from_set)
subset_genes_df['ENSG_ID'] = subset_genes_df['ENSG_ID'].apply(scgfx.remove_nan_from_set)
subset_genes_df['HGNC_ID'] = subset_genes_df['HGNC_ID'].apply(scgfx.remove_nan_from_set)

subset_genes_df.to_hdf("../output/subset_genes_df.h5", key='df', mode='w')
subset_genes_df.to_csv("../output/subset_genes_df.csv", index=True)