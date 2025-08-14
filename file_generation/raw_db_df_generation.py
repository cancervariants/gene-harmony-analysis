import pandas as pd
import numpy as np
import re
from datetime import datetime

#Ensembl
raw_ensg_df = pd.read_csv(
    "input/ensg_biomart_gene20250625.txt", sep="\t",dtype={"NCBI gene (formerly Entrezgene) ID": pd.Int64Dtype()}
)
raw_ensg_df = raw_ensg_df.rename(
    columns={
        "HGNC ID": "HGNC_ID",
        "Gene Synonym": "alias_symbol",
        "Gene name": "gene_symbol",
        "Gene stable ID": "ENSG_ID",
        "NCBI gene (formerly Entrezgene) ID": "NCBI_ID",
    }
)

#Format the ID values as strings. Cast to int before formatting to remove the .0
raw_ensg_df["NCBI_ID"] = raw_ensg_df["NCBI_ID"].apply(
lambda x: f"GENE ID:{int(x)}" if pd.notna(x) and x == int(x) else f"GENE ID:{x}" if pd.notna(x) else x
) 

#Replace empty strings with nans
raw_ensg_df = raw_ensg_df.replace(" ", np.nan)
raw_ensg_df = raw_ensg_df.replace("", np.nan)
raw_ensg_df = raw_ensg_df.replace("-", np.nan)
raw_ensg_df = raw_ensg_df.replace("<NA>", np.nan)


raw_ensg_df.to_csv('output/raw_ensg_df.csv', index=False) 

#HGNC
hgnc_file_path = "input/hgnc_biomart_gene20250625.txt"

raw_hgnc_df = pd.read_csv(
    hgnc_file_path, sep="\t"
)

# Rename columns
raw_hgnc_df = raw_hgnc_df.rename(
    columns={
        "HGNC ID": "HGNC_ID",
        "Approved symbol": "gene_symbol",
        "Ensembl gene ID": "ENSG_ID",
    }
)

#structure and labeling in HGNC download files changed between 2024 amd 2025
if "Alias symbol" in raw_hgnc_df.columns:
    raw_hgnc_df = raw_hgnc_df.rename(columns={"Alias symbol": "alias_symbol"})
elif "Alias symbols" in raw_hgnc_df.columns:
    raw_hgnc_df = raw_hgnc_df.rename(columns={"Alias symbols": "alias_symbol"})
else:
    raw_hgnc_df["alias_symbol"] = pd.NA

if "NCBI gene ID" in raw_hgnc_df.columns:
    raw_hgnc_df = raw_hgnc_df.rename(columns={"NCBI gene ID": "NCBI_ID"})
elif "NCBI Gene ID" in raw_hgnc_df.columns:
    raw_hgnc_df = raw_hgnc_df.rename(columns={"NCBI Gene ID": "NCBI_ID"})
else:
    raw_hgnc_df["NCBI_ID"] = pd.NA   

raw_hgnc_df["NCBI_ID"] = raw_hgnc_df["NCBI_ID"].astype(pd.Int64Dtype())

# Extract date from filename and check if it is before June 25, 2025
match = re.search(r'(\d{8})', hgnc_file_path)
file_date = datetime.strptime(match.group(1), "%Y%m%d") if match else None
cutoff_date = datetime.strptime("20250625", "%Y%m%d")

# Apply list-splitting logic for newer files
if file_date and file_date >= cutoff_date:
    raw_hgnc_df['alias_symbol'] = (
        raw_hgnc_df['alias_symbol']
        .fillna('')
        .str.split(',')
        .apply(lambda x: [a.strip() for a in x if a.strip()])
    )
    raw_hgnc_df = raw_hgnc_df.explode('alias_symbol')

#Format the ID values as strings. Cast to int before formatting to remove the .0
raw_hgnc_df["NCBI_ID"] = raw_hgnc_df["NCBI_ID"].apply(
lambda x: f"GENE ID:{int(x)}" if pd.notna(x) and x == int(x) else f"GENE ID:{x}" if pd.notna(x) else x
) 

#Replace empty strings with nans
raw_hgnc_df = raw_hgnc_df.replace(" ", np.nan)
raw_hgnc_df = raw_hgnc_df.replace("", np.nan)
raw_hgnc_df = raw_hgnc_df.replace("-", np.nan)

#Export as csv
raw_hgnc_df.to_csv('output/raw_hgnc_df.csv', index=False) 

#NCBI Gene
ncbi_file_path = "input/Homo_sapiens.gene_info20250625"
raw_ncbi_df = pd.read_csv(ncbi_file_path, sep="\t")

# Drop all columns besides ENSG_ID, gene_symbol, alias_symbol, and type_of_gene
raw_ncbi_df = raw_ncbi_df[
["GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene"]
]
raw_ncbi_df = raw_ncbi_df.rename(
    columns={"GeneID": "NCBI_ID", "Symbol": "gene_symbol", "Synonyms": "alias_symbol", "type_of_gene": "gene_type"})

#SplitdbXrefs into individual columns
for col in ["MIM", "HGNC_ID", "ENSG_ID", "AllianceGenome", "MIRbase", "IMGTgene_db", "dash", "unknown"]:
    raw_ncbi_df[col] = pd.Series(dtype="object")

index_pos = 0

while index_pos < len(raw_ncbi_df):
    xrefs = raw_ncbi_df["dbXrefs"][index_pos].split("|")

    for xref in xrefs:
        xref = xref.lower()
        if xref.startswith("mim:"):
            xref = xref.replace("mim:", "")
            raw_ncbi_df.loc[index_pos, "MIM"] = xref
        elif xref.startswith("hgnc:hgnc:"):
            xref = xref.replace("hgnc:hgnc:", "")
            raw_ncbi_df.loc[index_pos, "HGNC_ID"] = xref
        elif xref.startswith("ensembl:"):
            xref = xref.replace("ensembl:", "")
            raw_ncbi_df.loc[index_pos, "ENSG_ID"] = xref
        elif xref.startswith("alliancegenome:"):
            xref = xref.replace("alliancegenome:", "")
            raw_ncbi_df.loc[index_pos, "AllianceGenome"] = xref
        elif xref.startswith("mirbase"):
            xref = xref.replace("mirbase:", "")
            raw_ncbi_df.loc[index_pos, "MIRbase"] = xref
        elif xref.startswith("imgt/gene-db:"):
            xref = xref.replace("imgt/gene-db:", "")
            raw_ncbi_df.loc[index_pos, "IMGTgene_db"] = xref
        elif xref.startswith("-"):
            raw_ncbi_df.loc[index_pos, "dash"] = xref
        else:
            raw_ncbi_df.loc[index_pos, "unknown"] = xref

    index_pos += 1
    pass

#Alter gene ID to be consistent
raw_ncbi_df["ENSG_ID"] = raw_ncbi_df["ENSG_ID"].str.replace("ensg", "ENSG", 1)

#Drop unused columns
raw_ncbi_df = raw_ncbi_df.drop(
    [
        "AllianceGenome",
        "MIRbase",
        "IMGTgene_db",
        "dash",
        "unknown",
        "dbXrefs",
        "MIM",
    ],
    axis=1,
)

# Extract date from filename and check if it is before June 25, 2025
match = re.search(r'(\d{8})', ncbi_file_path)
file_date = datetime.strptime(match.group(1), "%Y%m%d") if match else None
cutoff_date = datetime.strptime("20250625", "%Y%m%d")

# Apply list-splitting logic for newer files
if file_date and file_date >= cutoff_date:
    raw_ncbi_df['alias_symbol'] = (
        raw_ncbi_df['alias_symbol']
        .fillna('')
        .str.split('|')
        .apply(lambda x: [a.strip() for a in x if a.strip()])
    )
    raw_ncbi_df = raw_ncbi_df.explode('alias_symbol')

#Format the ID values as strings. Cast to int before formatting to remove the .0
raw_ncbi_df["NCBI_ID"] = raw_ncbi_df["NCBI_ID"].apply(
lambda x: f"GENE ID:{int(x)}" if pd.notna(x) and x == int(x) else f"GENE ID:{x}" if pd.notna(x) else x
)

raw_ncbi_df["HGNC_ID"] = raw_ncbi_df["HGNC_ID"].apply(
lambda x: f"HGNC:{int(x)}" if pd.notna(x) and x == int(x) else f"HGNC:{x}" if pd.notna(x) else x
)

#Replace empty strings with nans
raw_ncbi_df = raw_ncbi_df.replace(" ", np.nan)
raw_ncbi_df = raw_ncbi_df.replace("", np.nan)
raw_ncbi_df = raw_ncbi_df.replace("-", np.nan)

#Export as csv
raw_ncbi_df.to_csv('output/raw_ncbi_df.csv', index=False) 