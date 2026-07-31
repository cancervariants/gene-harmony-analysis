###Automatic Alternate Abbreviation Annotation:** Adding and processing gene-alias pairs
# 
# Running this module is optional. It adds samples to the manually annotated Alternate Abbreviations dataset using a partially random selection process that preserves the existing proportion of captured versus non-captured alias symbols. For more information on Alternate Abbreviations and workflow please see [README](alternate_abbreviation_README).
###

import ast
from pathlib import Path

import polars as pl
from tqdm.notebook import tqdm
from alt_abbrev_llm_functions import *

# %%
ALT_ABBREV_ROOT = Path.cwd().resolve()
ALT_ABBREV_OUTPUT_PATH = ALT_ABBREV_ROOT / "output"

# %% [markdown]
# ## Define the number of new samples to add to the existing dataset

# %%
NUMBER_OF_NEW_SAMPLES_TO_ADD = 2

# %% [markdown]
# ## Load dataset with gene-alias pairs manually curated for Alternate Abbreviation alias symbols

# %%
curated_df = pl.read_excel(
    ALT_ABBREV_OUTPUT_PATH / "alt_abbrev_annotation_manually_annotated_df.xlsx"
)

# %%
dedup_keys = ["primary_gene_symbol", "gene_symbol"]
curated_keys = curated_df.select(dedup_keys).unique()

# %% [markdown]
# ## Load and clean dataset where to get new gene-alias pairs from

# %%
## This is the file with alias and primary gene symbol pairs and what categories the aliases were captured as
## Generated in the 5_symbol_capture_analysis.ipynb
## Clean it up with converting to booleans and renaming columns for clarity

ID_COLS = ["HGNC_ID", "ENSG_ID", "NCBI_ID"]

capture_df = (
    pl.read_csv(ALT_ABBREV_ROOT.parent.parent / "output" / "summary_df.csv")
    .rename(
        {
            "captured": "captured_status",
            "captured as:": "captured_category_list",
        }
    )
    .drop("")
    .with_columns(
        pl.when(pl.col("captured_status") == "T")
        .then(True)
        .when(pl.col("captured_status") == "F")
        .then(False)
        .otherwise(None)
        .alias("captured_status"),
        pl.col(ID_COLS).map_elements(
            lambda x: ", ".join(sorted(ast.literal_eval(x))) if x is not None else None,
            return_dtype=pl.String,
        ),
    )
)

# %% [markdown]
# ## Remove samples already in dataset

# %%
eligible_df = capture_df.join(curated_keys, on=dedup_keys, how="anti")

# %% [markdown]
# ## Select sample of new gene-alias pairs

# %%
## To add to the sample set for manual annotation

new_samples_df = eligible_df.group_by("captured_status").map_groups(
    lambda g: g.sample(
        n=min(len(g), NUMBER_OF_NEW_SAMPLES_TO_ADD // 2),
        seed=41,
    )
)

# %%
remaining = NUMBER_OF_NEW_SAMPLES_TO_ADD - new_samples_df.height

if remaining > 0:
    extra = eligible_df.join(
        new_samples_df.select(dedup_keys).unique(),
        on=dedup_keys,
        how="anti",
    ).sample(n=remaining, seed=41)

    new_samples_df = pl.concat([new_samples_df, extra])

# %% [markdown]
# ## Add gene names to the new samples of gene-alias pairs

# %%
# Add gene_name column safely to new samples df

if "gene_name" not in new_samples_df.columns:
    new_samples_df = new_samples_df.with_columns(pl.lit(None).alias("gene_name"))

# %%
# Add gene names, by HGNC ID, to the annotation set for easier manual review

missing_name_ids = (
    new_samples_df.filter(pl.col("gene_name").is_null())
    .select("HGNC_ID")
    .unique()
    .to_series()
    .to_list()
)

gene_map = {hgnc_id: get_gene_name(hgnc_id) for hgnc_id in tqdm(missing_name_ids)}

new_samples_df = new_samples_df.with_columns(
    pl.col("gene_name").fill_null(pl.col("HGNC_ID").replace(gene_map))
)

# %% [markdown]
# ## Combine new samples with those already annotated

# %%
# Add blank column for manual annotation
new_samples_df = new_samples_df.with_columns(
    pl.lit(None).alias("alternate_abbreviation_status")
)
# Reorder columns to match curated set for concatenation
new_samples_df = new_samples_df.select(curated_df.columns)

# %%
df = pl.concat(
    [curated_df, new_samples_df],
    how="vertical",
)

# %% [markdown]
# ## Export newly combined file for manually annotation in Excel

# %%
df.write_excel(ALT_ABBREV_OUTPUT_PATH / "alt_abbrev_annotation_to_annotate_df.xlsx")
f"The original set had {curated_df.height} samples, the new set has {df.height} samples after adding {new_samples_df.height} new samples"
# TODO: After manually annotating the new rows in alt_abbrev_annotation_to_annotate_df.xlsx, it needs to be renamed and saved to alt_abbrev_annotation_**manually_annotated**_df.xlsx


