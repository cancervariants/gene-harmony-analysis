import marimo

__generated_with = "0.23.9"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # **Automatic Alternate Abbreviation Annotation:** Adding and processing gene-alias pairs
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Running this notebook is optional. It adds samples to the manually annotated Alternate Abbreviations dataset using a partially random selection process that preserves the existing proportion of captured versus non-captured alias symbols. For more information on Alternate Abbreviations and workflow please see [README](alternate_abbreviation_README).
    """)
    return


@app.cell
def _():
    import ast
    import time
    from pathlib import Path

    import polars as pl
    import requests
    from tqdm.notebook import tqdm

    return Path, ast, pl, requests, time, tqdm


@app.cell
def _(requests, time):
    def get_gene_name(hgnc_id: str) -> str:
        """Retrieve official gene name from HGNC for each sample in the set.

        :param hgnc_id: The HGNC ID to retrieve the official gene name for
        :return: The official HGNC gene name or an error message
        """
        url = f"https://rest.genenames.org/fetch/hgnc_id/{hgnc_id}"
        headers = {"Accept": "application/json"}  # Request JSON format
        time.sleep(0.5)  # Sleep to avoid hitting API rate limits
        response = requests.get(url, headers=headers, timeout=10)

        if response.status_code != 200:
            return f"Error: {response.status_code}"

        data = response.json()

        try:
            return data["response"]["docs"][0]["name"]
        except IndexError:
            return f"No gene found for HGNC ID {hgnc_id}"

    return (get_gene_name,)


@app.cell
def _(Path):
    ALT_ABBREV_ROOT = Path.cwd().resolve()
    ALT_ABBREV_OUTPUT_PATH = ALT_ABBREV_ROOT / "output"
    GENE_HARMONY_OUTPUT_PATH = ALT_ABBREV_ROOT.parent.parent / "output"
    return ALT_ABBREV_OUTPUT_PATH, GENE_HARMONY_OUTPUT_PATH


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Define the number of new samples to add to the existing dataset
    """)
    return


@app.cell
def _():
    NUMBER_OF_NEW_SAMPLES_TO_ADD = 2
    return (NUMBER_OF_NEW_SAMPLES_TO_ADD,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load dataset with gene-alias pairs manually curated for Alternate Abbreviation alias symbols
    """)
    return


@app.cell
def _(ALT_ABBREV_OUTPUT_PATH, pl):
    curated_df = pl.read_excel(
        ALT_ABBREV_OUTPUT_PATH / "alt_abbrev_annotation_manually_annotated_df.xlsx"
    )
    return (curated_df,)


@app.cell
def _(curated_df):
    dedup_keys = ["primary_gene_symbol", "gene_symbol"]
    curated_keys = curated_df.select(dedup_keys).unique()
    return curated_keys, dedup_keys


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load and clean dataset where to get new gene-alias pairs from
    """)
    return


@app.cell
def _(GENE_HARMONY_OUTPUT_PATH, ast, pl):
    ## This is the file with alias and primary gene symbol pairs and what categories the aliases were captured as
    ## Generated in the 5_symbol_capture_analysis.ipynb
    ## Clean it up with converting to booleans and renaming columns for clarity

    ID_COLS = ["HGNC_ID", "ENSG_ID", "NCBI_ID"]

    capture_df = (
        pl.read_csv(GENE_HARMONY_OUTPUT_PATH / "summary_df.csv")
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
    return (capture_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Remove samples already in dataset
    """)
    return


@app.cell
def _(capture_df, curated_keys, dedup_keys):
    eligible_df = capture_df.join(curated_keys, on=dedup_keys, how="anti")
    return (eligible_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Select sample of new gene-alias pairs
    """)
    return


@app.cell
def _(NUMBER_OF_NEW_SAMPLES_TO_ADD, eligible_df):
    ## To add to the sample set for manual annotation

    new_samples_df = eligible_df.group_by("captured_status").map_groups(
        lambda g: g.sample(
            n=min(len(g), NUMBER_OF_NEW_SAMPLES_TO_ADD // 2),
            seed=41,
        )
    )
    return (new_samples_df,)


@app.cell
def _(
    NUMBER_OF_NEW_SAMPLES_TO_ADD,
    dedup_keys,
    eligible_df,
    new_samples_df,
    pl,
):
    remaining = NUMBER_OF_NEW_SAMPLES_TO_ADD - new_samples_df.height
    if remaining > 0:
        extra = eligible_df.join(new_samples_df.select(dedup_keys).unique(), on=dedup_keys, how='anti').sample(n=remaining, seed=41)
        new_samples_df_1 = pl.concat([new_samples_df, extra])
    return (new_samples_df_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Add gene names to the new samples of gene-alias pairs
    """)
    return


@app.cell
def _(new_samples_df_1, pl):
    # Add gene_name column safely to new samples df
    if 'gene_name' not in new_samples_df_1.columns:
        new_samples_df_2 = new_samples_df_1.with_columns(pl.lit(None).alias('gene_name'))
    return (new_samples_df_2,)


@app.cell
def _(get_gene_name, new_samples_df_2, pl, tqdm):
    # Add gene names, by HGNC ID, to the annotation set for easier manual review
    missing_ids = new_samples_df_2.filter(pl.col('gene_name').is_null()).select('HGNC_ID').unique().to_series().to_list()
    gene_map = {hgnc_id: get_gene_name(hgnc_id) for hgnc_id in tqdm(missing_ids)}
    new_samples_df_3 = new_samples_df_2.with_columns(pl.col('gene_name').fill_null(pl.col('HGNC_ID').replace(gene_map)))
    return (new_samples_df_3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Combine new samples with those already annotated
    """)
    return


@app.cell
def _(curated_df, new_samples_df_3, pl):
    # Add blank column for manual annotation
    new_samples_df_4 = new_samples_df_3.with_columns(pl.lit(None).alias('alternate_abbreviation_status'))
    # Reorder columns to match curated set for concatenation
    new_samples_df_4 = new_samples_df_4.select(curated_df.columns)
    return (new_samples_df_4,)


@app.cell
def _(curated_df, new_samples_df_4, pl):
    df = pl.concat([curated_df, new_samples_df_4], how='vertical')
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Export newly combined file for manually annotation in Excel
    """)
    return


@app.cell
def _(ALT_ABBREV_OUTPUT_PATH, curated_df, df, new_samples_df_4):
    df.write_excel(ALT_ABBREV_OUTPUT_PATH / 'alt_abbrev_annotation_to_annotate_df.xlsx')
    f'The original set had {curated_df.height} samples, the new set has {df.height} samples after adding {new_samples_df_4.height} new samples'
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### After manually annotating the new rows in alt_abbrev_annotation_to_annotate_df.xlsx, it needs to be renamed and saved to alt_abbrev_annotation_**manually_annotated**_df.xlsx
    """)
    return


if __name__ == "__main__":
    app.run()

