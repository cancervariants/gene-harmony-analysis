"""Functions used in the Alias-Alias Collision Analysis Notebook"""

import pandas as pd


def create_aa_collision_df(subset_genes_xxxx_df: pd.DataFrame, merged_alias_xxxx_df: pd.DataFrame, source: str, case_sensitive: bool = False) -> pd.DataFrame:
    """Create a DataFrame of alias-alias collision symbols.

    :param subset_genes_xxxx_df: Processed df of gene records (one alias per row, no records without aliases)
    :param merged_alias_xxxx_df: Merged alias df for post-processing enrichment
    :param source: Representation of the gene source ("HGNC", "NCBI", "ENSG")
    :param case_sensitive: Whether alias matching is case-sensitive
    :return: A df of genes that share an alias with another gene
    """
    # Normalize alias symbols if not case-sensitive
    if not case_sensitive:
        subset_genes_xxxx_df["alias_symbol_upper"] = subset_genes_xxxx_df["alias_symbol"].str.upper()
        alias_col = "alias_symbol_upper"
    else:
        alias_col = "alias_symbol"

    # Find alias symbols used by multiple different genes
    dup_alias = (
        subset_genes_xxxx_df.groupby(alias_col)["primary_gene_symbol"]
        .nunique()
        .reset_index()
    )
    dup_alias = dup_alias[dup_alias["primary_gene_symbol"] > 1][alias_col]

    # Filter for rows with those duplicate aliases
    aa_collision_xxxx_df = subset_genes_xxxx_df[
        subset_genes_xxxx_df[alias_col].isin(dup_alias)
    ].copy()

    # Rename column to "collision" and clean up temp column
    aa_collision_xxxx_df = aa_collision_xxxx_df.rename(columns={"alias_symbol": "collision"})
    if not case_sensitive:
        aa_collision_xxxx_df = aa_collision_xxxx_df.drop("alias_symbol_upper", axis=1)

    # Sort and tag by source
    aa_collision_xxxx_df = aa_collision_xxxx_df.sort_values("collision")
    aa_collision_xxxx_df["source"] = str(source)

    # Save single alias collision output
    aa_collision_xxxx_df.to_csv(f"../output/single_alias_aa_collision_{source.lower()}_df.csv", index=False)

    # Prepare for merged output
    columns_map = {
        "ENSG": ["NCBI_ID", "HGNC_ID"],
        "HGNC": ["NCBI_ID", "ENSG_ID"],
        "NCBI": ["HGNC_ID", "ENSG_ID"]
    }
    cols_of_interest = columns_map.get(source, [])

    merged_alias_aa_collision_xxxx_df = aa_collision_xxxx_df.drop(columns=cols_of_interest)

    merged_alias_aa_collision_xxxx_df = pd.merge(
        merged_alias_aa_collision_xxxx_df,
        merged_alias_xxxx_df[
            [f"{source}_ID", *cols_of_interest, "alias_symbol"]
        ],
        on=f"{source}_ID",
        how="left"
    )

    # Save merged alias output
    merged_alias_aa_collision_xxxx_df.to_csv(
        f"../output/merged_alias_aa_collision_{source.lower()}_df.csv", index=False
    )

    return merged_alias_aa_collision_xxxx_df.head()
