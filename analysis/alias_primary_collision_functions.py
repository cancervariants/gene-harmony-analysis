"""Functions used in the Alias-Primary Collision Analysis Notebook"""

import numpy as np
import pandas as pd


def create_ap_collision_df(mini_xxxx_df: pd.DataFrame, source: str, case_sensitive: bool = False) -> pd.DataFrame:
    """Create a DataFrame of alias-primary collision symbols.
    :param mini_xxxx_df: Processed df of gene records (blank values converted to nans)
    :param source: Representation of the gene source ("HGNC", "NCBI", "ENSG")
    :param case_sensitive: Whether symbol matching is case-sensitive
    :return: A df of genes that share an alias with another gene"s primary gene symbol
    """
    # Remove placeholder gene records
    mini_xxxx_df = mini_xxxx_df.dropna(subset=["primary_gene_symbol"])
    mini_xxxx_df.to_csv(f"../output/mini_{source.lower()}_df.csv", index=False)

    # Create merged alias df (for later verification)
    merged_alias_xxxx_df = (
        mini_xxxx_df
        .copy()
        .fillna("")
        .groupby(f"{source}_ID", dropna=False)
        .agg(lambda col: ",".join(sorted(set(map(str, col)))))
        .reset_index()
    )
    merged_alias_xxxx_df.to_csv(f"../output/merged_alias_{source}_df.csv", index=False)

    # Case normalization logic
    def normalize(val: object) -> object:
        return val if case_sensitive or pd.isna(val) else str(val).upper()

    # Normalize primary gene symbols
    mini_xxxx_df["normalized_primary"] = mini_xxxx_df["primary_gene_symbol"].map(normalize)
    xxxx_gene_symbol_set = set(mini_xxxx_df["normalized_primary"])

    subset_genes_xxxx_df = mini_xxxx_df.copy()

    # Normalize alias and remove if alias matches primary
    subset_genes_xxxx_df["normalized_alias"] = subset_genes_xxxx_df["alias_symbol"].map(normalize)
    subset_genes_xxxx_df["alias_symbol"] = subset_genes_xxxx_df.apply(
        lambda row: np.nan if row["normalized_alias"] == row["normalized_primary"] else row["alias_symbol"], axis=1
    )
    subset_genes_xxxx_df = subset_genes_xxxx_df.drop(columns=["normalized_alias", "normalized_primary"])
    subset_genes_xxxx_df.to_csv(f"../output/subset_genes_{source}_df.csv", index=False)

    # Create alias-primary collision DataFrame
    ap_collision_xxxx_df = subset_genes_xxxx_df.dropna(subset=["alias_symbol"]).copy()
    ap_collision_xxxx_df["normalized_alias"] = ap_collision_xxxx_df["alias_symbol"].map(normalize)
    ap_collision_xxxx_df["alias_symbol_set"] = ap_collision_xxxx_df["normalized_alias"].apply(lambda x: {x})

    ap_collision_xxxx_df["collision"] = ap_collision_xxxx_df["alias_symbol_set"].apply(
        lambda x: x & xxxx_gene_symbol_set
    )
    ap_collision_xxxx_df = ap_collision_xxxx_df[ap_collision_xxxx_df["collision"].apply(lambda x: len(x) > 0)]
    ap_collision_xxxx_df["collision"] = ap_collision_xxxx_df["collision"].apply(lambda x: ", ".join(x))
    ap_collision_xxxx_df = ap_collision_xxxx_df.sort_values("collision")
    ap_collision_xxxx_df["source"] = str(source.upper())
    ap_collision_xxxx_df.drop(columns=["normalized_alias", "alias_symbol_set"], inplace=True)
    ap_collision_xxxx_df.to_csv(f"../output/single_alias_ap_collision_{source.lower()}_df.csv", index=False)

    # Merge alias list back for verification
    columns_map = {
        "ENSG": ["NCBI_ID", "HGNC_ID"],
        "HGNC": ["NCBI_ID", "ENSG_ID"],
        "NCBI": ["HGNC_ID", "ENSG_ID"]
    }
    cols_of_interest = columns_map.get(source, [])

    merged_alias_ap_collision_xxxx_df = ap_collision_xxxx_df.drop(
        columns=[*cols_of_interest, "alias_symbol"],
        errors="ignore",
    )

    merged_alias_ap_collision_xxxx_df = pd.merge(
        merged_alias_ap_collision_xxxx_df,
        merged_alias_xxxx_df[
            [f"{source}_ID", *cols_of_interest, "alias_symbol"]
        ],
        on=f"{source}_ID",
        how="left"
    ).drop_duplicates(subset=[f"{source}_ID"])

    # Test if all collisions are present in alias list
    test_df = merged_alias_ap_collision_xxxx_df.copy()
    test_df["alias_symbol_list"] = test_df["alias_symbol"].str.split(",")
    test_df["alias_symbol_list"] = test_df["alias_symbol_list"].apply(lambda x: [normalize(i.strip()) for i in x] if isinstance(x, list) else [])

    test_df["collision_list"] = test_df["collision"].str.split(", ")
    test_df["collision_in_alias"] = test_df.apply(
        lambda row: all(c in row["alias_symbol_list"] for c in row["collision_list"]),
        axis=1
    )

    if test_df["collision_in_alias"].all():
        merged_alias_ap_collision_xxxx_df.to_csv(f"../output/merged_alias_ap_collision_{source.lower()}_df.csv", index=False)
        print("All collisions are present in gene alias lists.")
    else:
        print("Some collisions are not present in gene alias lists.")

    return mini_xxxx_df.head()
