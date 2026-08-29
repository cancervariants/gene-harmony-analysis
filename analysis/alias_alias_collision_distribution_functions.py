"""Functions used in the Alias-Alias Collision Distribution Analysis Notebook"""

import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure


def create_aa_collision_histogram(
    aa_collision_gene_xxxx_df: pd.DataFrame,
    source: str,
    xxxx_aa_collision_count: int,
) -> Figure:
    """Create a histogram of the frequencies at which aliases are shared

    :param aa_collision_gene_xxxx_df: A df of alias-alias collisions organized by primary gene symbol
    :param source: Representation of the source of the gene records
    :param xxxx_alias_count: Number of aliases total in the source
    :return: A histogram of the percentage of aliases that are shared between 2 genes, 3 genes, and so on
    """
    #Count the number of times each shared alias is used
    aa_collision_xxxx_count_df = aa_collision_gene_xxxx_df.pivot_table(
    index=["collision"], aggfunc="size"
    )
    aa_collision_xxxx_count_df = aa_collision_xxxx_count_df.reset_index()
    aa_collision_xxxx_count_df.rename(columns={0: "num_gene_records"}, inplace=True)
    aa_collision_xxxx_count_df = aa_collision_xxxx_count_df.sort_values(
        "num_gene_records", ascending=False)

    #Convert to csv
    aa_collision_xxxx_count_df.to_csv(f"../output/aa_collision_{source}_count_df.csv", index=True)

    #Count the frequency at which aliases are shared 
    aa_collision_xxxx_distribution_df = aa_collision_xxxx_count_df.pivot_table(
    index=["num_gene_records"], aggfunc="size"
    )
    aa_collision_xxxx_distribution_df = aa_collision_xxxx_distribution_df.reset_index()
    aa_collision_xxxx_distribution_df.rename(columns={0: "num_collision_symbol"}, inplace=True)
    aa_collision_xxxx_distribution_df["percent_collision_symbol"] = (
        aa_collision_xxxx_distribution_df["num_collision_symbol"] / xxxx_aa_collision_count
    ) * 100

    #Convert to csv
    aa_collision_xxxx_distribution_df.to_csv(f"../output/aa_collision_{source}_distribution_df.csv", index=True)

    #Create histogram df 
    xxxx_alias_count_histogram_df = aa_collision_xxxx_distribution_df.drop(
    "num_collision_symbol", axis=1)

    #Convert to csv
    xxxx_alias_count_histogram_df.to_csv(f"../output/{source}_alias_count_histogram_df.csv", index=True)

    return px.bar(xxxx_alias_count_histogram_df, x="num_gene_records", y="percent_collision_symbol")

