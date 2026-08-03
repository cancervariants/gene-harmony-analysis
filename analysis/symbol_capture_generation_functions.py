"""Functions used in the Symbol Capture Generation Notebook

The sections that the functions are used in are highlighted in this file as headings
"""
import ast
import re
import time
from collections import deque
from collections.abc import Hashable
from io import StringIO

import nest_asyncio
import pandas as pd
import requests
from pybiomart import Server
from requests.exceptions import RequestException
from tqdm import tqdm

tqdm.pandas()

nest_asyncio.apply()


# Download gene records from ENSG, HGNC, and NCBI
def read_subset_genes_csv(location: str) -> pd.DataFrame:
    """Create a df of primary gene symbol- alias symbol pairs

    :param location: file location
    return: a df of gene records
    """
    subset_genes_xxxx_df = pd.read_csv(
        location, dtype={"NCBI_ID": str,"HGNC_ID":str})
    return subset_genes_xxxx_df

def remove_nan_from_set(s: set[Hashable]) -> set[Hashable]:
    """Remove null instances from set

    :param s: selected set
    :return: set with no null values
    """
    return {x for x in s if pd.notna(x)}

# Combine all three resources
def combine_rows(row: pd.Series) -> tuple[bool, list[str] | None]:
    """Combine identifier matches and their associated sources from a row.

    :param row: a DataFrame row containing gene identifier match columns and
                their associated source columns.
    :return: a tuple containing:
             - a boolean indicating whether any gene identifier match exists.
             - a list of unique associated sources for matching identifiers, or
               None if no matches are found.
    """
    a1, a2, a3 = row.get("Gene Identifier Match_1"), row.get("Gene Identifier Match_2"), row.get("Gene Identifier Match_3")
    b1, b2, b3 = row.get("matching_attribute"), row.get("associated_source"), row.get("associated_source_3")

    combined = []
    match = False

    for a, b in [(a1, b1), (a2, b2), (a3, b3)]:
        if a:
            match = True
            if isinstance(b, list):
                combined.extend(b)
            elif b and b != "<NA>":
                combined.append(b)

    if not combined:
        return (False, None)

    return (match, list(set(combined)))

# Ortholog Symbol Capture
def make_col_ortholog_match(
    recording_df: pd.DataFrame,
    source_df: pd.DataFrame,
    animal: str,
) -> pd.DataFrame:
    """Check for ortholog matches in the primary gene symbol- alias symbol pairs.
    Adds a T/F column for each pair. T if the alias is an ortholog from the specified animal and F if not

    :param recording_df: df that contains the primary gene symbol- alias symbol pairs
    :param source_df: df that contains the orthologs and their associated human genes
    :param animal: the animal from with the orthologs are being checked
    return: the number of primary gene symbol- alias symbol pairs where the alias is an ortholog from the specified animal
    """
    source_df["Gene name upper"] = source_df["Gene name"].str.upper()
    source_df[f"{animal} gene name upper"] = source_df[f"{animal} gene name"].str.upper()

    # Build a set of (gene, ortholog) pairs for this animal
    ortholog_pairs = set(zip(source_df["Gene name upper"], source_df[f"{animal} gene name upper"]))
    source_df = source_df.drop(["Gene name upper", f"{animal} gene name upper"], axis=1)

    # Match using vectorized lookup
    df = recording_df.copy()
    df["primary_gene_symbol_upper"] = df["primary_gene_symbol"].str.upper()
    df["alias_symbol_upper"] = df["alias_symbol"].str.upper()
    df[f"{animal} Match"] = df.apply(
        lambda row: (row["primary_gene_symbol_upper"], row["alias_symbol_upper"],) in ortholog_pairs,
        axis=1
    )
    df = df.drop(["primary_gene_symbol_upper", "alias_symbol_upper"], axis=1)
    print(f"Added column: {animal} Match")
    return df

def match_alias_to_ortholog(
    og_recording_df: pd.DataFrame,
    source_df: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Apply the make_col_ortholog_match function to all animal columns in the DataFrame.

    :param og_recording_df: DataFrame containing the primary gene symbol- alias symbol pairs
    :param source_df: DataFrame containing the orthologs and their associated human genes
    :return: a DataFrame with all match columns added
    """
    source_df = source_df.groupby("Gene name", as_index=False).agg(combine_rows)

    recording_df = og_recording_df.copy()
    recording_df.columns = recording_df.columns.str.strip().str.replace(r"\s+", " ", regex=True)
    source_df.columns = source_df.columns.str.strip().str.replace(r"\s+", " ", regex=True)

    animal_columns = [col for col in source_df.columns if "gene name" in col and col != "Gene name"]

    true_counts = {}

    for animal in animal_columns:
        # Extract the animal name
        animal_name = animal.replace(" gene name", "")

        # Keep only rows with values in both columns
        source_df = source_df.dropna(subset= animal)

        recording_df = make_col_ortholog_match(recording_df, source_df, animal_name)

        true_count = recording_df[f"{animal_name} Match"].sum()
        true_counts[animal_name] = true_count
    return recording_df, true_counts

# Disorder/Disease Symbol Capture
def split_name_symbol_pairs(value: str) -> list[list[str]]:
    """Split a string of name-symbol pairs into a list of pairs.

    :param value: string containing name-symbol pairs separated by ``;;``,
                  with each pair separated by ``;``.
    :return: a list of two-element lists, where each inner list contains a
             name and its associated symbol.
    """
    pairs = value.split(";;")

    result = []
    for pair in pairs:
        parts = pair.strip().split(";")
        if len(parts) == 2:
            result.append(parts)
    return result


def split_phenotype_on_last_comma(phenotype: str) -> tuple[str, str]:
    """Split a phenotype string at the last comma.

    :param phenotype: phenotype string containing a description and location.
    :return: a tuple containing the text before the last comma and the text
             after the last comma. If no comma is present, the second value is
             an empty string.
    """
    parts = phenotype.rsplit(",", 1)
    return parts[0], parts[1] if len(parts) > 1 else ""


def split_phenotype_description(description: str) -> tuple[str, str]:
    """Split a phenotype description from its OMIM identifier.

    :param description: phenotype description that may end with an OMIM
                        identifier in parentheses.
    :return: a tuple containing the phenotype description and the OMIM
             identifier. If no identifier is present, the second value is an
             empty string.
    """
    match = re.match(r"(.*?)(\s*\(\d+\))?$", description)
    if match:
        return match.group(1).strip(), match.group(2).strip() if match.group(2) else ""
    else:
        return description, ""


def check_prefix_disease_match(
    row: pd.Series,
    pheno_map: dict[str, list[str]],
) -> bool:
    """Check whether an alias begins with a phenotype symbol.

    :param row: a DataFrame row containing the uppercase alias symbol and
                primary gene symbol.
    :param pheno_map: mapping of primary gene symbols to phenotype symbols.
    :return: True if the alias begins with a phenotype symbol associated with
             the primary gene symbol; otherwise, False.
    """
    alias = row["alias_symbol_upper"]
    primary_gene_symbol = row["primary_gene_symbol_upper"]

    if not isinstance(alias, str):
        return False

    phenotypes = pheno_map.get(primary_gene_symbol, [])
    return any(alias.startswith(phenotype) for phenotype in phenotypes)


def get_matching_phenotype_symbol_fast(
    row: pd.Series,
    pheno_map: dict[str, list[str]],
) -> str:
    """Find the matching phenotype prefix for an alias.

    :param row: a DataFrame row containing the uppercase alias symbol and
                primary gene symbol.
    :param pheno_map: mapping of primary gene symbols to phenotype symbols.
    :return: the matching phenotype prefix if the alias begins with one;
             otherwise, an empty string.
    """
    alias = row["alias_symbol_upper"]
    primary_gene_symbol = row["primary_gene_symbol_upper"]

    if not isinstance(alias, str):
        return ""

    phenotypes = pheno_map.get(primary_gene_symbol, [])
    return next(
        (symbol for symbol in phenotypes if alias.startswith(symbol)),
        "",
    )

# Previous Symbol Capture
def match_alias_to_previous_symbol(
    analysis_df: pd.DataFrame,
    sources: list[tuple[pd.DataFrame, str, str, str, str]],
) -> pd.DataFrame:
    """Identify previous symbol matches across one or more reference sources.

    :param analysis_df: DataFrame containing alias symbols and identifier columns.
    :param sources: List of tuples containing
        (reference_df, analysis_id_col, reference_id_col,
        reference_symbol_col, source_name).
    :return: DataFrame containing the columns
        ``Previous Symbol Match`` and ``Previous Symbol Source``.
    """
    # Precompute lookup sets once for each source
    source_lookup = []
    for (
        reference_df,
        analysis_id_col,
        reference_id_col,
        reference_symbol_col,
        source_name,
    ) in sources:
        reference_pairs = set(
            zip(
                reference_df[reference_id_col].str.upper(),
                reference_df[reference_symbol_col].str.upper(),
            )
        )

        source_lookup.append(
            (
                analysis_id_col,
                reference_pairs,
                source_name,
            )
        )

    matches = []
    source_labels = []

    for _, row in tqdm(
        analysis_df.iterrows(),
        total=len(analysis_df),
        desc="Matching previous symbols",
    ):
        alias = row["alias_symbol"]

        if pd.isna(alias):
            matches.append(False)
            source_labels.append(pd.NA)
            continue

        alias = alias.upper()
        row_sources = []

        for analysis_id_col, reference_pairs, source_name in source_lookup:
            ids = row[analysis_id_col]

            if isinstance(ids, (list, set)):
                if any(
                    (gene_id.upper(), alias) in reference_pairs
                    for gene_id in ids
                    if isinstance(gene_id, str)
                ):
                    row_sources.append(source_name)

        matches.append(bool(row_sources))
        source_labels.append(", ".join(row_sources) if row_sources else pd.NA)

    return pd.DataFrame(
        {
            "Previous Symbol Match": matches,
            "Previous Symbol Source": source_labels,
        },
        index=analysis_df.index,
    )

# Gene Family Symbol Capture
def expand_gene_groups(
    gene_df: pd.DataFrame,
    hierarchy_df: pd.DataFrame,
) -> pd.DataFrame:
    """Expand gene group IDs to include parent and grandparent groups.

    :param gene_df: DataFrame containing HGNC IDs, approved symbols, and gene
                    group IDs.
    :param hierarchy_df: DataFrame containing parent and child gene group IDs.
    :return: DataFrame containing the original gene groups and all associated
             parent gene groups.
    """
    child_to_parents = (
        hierarchy_df.groupby("child_fam_id")["parent_fam_id"]
        .apply(list)
        .to_dict()
    )

    expanded_rows = []
    queue = deque()
    seen = set()

    for _, row in gene_df.iterrows():
        hgnc_id = row["HGNC ID"]
        approved_symbol = row["Approved symbol"]
        group_id = row["Gene group ID"]

        queue.append((hgnc_id, approved_symbol, group_id))
        seen.add((hgnc_id, group_id))
        expanded_rows.append(
            {
                "HGNC ID": hgnc_id,
                "Approved symbol": approved_symbol,
                "Gene group ID": group_id,
            }
        )

    while queue:
        hgnc_id, symbol, group_id = queue.popleft()

        parent_ids = child_to_parents.get(group_id, [])
        for parent_id in parent_ids:
            key = (hgnc_id, parent_id)

            if key not in seen:
                seen.add(key)
                queue.append((hgnc_id, symbol, parent_id))
                expanded_rows.append(
                    {
                        "HGNC ID": hgnc_id,
                        "Approved symbol": symbol,
                        "Gene group ID": parent_id,
                    }
                )

    return pd.DataFrame(expanded_rows)


def get_matching_abbreviation(
    row: pd.Series,
    abbrev_map: dict[str, list[str]],
) -> str:
    """Find the gene group abbreviation that matches an alias prefix.

    :param row: DataFrame row containing uppercase alias and primary gene
                symbols.
    :param abbrev_map: mapping of primary gene symbols to abbreviations.
    :return: the first matching abbreviation, or an empty string if no match
             is found.
    """
    alias = row["alias_symbol_upper"]
    primary_gene_symbol = row["primary_gene_symbol_upper"]

    if not isinstance(alias, str):
        return ""

    abbreviations = abbrev_map.get(primary_gene_symbol, [])

    for abbreviation in abbreviations:
        if alias.startswith(abbreviation):
            return abbreviation

    return ""

# Protein Mass Symbol Capture
def get_uniprot_data(
    uniprot_ids: list[str],
    fields: list[str] | None = None,
) -> pd.DataFrame:
    """Retrieve UniProt attributes for a list of UniProt accessions.

    :param uniprot_ids: list of UniProt accession identifiers to query.
    :param fields: list of UniProt fields to retrieve. If None, retrieves the
                   accession and protein mass.
    :return: a DataFrame containing the requested UniProt attributes. Returns
             an empty DataFrame if the request fails.
    """
    if fields is None:
        fields = ["accession", "mass"]

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = " OR ".join(f"accession:{uid}" for uid in uniprot_ids)
    field_str = ",".join(fields)

    params = {
        "query": query,
        "fields": field_str,
        "format": "tsv",
        "size": len(uniprot_ids),
    }

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        df = pd.read_csv(StringIO(response.text), sep="\t")
        return df
    else:
        print(f"Failed chunk with status {response.status_code}")
        return pd.DataFrame()


def digit_prefix_matches_with_alias(
    row: pd.Series,
    df2: pd.DataFrame,
    uniprot_associated_gene_id_col: str,
    subset_df_database_id_col: str,
) -> tuple[bool, str | None]:
    """Check whether a protein mass matches the numeric prefix of an alias.

    :param row: a DataFrame row containing UniProt protein mass information.
    :param df2: DataFrame containing alias symbols and associated gene
                identifiers.
    :param uniprot_associated_gene_id_col: name of the column containing the
                                           UniProt-associated gene identifier.
    :param subset_df_database_id_col: name of the column containing the gene
                                      identifier collection in ``df2``.
    :return: a tuple containing a boolean indicating whether a match was found
             and the matching alias symbol. Returns ``None`` for the alias if
             no match is found.
    """
    val1 = row[uniprot_associated_gene_id_col]
    val2_df1 = str(row["Mass"])
    matching_rows = df2[df2[subset_df_database_id_col].apply(lambda s: val1 in s)]
    pattern = re.compile(r"^[Pp](\d{1,3})")

    for _, r2 in matching_rows.iterrows():
        alias = str(r2["alias_symbol"])
        m = pattern.match(alias)
        if m:
            digits = m.group(1)
            if val2_df1.startswith(digits):
                return True, alias

    return False, None


def alias_match_for_row(
    row: pd.Series,
    source_df: pd.DataFrame,
    subset_df_database_id_col: str,
    uniprot_associated_gene_id_col: str,
) -> bool:
    """Check whether an alias matches a UniProt-associated gene identifier.

    :param row: a DataFrame row containing an alias symbol and associated gene
                identifiers.
    :param source_df: DataFrame containing matched alias symbols and
                      UniProt-associated gene identifiers.
    :param subset_df_database_id_col: name of the column containing the gene
                                      identifier collection.
    :param uniprot_associated_gene_id_col: name of the column containing the
                                           UniProt-associated gene identifier.
    :return: True if the alias symbol matches a record in ``source_df`` and
             the associated gene identifier is present; otherwise, False.
    """
    alias_symbol_raw = row["alias_symbol"]
    if pd.isna(alias_symbol_raw):
        return False

    alias_symbol = alias_symbol_raw.casefold()
    ens_set = row[subset_df_database_id_col]

    if not isinstance(ens_set, (set, list)):
        return False

    matched_rows = source_df[
        source_df["Matched Alias Symbol"].str.casefold() == alias_symbol
    ]

    for _, r2 in matched_rows.iterrows():
        if r2[uniprot_associated_gene_id_col] in ens_set:
            return True

    return False

# Gene Neighbor Symbol Capture
def optimized_split_and_explode(
    df: pd.DataFrame,
    col: str,
) -> pd.DataFrame:
    """Split pipe-delimited values in a column and expand them into separate rows.

    :param df: DataFrame containing the column to split and expand.
    :param col: name of the column containing pipe-delimited values.
    :return: a DataFrame in which each pipe-delimited value has been expanded
             into its own row. Rows without pipe-delimited values are retained
             unchanged.
    """
    has_pipe = df[col].str.contains(r"\|", na=False)
    df_to_explode = df[has_pipe].copy()
    df_to_explode[col] = df_to_explode[col].str.split("|")
    exploded = df_to_explode.explode(col)

    # Add back rows that didn't need exploding
    df_no_pipe = df[~has_pipe]
    return pd.concat([df_no_pipe, exploded], ignore_index=True)

# Gene Identifier Symbol Capture
def safe_query(
    dataset,
    attributes: list[str],
    filters: dict[str, list[str]],
    retries: int = 3,
    wait: int = 5,
) -> pd.DataFrame:
    """Query BioMart with retry logic.

    :param dataset: BioMart dataset to query.
    :param attributes: list of BioMart attributes to retrieve.
    :param filters: dictionary of query filters.
    :param retries: maximum number of retry attempts.
    :param wait: seconds to wait between retries.
    :return: DataFrame containing the query results.
    """
    for attempt in range(1, retries + 1):
        try:
            return dataset.query(attributes=attributes, filters=filters)
        except (
            ConnectionResetError,
            ConnectionAbortedError,
            RequestException,
            TimeoutError,
        ) as e:
            print(f"Connection error on attempt {attempt}/{retries}: {e}")
            if attempt == retries:
                raise
            time.sleep(wait * attempt)
        except Exception as e:
            print(f"Unexpected error on attempt {attempt}/{retries}: {e}")
            if attempt == retries:
                raise
            time.sleep(wait * attempt)

    raise RuntimeError("safe_query exited without returning a result")

def get_attributes_for_ensembl_ids(
    ensembl_ids: list[str] | set[str],
    attributes: list[str],
    dataset=None,
) -> dict[str, dict[str, list[str]]]:
    """Query Ensembl BioMart for MULTIPLE Ensembl gene IDs and return attribute values.

    :param ensembl_ids: Set or list of Ensembl gene IDs (e.g., {"ENSG00000141510", "ENSG00000155657"}).
    :param attributes: List of attribute keys to retrieve (e.g., ["external_gene_name", "hgnc_symbol"]).
    :param dataset: (pybiomart.dataset.Dataset, optional) Pass a pre-initialized dataset object to reuse the connection.

    :return: dict in the form {ensembl_id: {attribute_key: [values]}} 
             where each value is a list of unique non-null entries returned from BioMart.
    """
    if dataset is None:
        server = Server(host="http://www.ensembl.org")
        dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

    display_names = {attr.name: attr.display_name for attr in dataset.attributes.values()}

    result = {ensembl_id: {attr: [] for attr in attributes} for ensembl_id in ensembl_ids}

    for attr in attributes:
        display_name = display_names.get(attr, attr)

        try:
            query_result = safe_query(
                dataset,
                attributes=["ensembl_gene_id", attr],
                filters={"link_ensembl_gene_id": list(ensembl_ids)}
            )

            if display_name in query_result.columns and "Gene stable ID" in query_result.columns:
                grouped = query_result.dropna().groupby("Gene stable ID")[display_name]

                for ensembl_id, values in grouped:
                    if ensembl_id in result:
                        result[ensembl_id][attr] = list(values.dropna().unique())
        except Exception as e:
            print(f"Error querying attribute '{attr}': {e}")
            # All Ensembl IDs already default to [] for this attribute

    return result

def batch_query_ensembl(
    ids: list[str] | set[str],
    attributes: list[str],
    dataset,
    batch_size: int = 400,
) -> dict[str, dict[str, list[str]]]:
    """Query Ensembl BioMart in batches to avoid exceeding URL length or server limits.

    :param ids: Iterable of Ensembl gene IDs (e.g., list or set of "ENSG..." strings).
    :param attributes: List of attribute keys to retrieve (e.g., ["external_gene_name", "hgnc_symbol"]).
    :param dataset: (pybiomart.dataset.Dataset) A pre-initialized dataset object for querying BioMart.
    :param batch_size: Number of Ensembl IDs to include per batch (default is 400).

    :return: dict in the form {ensembl_id: {attribute_key: [values]}} with results merged across all batches.
    """
    results = {}
    for i in tqdm(range(0, len(ids), batch_size), desc="Batch querying BioMart"):
        batch = list(ids)[i:i+batch_size]
        try:
            batch_results = get_attributes_for_ensembl_ids(batch, attributes, dataset=dataset)
            results.update(batch_results)
        except Exception as e:
            print(f"Error querying batch {i//batch_size+1}: {e}")
    return results

def match_alias_to_ensembl_attributes(
    df: pd.DataFrame,
    id_column: str,
    alias_column: str,
    attributes: list[str],
    dataset=None,
) -> pd.DataFrame:
    """Match alias values in a DataFrame to Ensembl attribute values using BioMart.

    For each row, checks whether the alias value appears in any of the specified Ensembl attributes
    for any of the associated Ensembl gene IDs. Stops at the first match found (attribute-level or ID-level)
    and returns match details.

    :param df: pandas.DataFrame with columns for Ensembl IDs and alias values.
    :param id_column: Name of the column in `df` containing a set or list of Ensembl gene IDs per row.
    :param alias_column: Name of the column in `df` containing the alias string to match.
    :param attributes: List of attribute keys to query from Ensembl (e.g., ["external_gene_name", "hgnc_symbol"]).
    :param dataset: (pybiomart.dataset.Dataset, optional) Reuse a pre-initialized dataset if provided.

    :return: pandas.DataFrame with additional columns:
        - "match_found": Boolean indicating if alias matched any attribute value.
        - "matching_value": The alias value that matched (or None).
        - "matching_attribute": Display name of the attribute that matched (or None).
    """
    if dataset is None:
        server = Server(host="http://www.ensembl.org")
        dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

    # Map internal attribute keys to display names
    display_name_map = {attr.name: attr.display_name for attr in dataset.attributes.values()}

    all_ids = set().union(*df[id_column])
    attribute_results = batch_query_ensembl(all_ids, attributes, dataset)

    match_flags = []
    match_values = []
    match_attributes = []

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Matching rows"):        
        alias = row[alias_column]
        id_set = row[id_column]

        match_found = False
        matched_value = None
        matched_attr_display = None

        for eid in id_set:
            attr_data = attribute_results.get(eid, {})
            for attr_key, values in attr_data.items():
                if pd.isna(alias):
                    continue  # skip if alias is missing

                alias_lower = str(alias).lower()

                for attr_key, values in attr_data.items():
                    # Filter out NaN and lowercase values for comparison
                    filtered_values = [str(v).lower() for v in values if pd.notna(v)]

                    if alias_lower in filtered_values:
                        match_found = True
                        matched_value = alias  # keep original alias casing
                        matched_attr_display = display_name_map.get(attr_key, attr_key)
                        break
            if match_found:
                break

        match_flags.append(match_found)
        match_values.append(matched_value)
        match_attributes.append(matched_attr_display)

    df = df.copy()
    df["Gene Identifier Match"] = match_flags
    df["matching_value"] = match_values
    df["matching_attribute"] = match_attributes

    return df

def match_alias_to_gene_identifier_by_prefix(
    row: pd.Series,
    boolean_match_col: str,
    alias_symbol_col: str,
    source_of_id_col: str,
    prefix: str,
    suffix_condition=None,
) -> tuple[bool, str | None]:
    """Check and update match status and source information based on alias prefix.

    If the boolean match column is already True, returns existing values.
    If the alias symbol starts with the given prefix, sets match to True and
    sets the source column to a descriptive prefix string.

    :param row: a single row from a DataFrame (Series)
    :param boolean_match_col: name of the column that holds the True/False match status
    :param alias_symbol_col: name of the column containing the alias symbol
    :param source_of_id_col: name of the column where the source label should be added
    :param prefix: string prefix to match against the alias symbol
    :param suffix_condition: function to apply to the suffix after prefix (optional)
    :return: tuple (updated match boolean, updated source string)
    """
    if row[boolean_match_col]:
        return row[boolean_match_col], row[source_of_id_col]

    value = row.get(alias_symbol_col)

    if isinstance(value, str) and value.upper().startswith(prefix.upper()):
        suffix = value[len(prefix):]

        if suffix_condition is None or suffix_condition(suffix):
            return True, f"{prefix} prefix"

    return False, row[source_of_id_col]

def match_alias_to_hgnc_associated_ids(
    analysis_df: pd.DataFrame,
    associated_id_df: pd.DataFrame,
) -> pd.DataFrame:
    """Match aliases to HGNC-associated identifiers.

    :param analysis_df: DataFrame containing aliases and HGNC IDs.
    :param associated_id_df: DataFrame of HGNC-associated identifiers.
    :return: DataFrame with HGNC identifier match results.
    """
    # Make copies to avoid modifying originals
    analysis_df = analysis_df.copy()
    associated_id_df = associated_id_df.copy()

    # Ensure strings and strip spaces
    analysis_df["HGNC_ID"] = analysis_df["HGNC_ID"].astype(str).str.strip()
    analysis_df["alias_symbol"] = analysis_df["alias_symbol"].astype(str).str.strip()
    associated_id_df["HGNC ID"] = associated_id_df["HGNC ID"].astype(str).str.strip()
    associated_id_df["associated_id"] = associated_id_df["associated_id"].astype(str).str.strip()
    associated_id_df["source"] = associated_id_df["source"].astype(str).str.strip()

    analysis_df["alias_symbol_upper"] = analysis_df["alias_symbol"].str.upper()
    associated_id_df["associated_id_upper"] = associated_id_df["associated_id"].str.upper()

    # Merge on HGNC_ID and alias_symbol matching associated_id
    merged = analysis_df.merge(
        associated_id_df,
        how="left",
        left_on=["HGNC_ID", "alias_symbol_upper"],
        right_on=["HGNC ID", "associated_id_upper"]
    )

    # Gene Identifier Match = True if a match occurred
    merged["Gene Identifier Match"] = merged["associated_id"].notna()

    # Rename to match expected output columns
    merged = merged.rename(columns={"source": "associated_source"})

    # Drop merge helper column
    merged = merged.drop(columns=["HGNC ID", "alias_symbol_upper", "associated_id_upper"], errors="ignore")

    return merged

def match_alias_to_ncbi_associated_ids(
    analysis_df: pd.DataFrame,
    associated_id_df: pd.DataFrame,
    associated_id_col: str,
) -> pd.DataFrame:
    """Match aliases to NCBI-associated identifiers.

    :param analysis_df: DataFrame containing aliases and NCBI IDs.
    :param associated_id_df: DataFrame of NCBI-associated identifiers.
    :param associated_id_col: name of the associated identifier column.
    :return: DataFrame with NCBI identifier match results.
    """
    # Ensure all values are strings and stripped
    analysis_df = analysis_df.copy()
    associated_id_df = associated_id_df.copy()

    analysis_df["NCBI_ID"] = analysis_df["NCBI_ID"].astype(str).str.strip()
    analysis_df["alias_symbol"] = analysis_df["alias_symbol"].astype(str).str.strip()
    associated_id_df["GeneID"] = associated_id_df["GeneID"].astype(str).str.strip()
    associated_id_df[associated_id_col] = associated_id_df[associated_id_col].astype(str).str.strip()

    # Rename for merge
    temp = associated_id_df.rename(columns={associated_id_col: f"{associated_id_col}_id"})

    # Perform the merge
    merged = analysis_df.merge(
        temp[["GeneID", f"{associated_id_col}_id"]],
        left_on=["NCBI_ID", "alias_symbol"],
        right_on=["GeneID", f"{associated_id_col}_id"],
        how="left"
    )

    # Determine match
    merged[f"{associated_id_col} Match"] = merged["GeneID"].notna()

    return merged.drop(columns=["GeneID"], errors="ignore")

# Withdrawn Ortholog Symbol Capture
def match_alias_to_mgi_withdrawn_symbols(
    analysis_df: pd.DataFrame,
    mgi_df: pd.DataFrame,
) -> pd.DataFrame:
    """Match aliases to withdrawn MGI symbols.

    :param analysis_df: DataFrame containing alias and primary gene symbols.
    :param mgi_df: DataFrame containing MGI marker symbols, names, and statuses.
    :return: DataFrame with withdrawn MGI symbol match information.
    """
    # Copy inputs to avoid modifying original data
    analysis_df = analysis_df.copy()
    mgi_df = mgi_df.copy()

    # Standardize and strip strings
    analysis_df["alias_symbol"] = (
        analysis_df["alias_symbol"].astype(str).str.strip().str.upper()
    )
    analysis_df["primary_gene_symbol"] = (
        analysis_df["primary_gene_symbol"].astype(str).str.strip().str.upper()
    )

    mgi_df["Marker Symbol"] = (
        mgi_df["Marker Symbol"].astype(str).str.strip().str.upper()
    )
    mgi_df["Marker Name"] = mgi_df["Marker Name"].astype(str).str.strip()
    mgi_df["Status"] = mgi_df["Status"].astype(str).str.strip()

    withdrawn_mgi_df = mgi_df[mgi_df["Status"] == "W"]

    temp = withdrawn_mgi_df.rename(
        columns={"Marker Symbol": "MGI Withdrawn Symbol"}
    )

    merged = analysis_df.merge(
        temp[["MGI Withdrawn Symbol", "Marker Name"]],
        left_on="alias_symbol",
        right_on="MGI Withdrawn Symbol",
        how="left",
    )

    primary_in_marker_name = merged.apply(
        lambda row: (
            pd.notna(row["Marker Name"])
            and row["primary_gene_symbol"] in row["Marker Name"].upper()
        ),
        axis=1,
    )

    merged["MGI Withdrawn Symbol Match"] = (
        merged["MGI Withdrawn Symbol"].notna()
        & primary_in_marker_name
    )

    return merged


# Alternate Abbreviation Symbol Capture
def parse_hgnc(x: str) -> list[str]:
    """Parse a string representation of an HGNC value list.

    :param x: string representation of a list of HGNC values.
    :return: parsed list of HGNC values, or an empty list if the input is empty.
    """
    return list(ast.literal_eval(x)) if x else []