"""Functions used in the ranked_category_resolver notebook."""

import ast
import importlib
import re
import time
from collections.abc import Iterable, Iterator
from typing import Any

import gene_ids_in_lit.functions as giilfn
import pandas as pd


importlib.reload(giilfn)


def find_aliases_in_document(
    document: dict[str, Any],
    aliases: list[str],
) -> list[str]:
    """Find gene aliases appearing in a PubTator document.

    Alias matching is case-insensitive and requires each alias to appear as a
    standalone token rather than as part of a larger word.

    :param document: PubTator document containing a ``passages`` field.
    :param aliases: Gene aliases to search for in the document text.
    :return: Sorted list of unique aliases found in the document.
    :rtype: list[str]
    """
    found = set()

    for passage in document.get("passages", []):
        text = str(passage.get("text", ""))

        for alias in aliases:
            pattern = re.compile(
                rf"(?<!\w){re.escape(alias)}(?!\w)",
                re.IGNORECASE,
            )

            if pattern.search(text):
                found.add(alias)

    return sorted(found)


def fetch_documents_by_pmids(
    pmids: Iterable[str | int],
    batch_size: int = 50,
) -> Iterator[dict[str, Any]]:
    """Fetch PubTator documents for a collection of PMIDs.

    PMIDs are deduplicated, sorted, and submitted to the PubTator API in
    batches. Several possible PubTator response formats are supported.

    :param pmids: Collection of PubMed identifiers.
    :param batch_size: Maximum number of PMIDs requested per API call.
    :yield: PubTator documents returned by the API.
    :rtype: Iterator[dict[str, Any]]
    """
    sorted_pmids = sorted({str(pmid) for pmid in pmids})

    for batch in giilfn.chunked(sorted_pmids, size=batch_size):
        response = giilfn.get_with_retry(
            f"{giilfn.BASE_URL}/publications/export/biocjson",
            params={
                "pmids": ",".join(batch),
                "full": "true",
            },
            timeout=180,
        )

        result = response.json()

        if isinstance(result, list):
            documents = result
        elif isinstance(result, dict) and "PubTator3" in result:
            documents = result["PubTator3"]
        elif isinstance(result, dict) and "documents" in result:
            documents = result["documents"]
        elif isinstance(result, dict) and "id" in result:
            documents = [result]
        else:
            documents = []

        yield from documents

        time.sleep(2)


def find_row_aliases(
    row: pd.Series,
    documents_by_pmid: dict[str, dict[str, Any]],
) -> list[str]:
    """Find ambiguous aliases associated with a DataFrame row.

    :param row: DataFrame row containing ``PMID`` and ``ambiguous_symbol``.
    :param documents_by_pmid: Mapping of PMID strings to PubTator documents.
    :return: Sorted list of aliases found in the associated document. Returns
        an empty list if no document is available for the PMID.
    :rtype: list[str]
    """
    pmid = str(row["PMID"])
    document = documents_by_pmid.get(pmid)

    if document is None:
        return []

    return find_aliases_in_document(
        document,
        row["ambiguous_symbol"],
    )


def document_has_full_text(
    document: dict[str, Any],
) -> bool:
    """Determine whether a PubTator document contains full-text passages.

    A document is considered to have full text if its passage types include
    content beyond only a title and abstract.

    :param document: PubTator document containing passage metadata.
    :return: True if the document appears to contain full text, otherwise
        False.
    :rtype: bool
    """
    passage_types = {
        passage.get("infons", {}).get("type", "unknown")
        for passage in document.get("passages", [])
    }

    return not passage_types <= {"title", "abstract"}


def analyze_row(
    row: pd.Series,
    documents_by_pmid: dict[str, dict[str, Any]],
) -> pd.Series:
    """Analyze PubTator content associated with a DataFrame row.

    Determine whether full text is available for the row's PMID and search
    the associated PubTator document for ambiguous gene aliases.

    :param row: DataFrame row containing ``PMID`` and ``ambiguous_symbol``.
    :param documents_by_pmid: Mapping of PMID strings to PubTator documents.
    :return: Series containing ``full_text_available`` and ``aliases_found``.
    :rtype: pandas.Series
    """
    pmid = str(row["PMID"])
    document = documents_by_pmid.get(pmid)

    if document is None:
        return pd.Series(
            {
                "full_text_available": False,
                "aliases_found": [],
            }
        )

    return pd.Series(
        {
            "full_text_available": document_has_full_text(document),
            "aliases_found": find_aliases_in_document(
                document,
                row["ambiguous_symbol"],
            ),
        }
    )


def to_set(value: Any) -> set:
    """Convert a string representation of a set to a Python set.

    Existing sets are returned unchanged. Missing values, invalid strings,
    and values that do not evaluate to sets are converted to an empty set.

    :param value: Value to convert.
    :return: Parsed set, or an empty set if conversion is not possible.
    :rtype: set
    """
    if pd.isna(value):
        return set()

    if isinstance(value, set):
        return value

    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return set()

        if isinstance(parsed, set):
            return parsed

    return set()


def clean_id(value: Any) -> str:
    """Normalize an identifier by removing its source prefix.

    For identifiers containing a colon, only the portion after the final
    colon is retained. Leading and trailing whitespace is removed.

    For example, ``HGNC:277`` becomes ``277`` and ``GENE ID:148`` becomes
    ``148``.

    :param value: Identifier to normalize.
    :return: Normalized identifier.
    :rtype: str
    """
    value = str(value).strip()
    return value.split(":")[-1].strip()


def id_in_set(
    value: Any,
    values: set,
) -> bool:
    """Determine whether an identifier is contained in a set of identifiers.

    Normalize both the single identifier and the identifiers in the set using
    ``clean_id`` before comparison.

    :param value: Identifier to search for.
    :param values: Set of identifiers to search.
    :return: True if the normalized identifier occurs in the normalized set,
        otherwise False.
    :rtype: bool
    """
    if pd.isna(value):
        return False

    if not isinstance(values, set):
        return False

    normalized_value = clean_id(value)

    normalized_values = {
        clean_id(item)
        for item in values
        if pd.notna(item)
    }

    return normalized_value in normalized_values


def check_rank_match(
    row: pd.Series,
    capture_df: pd.DataFrame,
    rank_map: dict[str, int],
) -> pd.Series:
    """Compare a gene row against its highest-ranked capture category.

    Match rows in ``capture_df`` to the gene name and rank their capture
    categories using ``rank_map``. Select the highest-ranked capture row and
    compare its HGNC, NCBI, and Ensembl identifiers with those in the input
    row.

    ``rank_match`` is set to ``pd.NA`` when matching capture rows exist but
    no usable capture category is available for ranking.

    :param row: Row from the ambiguous gene DataFrame containing ``name``,
        ``HGNC_ID``, ``NCBI_ID``, and ``ENSG_ID``.
    :param capture_df: DataFrame containing captured gene symbols, categories,
        and identifier sets.
    :param rank_map: Mapping of capture category names to numeric ranks. Lower
        numbers indicate higher priority.
    :return: Series containing ``rank_match``, ``rank_status``,
        ``winning_category``, and ``winning_HGNC_ID``.
    :rtype: pandas.Series
    """
    gene_name = str(row["name"]).casefold()

    matches = capture_df[
        capture_df["gene_symbol"].astype(str).str.casefold() == gene_name
    ].copy()

    if matches.empty:
        return pd.Series(
            {
                "rank_match": False,
                "rank_status": "no gene_symbol match",
                "winning_category": pd.NA,
                "winning_HGNC_ID": pd.NA,
            }
        )

    captured_as = (
        matches["captured as:"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    matches = matches[captured_as.ne("")].copy()

    if matches.empty:
        return pd.Series(
            {
                "rank_match": pd.NA,
                "rank_status": "no captured as",
                "winning_category": pd.NA,
                "winning_HGNC_ID": pd.NA,
            }
        )

    matches["captured_categories"] = (
        matches["captured as:"]
        .astype(str)
        .str.split(",")
        .apply(
            lambda categories: [
                category.strip()
                for category in categories
                if category.strip()
            ]
        )
    )

    def best_rank(categories: list[str]) -> int | float:
        """Return the highest-priority rank represented by the categories.

        :param categories: Capture categories associated with a capture row.
        :return: Lowest numeric rank in ``rank_map``, or positive infinity if
            none of the categories occur in ``rank_map``.
        :rtype: int | float
        """
        ranks = [
            rank_map[category]
            for category in categories
            if category in rank_map
        ]

        return min(ranks) if ranks else float("inf")

    matches["category_rank"] = (
        matches["captured_categories"].apply(best_rank)
    )

    ranked_matches = matches[
        matches["category_rank"] != float("inf")
    ].copy()

    if ranked_matches.empty:
        return pd.Series(
            {
                "rank_match": pd.NA,
                "rank_status": "captured as not in rank_map",
                "winning_category": pd.NA,
                "winning_HGNC_ID": pd.NA,
            }
        )

    # Find the highest-priority rank
    best_rank = ranked_matches["category_rank"].min()

    # Keep all rows tied at that rank
    best_candidates = ranked_matches[
        ranked_matches["category_rank"] == best_rank
    ].copy()

    # Count the number of relationships on each tied candidate row
    best_candidates["relationship_count"] = (
        best_candidates["captured_categories"].apply(len)
    )

    # Pick the tied row with the most relationships
    best_row = best_candidates.loc[
        best_candidates["relationship_count"].idxmax()
    ]

    winning_categories = [
        category
        for category in best_row["captured_categories"]
        if category in rank_map
    ]

    winning_category = min(
        winning_categories,
        key=lambda category: rank_map[category],
    )

    hgnc_match = id_in_set(
        row["normalized_HGNC_ID"],
        best_row["HGNC_ID"],
    )

    ncbi_match = id_in_set(
        row["normalized_NCBI_ID"],
        best_row["NCBI_ID"],
    )

    rank_match = hgnc_match or ncbi_match

    return pd.Series(
        {
            "rank_match": rank_match,
            "rank_status": "matched" if rank_match else "ID mismatch",
            "winning_category": winning_category,
            "winning_HGNC_ID": best_row["HGNC_ID"],
        }
    )