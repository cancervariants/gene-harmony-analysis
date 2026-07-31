"""Functions for analyzing gene alias mentions in PubTator 3 annotations.

This module provides utilities to retrieve PubTator documents, identify gene
alias mentions, evaluate normalization to HGNC, Ensembl, and NCBI Gene
identifiers, and summarize identifier usage for ambiguous gene symbols in the
biomedical literature.
"""

import hashlib
import json
import random
import re
import time
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from http.client import RemoteDisconnected
from pathlib import Path
from typing import Any

import requests
from requests import Response
from requests.exceptions import (
    ChunkedEncodingError,
    ReadTimeout,
    Timeout,
)
from requests.exceptions import (
    ConnectionError as RequestsConnectionError,
)

type JSONValue = (
    str
    | int
    | float
    | bool
    | None
    | list["JSONValue"]
    | dict[str, "JSONValue"]
)

BASE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

class CheckpointMismatchError(ValueError):
    """Raised when a checkpoint does not match the requested alias."""

    def __init__(
        self,
        checkpoint_query: str,
        alias: str,
    ) -> None:
        """Initialize the exception with the mismatched alias values.

        :param checkpoint_query: the alias stored in the checkpoint
        :param alias: the requested alias
        """
        super().__init__(
            f"Checkpoint query {checkpoint_query!r} "
            f"does not match {alias!r}."
        )
class RequestRetriesExceededError(RuntimeError):
    """Raised when a request fails after all retry attempts."""

    def __init__(
        self,
        max_retries: int,
        url: str,
        params: dict[str, Any],
    ) -> None:
        """Initialize the exception with request details.

        :param max_retries: the maximum number of request attempts
        :param url: the requested URL
        :param params: the query parameters used in the request
        """
        super().__init__(
            f"Request failed after {max_retries} attempts. "
            f"URL: {url}; params: {params}"
        )

@dataclass(frozen=True)
class GenePair:
    """Store an alias gene symbol and the identifiers for its associated approved gene.

    :param alias: the alias gene symbol being searched
    :param approved_symbol: the approved human gene symbol associated with the alias
    :param ncbi_gene_id: the NCBI Gene identifier for the approved gene
    :param hgnc_id: the HGNC identifier for the approved gene
    :param ensembl_gene_id: the Ensembl gene identifier for the approved gene
    """

    alias: str
    approved_symbol: str
    ncbi_gene_id: str
    hgnc_id: str
    ensembl_gene_id: str

    @property
    def alias_prefix(self) -> str:
        """Create a lowercase alias value for use in file and directory names.

        return: the lowercase alias gene symbol
        """
        return self.alias.lower()

    @property
    def pair_prefix(self) -> str:
        """Create a lowercase prefix for the alias-approved symbol pair.

        return: a string containing the alias and approved gene symbol
        """
        return (
            f"{self.alias.lower()}_"
            f"{self.approved_symbol.lower()}"
        )

    @property
    def alias_output_dir(self) -> Path:
        """Create and return the shared output directory for an alias.

        return: the path to the alias-level output directory
        """
        path = OUTPUT_DIR / self.alias_prefix
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def result_output_dir(self) -> Path:
        """Create and return the gene-specific result directory for an alias-gene pair.

        return: the path to the result directory for the approved gene
        """
        path = (
            self.alias_output_dir
            / "results"
            / self.approved_symbol.lower()
        )
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def checkpoint_file(self) -> Path:
        """Return the file used to save search progress for an alias.

        return: the path to the alias search checkpoint file
        """
        # Shared by every gene associated with this alias.
        return self.alias_output_dir / "search_checkpoint.json"

    @property
    def pmid_cache_file(self) -> Path:
        """Return the file used to cache PMIDs for an alias.

        return: the path to the alias PMID cache file
        """
        # Shared by every gene associated with this alias.
        return self.alias_output_dir / "pmids.json"

    @property
    def document_cache_dir(self) -> Path:
        """Create and return the directory used to cache PubTator documents.

        return: the path to the alias document cache directory
        """
        # Shared by every gene associated with this alias.
        path = self.alias_output_dir / "document_cache"
        path.mkdir(parents=True, exist_ok=True)
        return path

def get_batch_cache_file(
    pair: GenePair,
    batch: list[str],
) -> Path:
    """Create a stable cache filename based on the PMIDs in a batch.

    :param pair: the alias-gene pair whose document cache directory will be used
    :param batch: the PMIDs included in the document request batch
    return: the path to the JSON cache file for the batch
    """
    batch_key = ",".join(batch)

    batch_hash = hashlib.sha256(
        batch_key.encode("utf-8")
    ).hexdigest()[:16]

    return (
        pair.document_cache_dir
        / f"batch_{batch_hash}.json"
    )


def make_session() -> requests.Session:
    """Create a requests session configured for PubTator API calls.

    return: a requests session with the required request headers
    """
    new_session = requests.Session()
    new_session.headers.update({
        "User-Agent": "gene-identifier-analysis/1.0",
        "Accept": "application/json",
        "Connection": "close",
    })
    return new_session


session = make_session()


def make_identifier_patterns(
    pair: GenePair,
) -> dict[str, re.Pattern]:
    """Create regular expression patterns for accepted gene identifiers.

    :param pair: the alias-gene pair containing the identifiers to match
    return: a dictionary that maps each identifier namespace to its compiled pattern
    """
    hgnc_number = pair.hgnc_id.split(":")[-1]

    return {
        "NCBI Gene": re.compile(
            rf"\b(?:NCBI\s+Gene|Entrez\s+Gene|Gene\s+ID)"
            rf"\s*[:#]?\s*{re.escape(pair.ncbi_gene_id)}\b",
            re.IGNORECASE,
        ),
        "HGNC": re.compile(
            rf"\bHGNC\s*[:#]?\s*{re.escape(hgnc_number)}\b",
            re.IGNORECASE,
        ),
        "Ensembl": re.compile(
            rf"\b{re.escape(pair.ensembl_gene_id)}(?:\.\d+)?\b",
            re.IGNORECASE,
        ),
    }

def extract_pmids(obj: JSONValue) -> set[str]:
    """Recursively extract valid PubMed identifiers from a nested object.

    :param obj: a dictionary, list, or value returned by the PubTator API
    return: the unique PMIDs found in the object
    """
    pmids = set()

    if isinstance(obj, dict):
        for key, value in obj.items():
            if key.lower() in {"pmid", "pmids"}:
                if isinstance(value, (str, int)):
                    if str(value).isdigit():
                        pmids.add(str(value))

                elif isinstance(value, list):
                    pmids.update(
                        str(item)
                        for item in value
                        if str(item).isdigit()
                    )

            pmids.update(extract_pmids(value))

    elif isinstance(obj, list):
        for item in obj:
            pmids.update(extract_pmids(item))

    return pmids


def get_with_retry(
    url: str,
    *,
    params: dict[str, Any],
    timeout: int = 90,
    max_retries: int = 10,
) -> Response:
    """Send a GET request and retry temporary connection, rate-limit, and server failures.

    :param url: the URL to request
    :param params: the query parameters to include in the request
    :param timeout: the maximum number of seconds to wait for each request
    :param max_retries: the maximum number of request attempts
    return: the successful HTTP response
    """
    global session

    for attempt in range(max_retries):
        try:
            response = session.get(
                url,
                params=params,
                timeout=timeout,
            )

            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")

                if retry_after:
                    wait_seconds = float(retry_after)
                else:
                    wait_seconds = min(
                        120,
                        5 * (2**attempt),
                    )

                wait_seconds += random.uniform(0, 2)  # noqa: S311

                print(  # noqa: T201
                    "Rate limited. Waiting "
                    f"{wait_seconds:.1f} seconds..."
                )

                time.sleep(wait_seconds)
                continue

            if response.status_code in {500, 502, 503, 504}:
                wait_seconds = (
                    min(120, 3 * (2**attempt))
                    + random.uniform(0, 2)  # noqa: S311
                )

                print(  # noqa: T201
                    f"Server returned {response.status_code}. "
                    f"Waiting {wait_seconds:.1f} seconds..."
                )

                time.sleep(wait_seconds)
                continue

            response.raise_for_status()

        except (
            RequestsConnectionError,
            RemoteDisconnected,
            ReadTimeout,
            Timeout,
            ChunkedEncodingError,
        ) as error:
            wait_seconds = (
                min(120, 3 * (2**attempt))
                + random.uniform(0, 2)  # noqa: S311
            )

            print(  # noqa: T201
                f"Connection failed: {type(error).__name__}. "
                f"Waiting {wait_seconds:.1f} seconds..."
            )

            session.close()
            session = make_session()

            time.sleep(wait_seconds)

        else:
            return response

    raise RequestRetriesExceededError(
        max_retries=max_retries,
        url=url,
        params=params,
    )


def search_all_pmids_checkpointed(
    pair: GenePair,
    delay: float = 2.0,
) -> set[str]:
    """Search PubTator for every PMID associated with an alias and save progress.

    :param pair: the alias-gene pair whose alias will be searched
    :param delay: the number of seconds to wait between search result pages
    return: the unique PMIDs returned for the alias
    """
    checkpoint_path = pair.checkpoint_file

    if checkpoint_path.exists():
        checkpoint = json.loads(
            checkpoint_path.read_text()
        )

        if checkpoint.get("query") != pair.alias:
            raise CheckpointMismatchError(
                checkpoint.get("query", ""),
                pair.alias,
            )

        all_pmids = set(checkpoint["pmids"])
        page = checkpoint["next_page"]

        print(  # noqa: T201
            f"Resuming {pair.alias} at page {page} with "
            f"{len(all_pmids):,} cached PMIDs."
        )

    else:
        all_pmids = set()
        page = 1

    while True:
        response = get_with_retry(
            f"{BASE_URL}/search/",
            params={
                "text": pair.alias,
                "page": page,
            },
        )

        page_pmids = extract_pmids(response.json())

        if not page_pmids:
            break

        new_pmids = page_pmids - all_pmids

        if not new_pmids:
            break

        all_pmids.update(new_pmids)

        checkpoint_path.write_text(
            json.dumps(
                {
                    "query": pair.alias,
                    "approved_symbol": pair.approved_symbol,
                    "next_page": page + 1,
                    "pmids": sorted(all_pmids),
                },
                indent=2,
            )
        )

        print(  # noqa: T201
            f"{pair.alias} page {page}: "
            f"{len(new_pmids):,} new PMIDs "
            f"({len(all_pmids):,} total)"
        )

        page += 1
        time.sleep(delay)

    pair.pmid_cache_file.write_text(
        json.dumps(
            sorted(all_pmids),
            indent=2,
        )
    )

    print(  # noqa: T201
        f"Search complete for {pair.alias}: "
        f"{len(all_pmids):,} unique PMIDs"
    )

    return all_pmids


def load_or_search_pmids(
    pair: GenePair,
) -> set[str]:
    """Load cached PMIDs for an alias or search PubTator when no cache exists.

    :param pair: the alias-gene pair whose alias PMIDs are needed
    return: the cached or newly retrieved PMIDs for the alias
    """
    if pair.pmid_cache_file.exists():
        set(
            json.loads(
                pair.pmid_cache_file.read_text()
            )
        )

    return search_all_pmids_checkpointed(pair)


def chunked(
    values: Iterable[str],
    size: int = 100,
) -> Iterator[list[str]]:
    """Divide an iterable of strings into lists of a specified size.

    :param values: the string values to divide into batches
    :param size: the maximum number of values in each batch
    return: an iterator that yields lists of values
    """
    values = list(values)

    for start in range(0, len(values), size):
        yield values[start:start + size]


def fetch_documents(
    pair: GenePair,
    pmids: set[str],
    *,
    batch_size: int = 50,
    force_refresh: bool = False,
) -> Iterator[dict[str, Any]]:
    """Load PubTator documents from the local cache or download missing batches.

    :param pair: the alias-gene pair whose shared document cache will be used
    :param pmids: the PMIDs of the documents to retrieve
    :param batch_size: the maximum number of PMIDs requested in each batch
    :param force_refresh: whether to ignore cached batches and download them again
    return: an iterator that yields PubTator document dictionaries
    """
    pair.document_cache_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    sorted_pmids = sorted(pmids)
    total_pmids = len(sorted_pmids)  # noqa: F841
    processed_pmids = 0

    for batch_number, batch in enumerate(  # noqa: B007
        chunked(sorted_pmids, size=batch_size),
        start=1,
    ):
        cache_file = get_batch_cache_file(
            pair,
            batch,
        )

        documents = None

        if cache_file.exists() and not force_refresh:
            try:
                cached_result = json.loads(
                    cache_file.read_text()
                )

                cached_pmids = cached_result.get(
                    "requested_pmids",
                    [],
                )

                # Make sure the file contains this exact batch.
                if cached_pmids == batch:
                    documents = cached_result.get(
                        "documents",
                        [],
                    )

            except (
                json.JSONDecodeError,
                OSError,
            ) as error:
                print(  # noqa: T201
                    f"Could not read {cache_file}: "
                    f"{error}. Downloading again."
                )

        if documents is None:
            response = get_with_retry(
                f"{BASE_URL}/publications/export/biocjson",
                params={
                    "pmids": ",".join(batch),
                    "full": "true",
                },
                timeout=180,
            )

            result = response.json()

            if isinstance(result, list):
                documents = result

            elif (
                isinstance(result, dict)
                and "PubTator3" in result
            ):
                documents = result["PubTator3"]

            elif (
                isinstance(result, dict)
                and "documents" in result
            ):
                documents = result["documents"]

            elif (
                isinstance(result, dict)
                and "id" in result
            ):
                documents = [result]

            else:
                print("Unexpected response structure:")  # noqa: T201
                print(str(result)[:1000])  # noqa: T201
                documents = []

            cache_file.write_text(
                json.dumps(
                    {
                        "alias": pair.alias,
                        "requested_pmids": batch,
                        "documents": documents,
                    },
                    indent=2,
                )
            )

            time.sleep(2)

        processed_pmids += len(batch)

        yield from documents


def analyze_gene_pair(
    pair: GenePair,
    candidate_pmids: set[str],
) -> dict[str, Any]:
    """Analyze papers for an exact alias annotation and accepted identifiers for a gene.

    :param pair: the alias-gene pair and identifiers being analyzed
    :param candidate_pmids: the PMIDs of papers that may contain the alias
    return: a dictionary containing paper sets, counts, and the identifier percentage
    """
    identifier_patterns = make_identifier_patterns(pair)

    papers_with_alias_gene_annotation: set[str] = set()
    papers_with_identifier_in_text: set[str] = set()

    papers_by_namespace = {
        namespace: set()
        for namespace in identifier_patterns
    }

    papers_by_section: dict[str, set[str]] = {}
    annotation_identifier_counts = Counter()
    processed_documents = 0

    for document in fetch_documents(
        pair,
        candidate_pmids,
        # force_refresh=True,
    ):
        processed_documents += 1

        pmid = str(
            document.get("pmid")
            or document.get("id")
            or ""
        ).strip()

        if not pmid:
            continue

        passages = document.get("passages", [])
        has_exact_alias_gene_annotation = False

        # Check whether PubTator tagged the exact alias as a gene.
        for passage in passages:
            for annotation in passage.get("annotations", []):
                infons = annotation.get("infons", {})

                entity_type = (
                    str(infons.get("type", ""))
                    .strip()
                    .casefold()
                )

                mention = (
                    str(annotation.get("text", ""))
                    .strip()
                )

                if (
                    entity_type == "gene"
                    and mention.casefold()
                    == pair.alias.casefold()
                ):
                    has_exact_alias_gene_annotation = True

                    annotation_identifier = infons.get(
                        "identifier"
                    )

                    if annotation_identifier:
                        annotation_identifier_counts[
                            str(annotation_identifier)
                        ] += 1

        # Do not search identifiers in papers without
        # an exact alias gene annotation.
        if not has_exact_alias_gene_annotation:
            continue

        papers_with_alias_gene_annotation.add(pmid)

        # Search each passage in this qualifying document.
        for passage in passages:
            passage_type = (
                passage.get("infons", {})
                .get("type", "unknown")
            )

            passage_text = str(
                passage.get("text", "")
            )

            for namespace, pattern in (
                identifier_patterns.items()
            ):
                if pattern.search(passage_text):
                    papers_by_namespace[namespace].add(pmid)
                    papers_with_identifier_in_text.add(pmid)

                    papers_by_section.setdefault(
                        passage_type,
                        set(),
                    ).add(pmid)

    denominator = len(
        papers_with_alias_gene_annotation
    )

    numerator = len(
        papers_with_identifier_in_text
    )

    percentage = (
        100 * numerator / denominator
        if denominator
        else 0
    )

    return {
        "pair": pair,
        "processed_documents": processed_documents,
        "candidate_papers": len(candidate_pmids),
        "papers_with_alias_gene_annotation":
            papers_with_alias_gene_annotation,
        "papers_with_identifier_in_text":
            papers_with_identifier_in_text,
        "papers_by_namespace": papers_by_namespace,
        "annotation_identifier_counts":
            annotation_identifier_counts,
        "papers_by_section":
            papers_by_section,
        "denominator": denominator,
        "numerator": numerator,
        "percentage": percentage,
    }

def print_analysis_results(results: dict[str, Any]) -> None:
    """Print a readable summary of the gene-pair analysis results.

    :param results: the dictionary returned by analyze_gene_pair
    return: None
    """
    pair: GenePair = results["pair"]

    print(f"Candidate papers for {pair.alias}: {results['candidate_papers']:,}")  # noqa: T201

    print(f"Documents retrieved: {results['processed_documents']:,}")  # noqa: T201

    print(  # noqa: T201
        f"Papers where PubTator tagged the exact alias "
        f"{pair.alias} as a gene: {results['denominator']:,}"
    )

    print(  # noqa: T201
        f"Papers containing an accepted NCBI Gene, HGNC, or "
        f"Ensembl identifier for {pair.approved_symbol}: "
        f"{results['numerator']:,}"
    )

    print(  # noqa: T201
        f"Percentage containing an identifier: "
        f"{results['percentage']:.2f}%"
    )

    print("\nCounts by identifier namespace:")  # noqa: T201

    for namespace in ("NCBI Gene", "HGNC", "Ensembl"):
        print(  # noqa: T201
            f"  {namespace}: "
            f"{len(results['papers_by_namespace'][namespace]):,} papers"
        )

    print("\nIdentifier locations:")  # noqa: T201

    for section in sorted(results["papers_by_section"]):
        pmids = sorted(
            results["papers_by_section"][section],
            key=int,
        )

        print(  # noqa: T201
            f"  {section}: "
            f"{', '.join(pmids)}"
        )
