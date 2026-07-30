import hashlib
import re
import json
import random
import time

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from http.client import RemoteDisconnected
from pathlib import Path
from typing import Any

import requests
from requests import Response
from requests.exceptions import (
    ChunkedEncodingError,
    ConnectionError,
    ReadTimeout,
    Timeout,
)
BASE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
@dataclass(frozen=True)
class GenePair:
    alias: str
    approved_symbol: str
    ncbi_gene_id: str
    hgnc_id: str
    ensembl_gene_id: str

    @property
    def alias_prefix(self) -> str:
        return self.alias.lower()

    @property
    def pair_prefix(self) -> str:
        return (
            f"{self.alias.lower()}_"
            f"{self.approved_symbol.lower()}"
        )

    @property
    def alias_output_dir(self) -> Path:
        """Shared search and document data for this alias."""
        path = OUTPUT_DIR / self.alias_prefix
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def result_output_dir(self) -> Path:
        """Gene-specific output for this alias-gene pairing."""
        path = (
            self.alias_output_dir
            / "results"
            / self.approved_symbol.lower()
        )
        path.mkdir(parents=True, exist_ok=True)
        return path

    @property
    def checkpoint_file(self) -> Path:
        # Shared by every gene associated with this alias.
        return self.alias_output_dir / "search_checkpoint.json"

    @property
    def pmid_cache_file(self) -> Path:
        # Shared by every gene associated with this alias.
        return self.alias_output_dir / "pmids.json"

    @property
    def document_cache_dir(self) -> Path:
        # Shared by every gene associated with this alias.
        path = self.alias_output_dir / "document_cache"
        path.mkdir(parents=True, exist_ok=True)
        return path
    
def get_batch_cache_file(
    pair: GenePair,
    batch: list[str],
) -> Path:
    """
    Create a stable filename based on the PMIDs in the batch.
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


def extract_pmids(obj: Any) -> set[str]:
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
    params: dict,
    timeout: int = 90,
    max_retries: int = 10,
) -> Response:
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
                        5 * (2 ** attempt),
                    )

                wait_seconds += random.uniform(0, 2)

                print(
                    f"Rate limited. Waiting "
                    f"{wait_seconds:.1f} seconds..."
                )

                time.sleep(wait_seconds)
                continue

            if response.status_code in {500, 502, 503, 504}:
                wait_seconds = (
                    min(120, 3 * (2 ** attempt))
                    + random.uniform(0, 2)
                )

                print(
                    f"Server returned {response.status_code}. "
                    f"Waiting {wait_seconds:.1f} seconds..."
                )

                time.sleep(wait_seconds)
                continue

            response.raise_for_status()
            return response

        except (
            ConnectionError,
            RemoteDisconnected,
            ReadTimeout,
            Timeout,
            ChunkedEncodingError,
        ) as error:
            wait_seconds = (
                min(120, 3 * (2 ** attempt))
                + random.uniform(0, 2)
            )

            print(
                f"Connection failed: {type(error).__name__}. "
                f"Waiting {wait_seconds:.1f} seconds..."
            )

            session.close()
            session = make_session()

            time.sleep(wait_seconds)

    raise RuntimeError(
        f"Request failed after {max_retries} attempts. "
        f"URL: {url}; params: {params}"
    )


def search_all_pmids_checkpointed(
    pair: GenePair,
    delay: float = 2.0,
) -> set[str]:
    checkpoint_path = pair.checkpoint_file

    if checkpoint_path.exists():
        checkpoint = json.loads(
            checkpoint_path.read_text()
        )

        if checkpoint.get("query") != pair.alias:
            raise ValueError(
                f"Checkpoint query "
                f"{checkpoint.get('query')!r} does not match "
                f"{pair.alias!r}."
            )

        all_pmids = set(checkpoint["pmids"])
        page = checkpoint["next_page"]

        print(
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

        print(
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

    print(
        f"Search complete for {pair.alias}: "
        f"{len(all_pmids):,} unique PMIDs"
    )

    return all_pmids


def load_or_search_pmids(
    pair: GenePair,
) -> set[str]:
    if pair.pmid_cache_file.exists():
        pmids = set(
            json.loads(
                pair.pmid_cache_file.read_text()
            )
        )

        return pmids

    return search_all_pmids_checkpointed(pair)


def chunked(
    values: Iterable[str],
    size: int = 100,
):
    values = list(values)

    for start in range(0, len(values), size):
        yield values[start:start + size]


def fetch_documents(
    pair: GenePair,
    pmids: set[str],
    *,
    batch_size: int = 50,
    force_refresh: bool = False,
):
    """
    Load PubTator documents from the local cache when available.

    If a batch is not cached, download it from PubTator and save it.
    """
    pair.document_cache_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    sorted_pmids = sorted(pmids)
    total_pmids = len(sorted_pmids)
    processed_pmids = 0

    for batch_number, batch in enumerate(
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
                print(
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
                print("Unexpected response structure:")
                print(str(result)[:1000])
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
    identifier_patterns = make_identifier_patterns(pair)

    papers_with_alias_gene_annotation: set[str] = set()
    papers_with_identifier_in_text: set[str] = set()

    papers_by_namespace = {
        namespace: set()
        for namespace in identifier_patterns
    }

    identifier_section_counts = Counter()
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
                    identifier_section_counts[
                        passage_type
                    ] += 1

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
        "identifier_section_counts":
            identifier_section_counts,
        "denominator": denominator,
        "numerator": numerator,
        "percentage": percentage,
    }

def print_analysis_results(results: dict[str, Any]) -> None:
    pair: GenePair = results["pair"]

    print(f"Candidate papers for {pair.alias}: {results['candidate_papers']:,}")

    print(f"Documents retrieved: {results['processed_documents']:,}")

    print(
        f"Papers where PubTator tagged the exact alias "
        f"{pair.alias} as a gene: {results['denominator']:,}"
    )

    print(
        f"Papers containing an accepted NCBI Gene, HGNC, or "
        f"Ensembl identifier for {pair.approved_symbol}: "
        f"{results['numerator']:,}"
    )

    print(
        f"Percentage containing an identifier: "
        f"{results['percentage']:.2f}%"
    )

    print("\nCounts by identifier namespace:")

    for namespace in ("NCBI Gene", "HGNC", "Ensembl"):
        print(
            f"  {namespace}: "
            f"{len(results['papers_by_namespace'][namespace]):,} papers"
        )

    print("\nIdentifier locations:")

    for section, count in results["identifier_section_counts"].most_common():
        print(f"  {section}: {count:,}")