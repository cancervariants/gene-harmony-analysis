"""Define data models and prompt templates for alternate abbreviation annotation."""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

import polars as pl
from pydantic import BaseModel, ConfigDict
from wags_llm.prompts import BasePromptTemplate


class MatchedRule(StrEnum):
    """Rules that determine if an alias gene symbol is NOT an alternate abbreviation"""

    CONFLICT_CATEGORY = "conflict_category"
    # alias gene symbols that have been found to be gene identifiers or clone symbols are not alternate abbreviations
    # Example: HSA-MIR-21 is the miRBase identifier for MIR21.

    EXTRA_CHARACTERS = "extra_characters"
    # an alias symbol cannot be an alternate abbreviation if it has extra characters compared to the primary gene symbol or gene name

    LOW_LCS_SIMILARITY = "low_lcs_similarity"
    # the lower the LCS similarity score between the alias symbol and the primary gene symbol/gene name, the less likely the alias symbol
    # is an alternate abbreviation. The threshold is assigned based on the distribution of LCS similarity scores in the manually annotated dataset.

class RuleResult(BaseModel):
    """Result of rule-based evaluation for alias gene symbol resolution.

    Stores outputs from heuristic filters used to determine whether an alias
    gene symbol should be excluded from LLM-based evaluation.

    :param lcs_similarity_score: Normalized longest common subsequence similarity
        between alias and primary gene symbol (if computed)
    :param num_extra_characters: Number of characters in alias not present in
        the gene name (if computed)
    :param conflicting_category: List of captured annotation categories that
        indicate rule-based conflicts (if any)
    :param matched_rule: Rule that triggered exclusion, if any
    """

    lcs_similarity_score: float | None = None
    num_extra_characters: int | None = None
    conflicting_category: list[str] | None = None
    matched_rule: MatchedRule | None = None

class AlternateAbbreviationPredictionResult(BaseModel):
    """Model for LLM and human curator result for determining whether an
    alias symbol corresponds to the primary gene symbol.
    """

    model_config = ConfigDict(extra="forbid", use_enum_values=True)

    llm_annotation: bool | None = None
    matched_rule: MatchedRule | None = None
    error_message: str | None = None
    lcs_similarity_score: float | None = None
    conflicting_category: list[str] | None = None
    num_extra_characters: float | None = None

@dataclass
class RunResult:
    """Store performance metrics and summaries for a single evaluation run."""

    temperature: float
    run_idx: int

    llm_accuracy: float
    llm_coverage: float | None
    llm_summary: Any
    llm_metrics: pl.DataFrame

    system_accuracy: float | None = None
    system_summary: Any = None
    system_metrics: pl.DataFrame | None = None


@dataclass
class AlternateAbbreviationPrompt(BasePromptTemplate):
    """Version 1 prompt for predicting alternate abbreviation relationships."""

    name = "alternate_abbreviation_annotation"
    PROMPT_DIR = Path(__file__).parent / "alt_abbrev_prompts"

    def __init__(self, version: str):
        """Initialize prompt configuration with a specific prompt version.

        :param version: Prompt version identifier used to locate the corresponding prompt file
        """
        self.version = version

    @property
    def prompt_path(self) -> Path:
        """Path to the system prompt file for the configured version.

        Constructs the file path by combining the prompt directory with the
        selected prompt version filename.

        :return: Path object pointing to the prompt .txt file
        """
        return self.PROMPT_DIR / f"{self.version}.txt"

    def build_system_prompt(self) -> str:
        """Load the system prompt text for the configured prompt version.

        If the prompt file for the selected version does not exist, raises an error
        listing all available prompt versions in the prompt directory.

        :return: System prompt text loaded from the corresponding prompt file
        :raises FileNotFoundError: If the configured prompt version is not found
        """
        if not self.prompt_path.exists():
            available_versions = sorted(
                p.stem for p in self.PROMPT_DIR.glob("*.txt")
            )

            error_msg = (
                f"Prompt version '{self.version}' not found. "
                f"Available versions: {available_versions}"
            )

            raise FileNotFoundError(error_msg)

        return self.prompt_path.read_text()

    def build_user_prompt(
        self,
        payload: Mapping[str, Any],
    ) -> str:
        """Build the user prompt for a single alias symbol.

        :param payload: The alias symbol, HGNC ID, primary gene symbol, official gene name to be evaluated,
        :return: User prompt text
        """
        return (
            f"Alias Symbol: {payload['gene_symbol']}\n"
            f"Primary Gene Symbol: {payload['primary_gene_symbol']}\n"
            f"Official Gene Name: {payload['gene_name']}\n"
            f"HGNC ID: {payload['hgnc_id']}\n"
        )

    def build_payload(
        self,
        gene_symbol: str,
        primary_gene_symbol: str,
        gene_name: str,
        hgnc_id: str,
    ) -> dict:
        """Build payload for alternate abbreviation annotation request.

        :param gene_symbol: Alias gene symbol being evaluated
        :param primary_gene_symbol: Official HGNC primary gene symbol
        :param gene_name: Full gene name associated with the gene
        :param hgnc_id: HGNC identifier for the gene
        :return: Dictionary payload used for LLM annotation input
        """
        return {
            "gene_symbol": gene_symbol,
            "primary_gene_symbol": primary_gene_symbol,
            "gene_name": gene_name,
            "hgnc_id": hgnc_id,
        }
