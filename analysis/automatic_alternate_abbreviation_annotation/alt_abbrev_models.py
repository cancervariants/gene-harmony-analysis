import re
from collections.abc import Mapping
from enum import StrEnum
from typing import Any
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from pydantic import BaseModel, ConfigDict
from wags_llm.prompts import BasePromptTemplate


class SkipReason(StrEnum):
    """Reason for why LLM invocation was skipped"""

    HSA_PREFIX = "hsa_prefix"
    "alias symbols with 'HSA-' prefix are gene identifiers in miRBase"
    "Example: HSA-MIR-21 is the miRBase identifier for MIR21."
    EXTRA_CHARACTERS = "extra_characters"
    "an alias symbol cannot be an alternate abbreviation if it has extra characters compared to the primary gene symbol or gene name"
    LOW_LCS_SIMILARITY = "low_lcs_similarity"
    "the lower the LCS similarity score between the alias symbol and the primary gene symbol/gene name, the less likely the alias symbol"
    "is an alternate abbreviation. The threshold is assigned based on the distribution of LCS similarity scores in the manually annotated dataset."


class AlternateAbbreviationPredictionResult(BaseModel):
    """Model for LLM and human curator result for determining whether an
    alias symbol corresponds to the primary gene symbol.
    """

    model_config = ConfigDict(extra="forbid", use_enum_values=True)

    llm_annotation: bool | None = None
    llm_skip_reason: SkipReason | None = None
    error_message: str | None = None
    lcs_similarity_score: float | None = None

@dataclass
class RunResult:
    """Store performance metrics and summaries for a single evaluation run."""

    temperature: float
    run_idx: int

    llm_accuracy: float
    llm_coverage: float | None
    llm_summary: any
    llm_metrics: pl.DataFrame

    system_accuracy: float | None = None
    system_summary: any = None
    system_metrics: pl.DataFrame | None = None


@dataclass
class AlternateAbbreviationPrompt(BasePromptTemplate):
    """Version 1 prompt for predicting alternate abbreviation relationships."""

    name = "alternate_abbreviation_annotation"
    PROMPT_DIR = Path(__file__).parent / "alt_abbrev_prompts"

    def __init__(self, version: str):
        self.version = version

    @property
    def prompt_path(self) -> Path:
        return self.PROMPT_DIR / f"{self.version}.txt"

    def build_system_prompt(self) -> str:
        if not self.prompt_path.exists():
            available_versions = sorted(
                p.stem for p in self.PROMPT_DIR.glob("*.txt")
            )

            raise FileNotFoundError(
                f"Prompt version '{self.version}' not found. "
                f"Available versions: {available_versions}"
            )

        return self.prompt_path.read_text()

    def build_user_prompt(
        self,
        payload: Mapping[str, Any],
    ) -> str:
        """Build the user prompt for a single alias symbol.

        :param payload: The alias symbol, HGNC ID, primary gene symbol, official gene name to be evaluated, 
        :returns: User prompt text
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
        return {
            "gene_symbol": gene_symbol,
            "primary_gene_symbol": primary_gene_symbol,
            "gene_name": gene_name,
            "hgnc_id": hgnc_id,
        }