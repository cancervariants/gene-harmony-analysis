import re
from collections.abc import Mapping
from enum import StrEnum
from typing import Any
from dataclasses import dataclass

import polars as pl
from pydantic import BaseModel, ConfigDict
from wags_llm.prompts import BasePromptTemplate, build_empty_registry


class SkipReason(StrEnum):
    """Reason for why LLM invocation was skipped"""

    HSA_PREFIX = "hsa_prefix"
    EXTRA_CHARACTERS = "extra_characters"
    LOW_LCS_SIMILARITY = "low_lcs_similarity"


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

PROMPT_NAME = "alias_symbol_annotation:alternate_abbreviation"
PROMPT_VERSION = "v1"


class AlternateAbbreviationPromptV1(BasePromptTemplate):
    """Version 1 prompt for predicting alternate abbreviation relationships."""

    version = PROMPT_VERSION
    name = PROMPT_NAME

    def build_system_prompt(self) -> str:
        """Build the system prompt for predicting whether an alias is an alternate abbreviation of the primary gene symbol.

        :returns: System prompt text.
        """
        return (
            "Role: You are a biomedical gene nomenclature curator trained in HGNC gene identity and alias symbol resolution.\n"
            "You will get fired if your accuracy is less than 95%.\n"
            "Background:\n"
            "An alternate abbreviation is when an alias symbol represents the same official gene name or the\n"
            "primary gene symbol but with different letters. If the alias symbol is representing a different description or a previous name of\n"
            "the gene then it is not an alternate abbreviation. Be careful with alias symbols that seem to be shortened versions\n"
            "of the primary gene symbol, they may be a family name and therefore not an alternate abbreviation. Alias symbols\n"
            "have extra characters may be alternate abbreviations unless they have characters that are not present in the official\n"
            "gene name provided in the prompt. Keep in mind, a gene name with a number after the name is not the same gene as a gene name\n"
            "with no numbers or different numbers after the name.\n"
            "Task:\n"
            "Determine whether the alias symbol is an abbreviation variant of the official HGNC gene name.\n"
        )

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