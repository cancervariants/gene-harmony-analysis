import time
import re
from typing import Any
import math

import requests
import polars as pl
try:
    from .alt_abbrev_models import (
        AlternateAbbreviationPredictionResult,
        AlternateAbbreviationPrompt,
        MatchedRule,
        RuleResult,
        RunResult,
    )
except ImportError:
    from alt_abbrev_models import (
        AlternateAbbreviationPredictionResult,
        AlternateAbbreviationPrompt,
        MatchedRule,
        RuleResult,
        RunResult,
    )
from rapidfuzz.distance import LCSseq
from tqdm.notebook import tqdm
from wags_llm.cache import InMemoryCache
from wags_llm.client import BedrockClaudeJsonClient
from wags_llm.prompts import build_empty_registry
from wags_llm.services import StructuredTaskRunner

def get_gene_name(hgnc_id: str) -> str:
    """Retrieve official gene name from HGNC for each sample in the set.

    :param hgnc_id: The HGNC ID to retrieve the official gene name for
    :return: The official HGNC gene name or an error message
    """
    url = f"https://rest.genenames.org/fetch/hgnc_id/{hgnc_id}"
    headers = {"Accept": "application/json"}  # Request JSON format
    time.sleep(0.5)  # Sleep to avoid hitting API rate limits
    response = requests.get(url, headers=headers, timeout=10)

    if response.status_code != 200:
        return f"Error: {response.status_code}"

    data = response.json()

    try:
        return data["response"]["docs"][0]["name"]
    except IndexError:
        return f"No gene found for HGNC ID {hgnc_id}"
    
MODEL_ID = "us.anthropic.claude-sonnet-4-6"
REGION_NAME = "us-east-1"
PROFILE_NAME = "dev-account"
MAX_TOKENS = 150


def build_llm_task_runner(
    model_id: str,
    region_name: str,
    profile_name: str,
    max_tokens: int,
    temperature: float,
) -> StructuredTaskRunner:
    """Build LLM alternate abbreviation alias annotator

    :param model_id: Bedrock model identifier.
    :param region_name: AWS region for the Bedrock runtime client.
    :param profile_name: AWS profile name.
    :param max_tokens: Maximum number of tokens to request from the model.
    :param temperature: Sampling temperature.
    :return: Configured structured task runner instance.
    """
    registry = build_empty_registry()
    registry.register(AlternateAbbreviationPrompt(version="v1"))
    llm_client = BedrockClaudeJsonClient(
        model_id=model_id,
        region_name=region_name,
        profile_name=profile_name,
        max_tokens=max_tokens,
        temperature=temperature,
    )
    cache = InMemoryCache()
    return StructuredTaskRunner(
        client=llm_client, prompt_registry=registry, cache=cache
    )

def rule_based_evaluation(
        df: pl.DataFrame,
        alias_symbol: str,
        primary_gene_symbol: str,
        gene_name: str,
        threshold: float = 0.20,
    ) -> RuleResult:
    """Evaluate the alias gene symbol using rules. If the alias meets the criteria, it will be marked as NOT an alternate abbreviation and will NOT be sent to the LLM for evaluation.

    Applies a series of heuristic filters to the alias symbol
    - checks for captured categories "Gene Identifier Symbol" and "Clone Symbol", not an alternate abbreviation alias
    because it is an identifier being used as an alias
    - 3 or more extra characters in the alias, not in the gene name
    - longest common subsequence (LCS) similarity to the primary gene symbol to be at
    least whatever the threshold is set to

    :param df: DataFrame containing alias annotation metadata
    :param alias_symbol: The alias gene symbol being evaluated
    :param primary_gene_symbol: The official primary gene symbol
    :param gene_name: The full gene name associated with the gene
    :param threshold: Minimum LCS similarity score required to pass filtering
    :return: RuleResult containing:
        - lcs_similarity_score: float similarity score between alias and primary gene symbol
        - num_extra_characters: number of extra characters in alias not present in gene name
        - conflicting_category: list of matched captured categories if rule is triggered, otherwise None
        - matched_rule: MatchedRule if a rule is triggered, otherwise None
    """

    def calc_lcs_similarity(s1: str, s2: str) -> float:
        """Calculate the longest common subsequence (LCS) similarity

        :param s1: the first string, alias symbol
        :param s2: the second string, primary gene symbol
        :return: value between 0 and 1
        """
        return LCSseq.normalized_similarity(s1, s2)

    def count_extra_characters(alias: str, gene_name: str) -> int:
        """Denote if the alias symbol has characters not present in the gene name

        :param alias: gene alias symbol
        :param gene_name: HGNC designated official gene name
        :return: Number of extra characters
        """
        norm_alias = re.sub(r"[^A-Z0-9]", "", alias.upper())

        alias_set = set(norm_alias)
        name_set = set(gene_name)

        extra_chars = alias_set - name_set
        return len(extra_chars)

    def find_conflicting_category(
        df: pl.DataFrame,
        alias: str,
        primary: str,
    ) -> list[str] | None:
        """Return conflicting category if alias row matches rule, else None.

        :param df: input dataframe containing alias annotations
        :param alias: alias gene symbol (normalized)
        :param primary: primary gene symbol (normalized)
        :return: list of conflicting categories if present, otherwise None
        """
        normalized = df.with_columns(
            pl.col("captured_category_list")
            .cast(pl.Utf8)
            .str.replace_all(r"[\[\]']", "")
            .str.split(r",\s*")
            .alias("captured_category_list_norm")
        )

        result = normalized.filter(
            (pl.col("gene_symbol").str.to_uppercase() == alias) &
            (pl.col("primary_gene_symbol").str.to_uppercase() == primary) &
            (
                pl.col("captured_category_list_norm").list.contains("Gene Identifier Symbol") |
                pl.col("captured_category_list_norm").list.contains("Clone Symbol")
            )
        )

        if result.is_empty():
            return None

        return (
            result
            .select("captured_category_list_norm")
            .explode("captured_category_list_norm")
            .unique()
            .get_column("captured_category_list_norm")
            .to_list()
        )

    alias = alias_symbol.upper()
    primary = primary_gene_symbol.upper()
    name = gene_name.upper()

    extra_count = count_extra_characters(alias, name)

    lcs_score = calc_lcs_similarity(alias, primary)

    conflicting_category = find_conflicting_category(df, alias, primary)

    if conflicting_category:
        return RuleResult(
            conflicting_category=conflicting_category,
            matched_rule=MatchedRule.CONFLICT_CATEGORY
        )

    # 3 or more extra characters in the alias that is not present in the gene name indicates that the alias is referring to something other than what the gene name is representing.
    # 1 or 2 extra characters may be acceptable since there are cases where a word is implied but not present in the gene name
    # Ex: the gene TRBV11-3, T cell receptor beta variable 11-3, has an alias (TCRBV11S3) that IS an alternate abbreviation where S represents "segment"
    # The gene is the third segment in the 11th gene family
    # Another example could be that the name refers to is as ABC-A but the alternate abbreviation is ABC-1 where the difference is using an alphabetical character vs a numerical character to indicate the gene as the first
    if extra_count >= 3:
        return RuleResult(
            lcs_similarity_score=None,
            num_extra_characters=extra_count,
            conflicting_category=None,
            matched_rule=MatchedRule.EXTRA_CHARACTERS
        )

    if lcs_score < threshold:
        return RuleResult(
            lcs_similarity_score=lcs_score,
            num_extra_characters=None,
            conflicting_category=None,
            matched_rule=MatchedRule.LOW_LCS_SIMILARITY
        )

    return RuleResult(
        lcs_similarity_score=lcs_score,
        num_extra_characters=extra_count,
        conflicting_category=None,
        matched_rule=None
    )

def get_alt_abbreviation_annotation(
    df: pl.DataFrame,
    task_runner: StructuredTaskRunner,
    prompt: AlternateAbbreviationPrompt,
    gene_symbol: str,
    primary_gene_symbol: str,
    gene_name: str,
    hgnc_id: str,
) -> AlternateAbbreviationPredictionResult:
    """Determine whether a gene symbol is an alternate abbreviation of a
    primary HGNC gene symbol.

    Alias symbols are first evaluated using rule-based criteria. If a rule
    determines that the alias is not an alternate abbreviation, a result
    containing the matched rule and similarity score is returned. Otherwise,
    the configured prompt is executed and the resulting annotation is returned
    along with the computed similarity score.

    :param df: input dataframe containing alias annotations
    :param task_runner: Structured task runner used to execute the LLM prompt.
    :param prompt: Configured prompt instance for alternate abbreviation annotation.
    :param gene_symbol: Candidate alias symbol to evaluate.
    :param primary_gene_symbol: Primary HGNC-approved gene symbol being
        compared against.
    :param gene_name: Full gene name associated with the primary gene symbol.
    :param hgnc_id: HGNC identifier for the primary gene.
    :return: Annotation result containing the rule-based decision or LLM
        prediction, similarity score, and any error information.
    """
    rule_result = rule_based_evaluation(
        df,
        gene_symbol,
        primary_gene_symbol,
        gene_name,
    )

    if rule_result.matched_rule is not None:
        return AlternateAbbreviationPredictionResult(
            matched_rule=rule_result.matched_rule,
            lcs_similarity_score=rule_result.lcs_similarity_score,
            num_extra_characters=rule_result.num_extra_characters,
            conflicting_category=rule_result.conflicting_category,
            llm_annotation=None,
            error_message=None,
        )

    payload = prompt.build_payload(
        gene_symbol=gene_symbol,
        primary_gene_symbol=primary_gene_symbol,
        gene_name=gene_name,
        hgnc_id=hgnc_id,
    )

    try:
        task_result = task_runner.execute(
            prompt_name=prompt.name,
            prompt_version=prompt.version,
            payload=payload,
            response_model=AlternateAbbreviationPredictionResult,
        )

        return AlternateAbbreviationPredictionResult(
            llm_annotation=task_result.llm_annotation,
            lcs_similarity_score=rule_result.lcs_similarity_score,
        )

    except Exception as e:
        return AlternateAbbreviationPredictionResult(
            error_message=str(e),
            lcs_similarity_score=rule_result.lcs_similarity_score,
        )
    
def run_experiments(
    df: pl.DataFrame,
    temperatures: list[float],
    num_runs: int,
    prompt_version: str,
) -> list[dict[str, Any]]:
    """Run LLM annotation experiments across a range of temperatures and multiple runs per temperature, storing results for analysis.

    :param df: Input dataframe containing gene symbols and associated information.
    :param temperatures: List of temperature values to experiment with.
    :param num_runs: Number of runs to execute for each temperature setting.
    :param prompt_version: Version of the prompt template to use for annotation.
    :return: List of stored runs with LLM outputs and diagnostics for each temperature and run
    """
    stored_runs = []
    prompt = AlternateAbbreviationPrompt(version=prompt_version)

    for temp in temperatures:
        for run_idx in range(num_runs):
            task_runner = build_llm_task_runner(
                MODEL_ID,
                REGION_NAME,
                PROFILE_NAME,
                MAX_TOKENS,
                temp,
            )

            print(f"Running temp={temp}, run={run_idx + 1}")

            results = []

            for row in tqdm(
                df.iter_rows(named=True),
                total=df.height,
                desc=f"T={temp}, run={run_idx + 1}",
                leave=False,
            ):
                result = get_alt_abbreviation_annotation(
                    df,
                    task_runner=task_runner,
                    prompt=prompt,
                    gene_symbol=row["gene_symbol"],
                    primary_gene_symbol=row["primary_gene_symbol"],
                    gene_name=row["gene_name"],
                    hgnc_id=row["HGNC_ID"],
                )

                if result.error_message:
                    error_msg = (
                        f"Error in temp={temp}, run={run_idx + 1}\n"
                        f"gene_symbol={row['gene_symbol']}\n"
                        f"primary_gene_symbol={row['primary_gene_symbol']}\n"
                        f"HGNC_ID={row['HGNC_ID']}\n"
                        f"error={result.error_message}"
                    )

                    raise RuntimeError(error_msg)

                results.append(
                    {
                        "llm_annotation": result.llm_annotation,
                        "matched_rule": result.matched_rule,
                        "error_message": result.error_message,
                        "lcs_similarity_score": result.lcs_similarity_score,
                        "gt": row["alternate_abbreviation_status"],
                    }
                )
            print("prompt_version being stored:", prompt_version)

            stored_runs.append(
                {
                    "run_idx": run_idx,
                    "prompt_version": prompt_version,
                    "temperature": temp,
                    "results": results,
                }
            )

            print(f"Done temp={temp}, run={run_idx + 1}")

    return stored_runs

def build_experiment_key(
        sample_name: str,
        prompt_version: str,
        temperatures: list[float],
        num_runs: int) -> str:
    """Build a unique key for identifying an experiment configuration based on sample name, prompt version, temperatures, and number of runs.

    :param sample_name: Name of the sample or dataset being used.
    :param prompt_version: Version of the prompt template used in the experiment.
    :param temperatures: List of temperature values used in the experiment.
    :param num_runs: Number of runs executed for each temperature setting.
    :return: A string key uniquely identifying the experiment configuration.
    """
    temp_str = "-".join(str(t).replace(".", "p") for t in temperatures)
    return f"{sample_name}_p{prompt_version}_t{temp_str}_agg{num_runs}"

def create_analysis_summary(cm: pl.DataFrame) -> pl.DataFrame:
    """Compute per-class recall-style summary from a confusion matrix.

    :param cm: Confusion matrix (square DataFrame)
    :return: Summary DataFrame per class
    """
    rows = []

    classes = cm.columns[1:]

    for cls in classes:
        row = cm.filter(pl.col("gt") == cls)

        if row.height == 0:
            numerator = 0
            denominator = 0
        else:
            numerator = row.select(pl.col(cls)).item()

            denominator = row.select(pl.sum_horizontal(classes)).item()

        rows.append(
            {
                "consensus_w_curator": cls,
                "numerator": numerator,
                "denominator": denominator,
                "percentage": (
                    (numerator or 0) / denominator * 100
                    if denominator is not None and denominator > 0
                    else 0.0
                ),
            }
        )

    return pl.DataFrame(rows)

def compute_overall_accuracy(cm: pl.DataFrame) -> float:
    """Compute overall accuracy from a confusion matrix.

    :param cm: Confusion matrix DataFrame (square matrix)
    :return: Accuracy as a float between 0 and 1
    """
    classes = cm.columns[1:]

    total = cm.select(pl.sum_horizontal(classes)).to_series().sum()

    correct = 0

    for cls in classes:
        row = cm.filter(pl.col("gt") == cls)

        if row.height > 0:
            correct += row.select(pl.col(cls)).item() or 0

    return correct / total if total > 0 else math.nan

def boolean_confusion_matrix(
    df: pl.DataFrame,
    gt_col: str,
    pred_col: str
) -> pl.DataFrame:
    """Compute a boolean confusion matrix in the traditional matrix style.

    :param df: pandas DataFrame
    :param gt_col: str, name of the ground truth column
    :param pred_col: str, name of the predicted column
    :return: confusion_matrix- pandas DataFrame with row/column labels
    """
    df_clean = df.filter(
        pl.col(gt_col).is_not_null() & pl.col(pred_col).is_not_null()
    ).with_columns(
        [
            pl.col(gt_col).cast(pl.Boolean).alias("gt"),
            pl.col(pred_col).cast(pl.Boolean).alias("pred"),
        ]
    )

    return (
        df_clean.group_by(["gt", "pred"])
        .len()
        .pivot(
            values="len",
            index="gt",
            on="pred",
            aggregate_function="sum",
        )
        .fill_null(0)
    )

def create_system_level_predictions(
    df: pl.DataFrame,
    pred_col: str,
    system_pred_col: str = "llm_system",
) -> pl.DataFrame:
    """Create a system-level prediction column where null predictions
    are treated as False.

    :param df: Input dataframe
    :param pred_col: Original prediction column
    :param system_pred_col: Name of new system-level prediction column
    :return: DataFrame with added system-level prediction column
    """
    return df.with_columns(pl.col(pred_col).fill_null(False).alias(system_pred_col))

def compute_coverage(
    df: pl.DataFrame,
    pred_col: str,
) -> float:
    """Compute proportion of rows with non-null predictions. Low coverage is not inherently bad, as it just shows the ration of samples that were evaluated by rules based methods compared to by LLM.

    :param df: Input dataframe
    :param pred_col: Original prediction column
    :return: Proportion of rows with llm predictions
    """
    return df.select(pl.col(pred_col).is_not_null().mean()).item()

def analyze_results(
    df: pl.DataFrame,
) -> RunResult:
    """Evaluate LLM predictions against ground truth using:
    1. Conditional LLM-only evaluation
    2. System-level evaluation (null -> False)

    :param df: Input dataframe
    :return: Tuple containing:
    - LLM match analysis summary DataFrame
    - LLM overall accuracy (float)
    - LLM precision/recall/f1 metrics DataFrame
    - LLM coverage (float)
    - System-level match analysis summary DataFrame
    - System-level overall accuracy (float)
    - System-level precision/recall/f1 metrics DataFrame
    """
    analysis_df = df.clone()

    # CONDITIONAL LLM EVALUATION

    llm_cm = boolean_confusion_matrix(
        analysis_df,
        "gt",
        "llm",
    )

    tn = llm_cm.filter(~pl.col("gt"))["false"].item()
    fp = llm_cm.filter(~pl.col("gt"))["true"].item()
    fn = llm_cm.filter(pl.col("gt"))["false"].item()
    tp = llm_cm.filter(pl.col("gt"))["true"].item()

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    llm_metrics_df = pl.DataFrame(
        {
            "evaluation_type": ["conditional_llm"],
            "precision": [precision],
            "recall": [recall],
            "f1": [f1],
            "TP": [tp],
            "FP": [fp],
            "TN": [tn],
            "FN": [fn],
        }
    )

    llm_match_analysis_summary = create_analysis_summary(llm_cm)

    llm_accuracy = compute_overall_accuracy(llm_cm)

    llm_coverage = compute_coverage(
        analysis_df,
        "llm",
    )

    # SYSTEM-LEVEL EVALUATION (null -> False)

    system_df = create_system_level_predictions(
        analysis_df,
        pred_col="llm",
        system_pred_col="llm_system",
    )

    system_cm = boolean_confusion_matrix(
        system_df,
        "gt",
        "llm_system",
    )

    system_tn = system_cm.filter(~pl.col("gt"))["false"].item()
    system_fp = system_cm.filter(~pl.col("gt"))["true"].item()
    system_fn = system_cm.filter(pl.col("gt"))["false"].item()
    system_tp = system_cm.filter(pl.col("gt"))["true"].item()

    system_precision = (
        system_tp / (system_tp + system_fp) if (system_tp + system_fp) > 0 else 0.0
    )

    system_recall = (
        system_tp / (system_tp + system_fn) if (system_tp + system_fn) > 0 else 0.0
    )

    system_f1 = (
        2 * system_precision * system_recall / (system_precision + system_recall)
        if (system_precision + system_recall) > 0
        else 0.0
    )

    system_metrics_df = pl.DataFrame(
        {
            "evaluation_type": ["system_level"],
            "precision": [system_precision],
            "recall": [system_recall],
            "f1": [system_f1],
            "TP": [system_tp],
            "FP": [system_fp],
            "TN": [system_tn],
            "FN": [system_fn],
        }
    )

    system_match_analysis_summary = create_analysis_summary(system_cm)

    system_accuracy = compute_overall_accuracy(system_cm)

    return RunResult(
        llm_accuracy=llm_accuracy,
        llm_coverage=llm_coverage,
        llm_summary=llm_match_analysis_summary,
        llm_metrics=llm_metrics_df,
        system_accuracy=system_accuracy,
        system_summary=system_match_analysis_summary,
        system_metrics=system_metrics_df,
    )

def add_correctness_column(df: pl.DataFrame) -> pl.DataFrame:
    """Add a boolean column indicating whether the LLM prediction matches the ground truth.

    :param df: Input DataFrame with 'gt' and 'llm' columns
    :return: DataFrame with added 'llm_correct' column
    """
    return df.with_columns(
        (pl.col("llm") == pl.col("alternate_abbreviation_status")).alias("llm_correct")
    )

def build_eval_df(
    df: pl.DataFrame,
    stored_runs: list,
    temperature: float,
    run_idx: int,
) -> pl.DataFrame:
    """Build evaluation DataFrame by attaching LLM outputs and diagnostics to the original dataframe for a specific run.

    :param df: Original input dataframe
    :param stored_runs: List of stored runs containing LLM outputs
    :param temperature: Temperature of the run to select
    :param run_idx: Index of the run to use from stored_runs (default: -1 for last run)
    :return: DataFrame with LLM outputs and diagnostics attached, plus correctness column
    """
    run = next(
        (
            r
            for r in stored_runs
            if r["temperature"] == temperature and r["run_idx"] == run_idx
        ),
        None,
    )

    if run is None:
        available = [(r["temperature"], r["run_idx"]) for r in stored_runs]

        error_msg = (
            f"No run found for (temperature={temperature}, run_idx={run_idx}). "
            f"Available runs: {available}"
        )

        raise ValueError(error_msg)

    results = run["results"]

    eval_rows = []

    for row, res in zip(df.iter_rows(named=True), results, strict=True):
        eval_rows.append(
            {
                **row,
                "llm": res.get("llm_annotation"),
                "matched_rule": res.get("matched_rule"),
                "temperature": temperature,
                "run_idx": run_idx,
            }
        )

    eval_df = pl.DataFrame(eval_rows)

    return add_correctness_column(eval_df)


def annotate_alt_abbrev(
    df: pl.DataFrame,
    temperature: float,
    prompt_version: str,
) -> pl.DataFrame:

    prompt = AlternateAbbreviationPrompt(version=prompt_version)

    task_runner = build_llm_task_runner(
        MODEL_ID,
        REGION_NAME,
        PROFILE_NAME,
        MAX_TOKENS,
        temperature=temperature,
    )

    predictions = []
    skip_reasons = []

    for row in tqdm(df.iter_rows(named=True), total=df.height):

        result = get_alt_abbreviation_annotation(
            df=df,
            task_runner=task_runner,
            prompt=prompt,
            gene_symbol=row["gene_symbol"],
            primary_gene_symbol=row["primary_gene_symbol"],
            gene_name=row["gene_name"],
            hgnc_id=row["HGNC_ID"],
        )

        if result.error_message:
            raise RuntimeError(result.error_message)

        predictions.append(result.llm_annotation)

        if result.matched_rule is not None:
            skip_reasons.append(str(result.matched_rule))

        elif result.llm_annotation is None and result.matched_rule is None:
            skip_reasons.append("LLM_ERROR")

        else:
            skip_reasons.append(None)

    return df.with_columns([
        pl.Series("Alternate Abbreviation Symbol", predictions),
        pl.Series("skip_reason", skip_reasons),
    ])