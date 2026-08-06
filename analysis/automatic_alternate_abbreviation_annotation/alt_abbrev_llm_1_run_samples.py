"""Automatic Alternate Abbreviation Annotation:** LLM predictions on gene-alias pairs for Alternate Abbreviation alias symbols

This notebook runs the Alternate Abbreviation annotation workflow on a dataset of gene alias symbols.
Alias symbols are first filtered using heuristic rules, after which those requiring further review are
submitted to the LLM for annotation. Results are cached and stored for downstream analysis. For more
information on Alternate Abbreviations and the overall workflow, see [README](alternate_abbreviation_README).
"""

import argparse
import pickle
from pathlib import Path

import alt_abbrev_llm_functions as aalfx
import polars as pl

ALT_ABBREV_ROOT = Path(__file__).resolve().parent
ALT_ABBREV_OUTPUT_PATH = ALT_ABBREV_ROOT / "output"
# Load dataset with gene-alias pairs manually curated for Alternate Abbreviation alias symbols
DEFAULT_INPUT_PATH = (
    ALT_ABBREV_OUTPUT_PATH
    / "alt_abbrev_annotation_manually_annotated_df.xlsx"
)

def load_cached_runs(
    stored_runs_path: Path,
    expected_metadata: dict,
) -> list[dict] | None:
    """Load cached runs when their metadata matches the requested experiment."""
    if not stored_runs_path.exists():
        return None

    with stored_runs_path.open("rb") as file:
        saved_data = pickle.load(file)

    if saved_data.get("metadata") != expected_metadata:
        print(f"Metadata mismatch in {stored_runs_path}. Re-running experiments.")
        return None

    return saved_data["runs"]

def run_annotation(
    input_path: Path = DEFAULT_INPUT_PATH,
    subset_size: int | None = None,
    temperatures: list[float] | None = None,
    num_runs: int = 3,
    prompt_version: str = "v1",
) -> list[dict]:
    """Run Alternate Abbreviation annotation experiments.

    Args:
        input_path: Excel file containing gene-alias pairs.
        subset_size: Number of rows to process. Use None for the full dataset.
        temperatures: LLM temperatures to evaluate.
        num_runs: Number of runs per temperature.
        prompt_version: Prompt version passed to the annotation workflow.

    Returns:
        The stored experiment runs.
    """
    if subset_size is not None and subset_size <= 0:
        raise ValueError("subset_size must be greater than zero.")

    if temperatures is None:
        temperatures = [0.0]

    llm_runs_path = ALT_ABBREV_OUTPUT_PATH / "llm_runs"
    llm_runs_path.mkdir(parents=True, exist_ok=True)

    df = pl.read_excel(input_path)

    if subset_size is not None:
        df = df.head(subset_size)
        sample_name = f"subset_{input_path.stem}"
    else:
        sample_name = input_path.stem

    experiment_key = aalfx.build_experiment_key(
        sample_name,
        prompt_version,
        temperatures,
        num_runs,
    )

    stored_runs_path = llm_runs_path / f"stored_runs_{experiment_key}.pkl"

    metadata = {
        "sample_path": str(input_path),
        "sample_name": sample_name,
        "prompt_version": prompt_version,
        "temperatures": temperatures,
        "num_runs": num_runs,
    }

    stored_runs = load_cached_runs(stored_runs_path, metadata)

    if stored_runs is None:
        print(f"Running experiments on {df.height} rows...")

        stored_runs = aalfx.run_experiments(
            df,
            temperatures,
            num_runs,
            prompt_version,
        )

        with stored_runs_path.open("wb") as file:
            pickle.dump(
                {
                    "metadata": metadata,
                    "runs": stored_runs,
                },
                file,
            )
    else:
        print(f"Loading cached runs from {stored_runs_path}")

    return stored_runs

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Alternate Abbreviation LLM annotation experiments."
    )

    parser.add_argument(
        "--input-path",
        type=Path,
        default=DEFAULT_INPUT_PATH,
        help="Path to the manually annotated Excel dataset.",
    )

    parser.add_argument(
        "--subset-size",
        type=int,
        default=None,
        help="Process only the first N rows. Omit to process the full dataset.",
    )

    parser.add_argument(
        "--num-runs",
        type=int,
        default=3,
        help="Number of LLM runs per temperature.",
    )

    parser.add_argument(
        "--prompt-version",
        default="v1",
        help="Prompt version to use.",
    )

    return parser.parse_args()

def main() -> None:
    """Run the workflow from the command line."""
    args = parse_args()

    run_annotation(
        input_path=args.input_path,
        subset_size=args.subset_size,
        num_runs=args.num_runs,
        prompt_version=args.prompt_version,
    )


if __name__ == "__main__":
    main()