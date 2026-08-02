"""Automatic Alternate Abbreviation Annotation:** LLM predictions on gene-alias pairs for Alternate Abbreviation alias symbols

This notebook runs the Alternate Abbreviation annotation workflow on a dataset of gene alias symbols.
Alias symbols are first filtered using heuristic rules, after which those requiring further review are
submitted to the LLM for annotation. Results are cached and stored for downstream analysis. For more
information on Alternate Abbreviations and the overall workflow, see [README](alternate_abbreviation_README).
"""

import pickle
from pathlib import Path

import alt_abbrev_llm_functions as aalfx
import polars as pl

ALT_ABBREV_ROOT = Path(__file__).resolve().parent
ALT_ABBREV_OUTPUT_PATH = ALT_ABBREV_ROOT / "output"

# A test is recommended before using resources to run full dataset

# Change to True if running a test
RUN_SUBSET = False

# Load dataset with gene-alias pairs manually curated for Alternate Abbreviation alias symbols
df = pl.read_excel(
    ALT_ABBREV_OUTPUT_PATH / "alt_abbrev_annotation_manually_annotated_df.xlsx"
)

# Create a truncated version of the dataset to test
if RUN_SUBSET:
    test_df = df.head(30)
    SAMPLE_PATH = Path(
        ALT_ABBREV_OUTPUT_PATH
        / "subset_alt_abbrev_annotation_manually_annotated_df.xlsx"
    )
    test_df.write_excel(SAMPLE_PATH)
    df = test_df
else:
    SAMPLE_PATH = Path(
        ALT_ABBREV_OUTPUT_PATH / "alt_abbrev_annotation_manually_annotated_df.xlsx"
    )

# Run LLM with gene symbols, name, and prompt
TEMPERATURES = [0.8]
NUM_RUNS = 3
PROMPT_VERSION = "v1"

sample_name = SAMPLE_PATH.stem
temp_str = "-".join(str(t).replace(".", "p") for t in TEMPERATURES)

experiment_key = aalfx.build_experiment_key(
    sample_name,
    PROMPT_VERSION,
    TEMPERATURES,
    NUM_RUNS,
)

stored_runs_path = (
    ALT_ABBREV_OUTPUT_PATH / "llm_runs" / f"stored_runs_{experiment_key}.pkl"
)

metadata = {
    "sample_path": str(SAMPLE_PATH),
    "sample_name": sample_name,
    "prompt_version": PROMPT_VERSION,
    "temperatures": TEMPERATURES,
    "num_runs": NUM_RUNS,
}
if stored_runs_path.exists():
    with stored_runs_path.open("rb") as f:
        saved_data = pickle.load(f)

    if saved_data.get("metadata") == metadata:
        print(f"Loading {stored_runs_path}")
        stored_runs = saved_data["runs"]

    else:
        print(f"Metadata mismatch in {stored_runs_path}. Re-running experiments.")

        stored_runs = aalfx.run_experiments(df, TEMPERATURES, NUM_RUNS, PROMPT_VERSION)

        with stored_runs_path.open("wb") as f:
            pickle.dump(
                {
                    "metadata": metadata,
                    "runs": stored_runs,
                },
                f,
            )

else:
    print("Running experiments...")

    stored_runs = aalfx.run_experiments(df, TEMPERATURES, NUM_RUNS, PROMPT_VERSION)

    with stored_runs_path.open("wb") as f:
        pickle.dump(
            {
                "metadata": metadata,
                "runs": stored_runs,
            },
            f,
        )

# Summarizes the conditions of the runs stored
sample_names = set()
all_runs = []

for path in (ALT_ABBREV_OUTPUT_PATH / "llm_runs").glob("stored_runs_*.pkl"):
    sample_name = path.stem.replace("stored_runs_", "")
    sample_names.add(sample_name)

    with path.open("rb") as f:
        data = pickle.load(f)

    for run in data["runs"]:
        run["sample_name"] = sample_name
        all_runs.append(run)

unique_runs = {
    (run["sample_name"], run["prompt_version"], run["temperature"], run["run_idx"])
    for run in all_runs
}

print(f"Files used for runs:{sorted(sample_names)}")
print(f"Number of unique runs: {len(unique_runs)}")
print(f"Unique runs: {unique_runs}")

df.write_parquet(stored_runs_path.with_suffix(".parquet"))
