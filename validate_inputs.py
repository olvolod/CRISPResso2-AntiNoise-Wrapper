"""Preflight validation for FASTQ inputs and CSV tables.

Adjust the paths in main.py (or set the FASTQ_DIR environment variable)
to match your project layout.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


SUPPORTED_FASTQ_SUFFIXES = (
    ".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.gzip", ".fq.gzip"
)

SAMPLES_REQUIRED = {
    "sample_name", "guide_id", "group", "edit_mode", "project", "noise"
}
SEQUENCES_REQUIRED = {
    "guide_id",
    "spacer",
    "scaffold",
    "extension",
    "amplicon",
    "project_id",
}
PARAMETERS_REQUIRED = {"parameters", "values", "project", "change"}


def ask_yes_no(question: str) -> bool:
    """Ask until the user enters y/yes or n/no."""
    while True:
        answer = input(f"{question} (y/n): ").strip().lower()
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please enter 'y' or 'n'.")


def read_csv(path: Path) -> pd.DataFrame:
    """Read comma- or semicolon-separated CSV files as strings."""
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")
    # sep=None lets Python's CSV sniffer detect comma vs semicolon.
    df = pd.read_csv(
        path,
        sep=None,
        engine="python",
        dtype=str,
        keep_default_na=False,
        index_col=False,
    )
    return df


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names without removing characters."""
    df = df.copy()



    df.columns = [
        str(col).strip().lstrip("\ufeff")
        for col in df.columns
    ]




    for col in df.columns:
        df[col] = df[col].map(
            lambda value: value.strip()
            if isinstance(value, str)
            else value
        )

    return df


def check_table_structure(
    df: pd.DataFrame,
    path: Path,
    required_columns: set[str],
    *,
    allow_extra_columns: bool = False,
    known_optional_columns: Optional[set[str]] = None,
) -> bool:
    """Check required/extra columns and completely empty rows."""
    ok = True
    optional = known_optional_columns or set()

    missing = sorted(required_columns - set(df.columns))
    extra = sorted(set(df.columns) - required_columns - optional)

    if missing:
        print(f"ERROR: {path.name} is missing required column(s): {', '.join(missing)}.")
        ok = False

    if extra and not allow_extra_columns:
        print(
            f"WARNING: {path.name} contains unexpected column(s): "
            f"{', '.join(extra)}. Please check the table for accidental extra columns."
        )

    empty_rows = df.index[
        df.apply(lambda row: all(not str(v).strip() for v in row), axis=1)
    ].tolist()
    if empty_rows:
        row_numbers = ", ".join(str(index + 2) for index in empty_rows)
        print(f"WARNING: {path.name} contains completely empty row(s): {row_numbers}.")
        if not ask_yes_no("Have you checked and removed unintended empty rows?"):
            ok = False

    return ok


def find_fastq_files(directory: Path) -> list[Path]:
    """Find supported FASTQ files recursively, case-insensitively."""
    if not directory.is_dir():
        return []
    found = []
    for path in directory.rglob("*"):
        if path.is_file() and path.name.lower().endswith(SUPPORTED_FASTQ_SUFFIXES):
            found.append(path)
    return found


def validate_fastq_directory(fastq_dir: Path) -> bool:
    if not fastq_dir.exists():
        print(f"ERROR: FASTQ directory does not exist: {fastq_dir}")
        print("Please set FASTQ_DIR to the directory containing your FASTQ files.")
        return False
    if not fastq_dir.is_dir():
        print(f"ERROR: FASTQ path is not a directory: {fastq_dir}")
        return False

    files = find_fastq_files(fastq_dir)
    if not files:
        print(
            f"ERROR: No supported FASTQ files were found in '{fastq_dir}'. "
            "Supported extensions: .fastq, .fq, .fastq.gz, .fq.gz, "
            ".fastq.gzip and .fq.gzip."
        )
        return False

    print(f"FASTQ check passed: found {len(files)} supported file(s) in {fastq_dir}.")
    return True

PARAMETERS_REQUIRED = {"parameters", "values"}
PARAMETERS_OPTIONAL = {"project", "change"}
def validate_parameters(path: Path) -> bool:
    try:
        df = normalize_columns(read_csv(path))
        print("Detected columns:", [repr(c) for c in df.columns])
    except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
        print(f"ERROR: Could not read parameters table: {exc}")
        return False

    ok = check_table_structure(
        df, path, PARAMETERS_REQUIRED,
        known_optional_columns=PARAMETERS_OPTIONAL,
    )

    # Excel may store localized or English formula/name errors as text.
    error_pattern = re.compile(r"#(?:NAME\?|ИМЯ\?)", re.IGNORECASE)
    bad_cells = []
    for row_idx, row in df.iterrows():
        for column, value in row.items():
            if error_pattern.search(str(value)):
                bad_cells.append(f"row {row_idx + 2}, column '{column}'")

    if bad_cells:
        print(
            "ERROR: Excel error value(s) #NAME? found in parameters table "
            f"({'; '.join(bad_cells)}). Please check and correct the affected cells."
        )
        ok = False

    parameter_column = next(
        (column for column in df.columns if column.lower() in {"parameters", "parameter", "param", "name"}),
        None,
    )
    if parameter_column is None:
        print(
            "ERROR: Could not identify the parameter-name column in parameters.csv. "
            "Expected a column named 'parameters'."
        )
        return False

    # --- check --suppress_plots ---
    names = [str(value).strip() for value in df[parameter_column]]
    suppress_count = names.count("--suppress_plots")

    if suppress_count == 0:
        print("ERROR: '--suppress_plots' is missing from parameters.csv.")
        print(
            "This flag is REQUIRED: without it, CRISPResso2 may attempt to generate plots "
            "and crash in headless environments (e.g., Google Colab), stopping the entire pipeline."
        )
        if ask_yes_no("Add '--suppress_plots' to parameters.csv automatically?"):
            new_row = {column: "" for column in df.columns}
            new_row[parameter_column] = "--suppress_plots"
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            df.to_csv(
                path,
                index=False,
                sep=";" if path.read_text(errors="ignore").count(";") >= path.read_text(errors="ignore").count(",") else ",",
            )
            print("Added '--suppress_plots' to parameters.csv.")
        else:
            print("Please add '--suppress_plots' to parameters.csv before continuing.")
            ok = False
    elif suppress_count > 1:
        print(
            "WARNING: '--suppress_plots' appears more than once in parameters.csv. "
            "Duplicate flags may cause unexpected behavior. Please keep only one occurrence."
        )
        if not ask_yes_no("Have you checked and removed duplicate '--suppress_plots' entries?"):
            ok = False

    return ok


def validate_samples(path: Path) -> tuple[bool, Optional[pd.DataFrame]]:
    try:
        df = normalize_columns(read_csv(path))
    except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
        print(f"ERROR: Could not read samples table: {exc}")
        return False, None

    ok = check_table_structure(df, path, SAMPLES_REQUIRED)
    if not SAMPLES_REQUIRED.issubset(df.columns):
        return False, df

    extension_pattern = re.compile(
        r"(?:\.fastq|\.fq)(?:\.(?:gz|gzip))?$",
        re.IGNORECASE,
    )
    bad_names = [
        (index + 2, value)
        for index, value in df["sample_name"].items()
        if extension_pattern.search(str(value).strip())
    ]
    if bad_names:
        examples = ", ".join(
            f"row {row}: '{name}'"
            for row, name in bad_names[:8]
        )
        print(
            "ERROR: sample_name contains FASTQ file extensions "
            f"({examples})."
        )
        print(
            "Please remove file extensions from sample_name in samples.csv. "
            "If needed, regenerate the samples table and fill in the metadata."
        )
        print("Pipeline stopped. Please correct the table and run again.")
        return False, df


    # For group 0, a guide_id of literal 0 is treated as empty.
    group_values = df["group"].astype(str).str.strip()
    guide_values = df["guide_id"].astype(str).str.strip()
    control_zero_mask = group_values.eq("0") & guide_values.eq("0")
    if control_zero_mask.any():
        df.loc[control_zero_mask, "guide_id"] = ""
        df.to_csv(
            path,
            index=False,
            sep=";" if path.read_text(errors="ignore").count(";") >= path.read_text(errors="ignore").count(",") else ",",
        )
        print(
            f"Corrected {int(control_zero_mask.sum())} control row(s): "
            "guide_id was changed from 0 to an empty cell where group = 0."
        )

    # A control row (group 0) with a non-empty, non-zero guide ID may be intentional.
    group_values = df["group"].astype(str).str.strip()
    guide_values = df["guide_id"].astype(str).str.strip()
    unusual_controls = df[
        group_values.eq("0") & guide_values.ne("") & guide_values.ne("0")
    ]
    if not unusual_controls.empty:
        guide_ids = sorted(set(unusual_controls["guide_id"]))
        print(
            "WARNING: Control rows (group = 0) have non-empty guide_id value(s): "
            f"{', '.join(guide_ids)}."
        )
        if not ask_yes_no(
            "Do you intentionally want to run controls only for the listed guide ID(s)?"
        ):
            print("Please correct samples.csv. The pipeline will not be started.")
            return False, df

    # Check whether every project with experimental samples also has controls.
    nonzero_groups = df[~group_values.eq("0")]
    control_rows = df[group_values.eq("0")]
    experimental_projects = {
        value for value in nonzero_groups["project"].astype(str).str.strip() if value
    }
    control_projects = {
        value for value in control_rows["project"].astype(str).str.strip() if value
    }
    missing_control_projects = sorted(experimental_projects - control_projects)
    if missing_control_projects:
        print(
            "WARNING: No group = 0 control samples were found for project(s): "
            f"{', '.join(missing_control_projects)}."
        )
        if not ask_yes_no("Are controls intentionally missing for these project(s)?"):
            print("Please add the missing controls to samples.csv. The pipeline will not be started.")
            return False, df

    return ok, df

SEQUENCES_OPTIONAL = {"expected_seq"}


def validate_sequences(path: Path) -> tuple[bool, Optional[pd.DataFrame]]:
    """Validate sequences.csv, allowing optional fields where appropriate."""

    try:
        df = normalize_columns(read_csv(path))
    except (FileNotFoundError, pd.errors.ParserError, UnicodeError) as exc:
        print(f"ERROR: Could not read sequences table: {exc}")
        return False, None

    ok = check_table_structure(
        df, path, SEQUENCES_REQUIRED,
        known_optional_columns=SEQUENCES_OPTIONAL,
    )

    if not SEQUENCES_REQUIRED.issubset(df.columns):
        return False, df

    # Only check projects that actually have a non-empty project_id.
    # Rows with empty project_id (e.g. "free" guides) are allowed to have
    # empty amplicon / expected_seq.
    project_mask = df["project_id"].astype(str).str.strip().ne("")
    df_with_project = df[project_mask]

    if not df_with_project.empty:
        amplicon_by_project = df_with_project.groupby("project_id")["amplicon"].apply(
            lambda values: values.astype(str).str.strip().ne("").any()
        )
        projects_without_amplicon = amplicon_by_project[
            ~amplicon_by_project
        ].index.tolist()

        if projects_without_amplicon:
            project_names = ", ".join(str(p) for p in projects_without_amplicon)
            print(
                f"ERROR: No amplicon sequence was provided for project(s): "
                f"{project_names}."
            )
            print(
                "Each project must have at least one non-empty amplicon. "
                "Other amplicon cells within the same project may remain empty."
            )
            ok = False

    # expected_seq column is optional, but if it exists, empty cells are allowed
    # only if at least one expected_seq is actually filled in.
    if "expected_seq" not in df.columns:
        print("WARNING: The 'expected_seq' column is missing from sequences.csv.")
        if ask_yes_no("Do you need the expected_seq column for this analysis?"):
            print(
                "Please add the 'expected_seq' column to sequences.csv "
                "and run the program again."
            )
            return False, df
        print("Continuing without the optional 'expected_seq' column.")
    else:
        amplicon = df["amplicon"].astype(str).str.strip()
        expected = df["expected_seq"].astype(str).str.strip()

        # Warn only when amplicon is filled but expected_seq next to it is empty.
        missing_expected = amplicon.ne("") & expected.eq("")
        if missing_expected.any():
            row_numbers = ", ".join(
                str(index + 2) for index in df.index[missing_expected]
            )
            print(
                f"WARNING: expected_seq is empty while amplicon is provided "
                f"in row(s): {row_numbers} of {path.name}. "
                "This is allowed, but please confirm it is intentional."
            )
            if not ask_yes_no(
                "Is it intentional that expected_seq is empty for these "
                "amplicon-bearing rows?"
            ):
                ok = False

    print("Sequences table validation completed.")
    return ok, df



def compare_sample_and_sequence_ids(
    samples: pd.DataFrame,
    sequences: pd.DataFrame,
) -> bool:
    """Compare guide IDs and project IDs as independent sets.

    In sequences.csv, guide definitions (guide_id + spacer) and project
    definitions (project_id + amplicon) live independently and may even
    share a row. The guide↔project association exists only in samples.csv,
    so we validate the two dimensions separately.
    """
    required_samples = {"guide_id", "project", "group"}
    required_sequences = {"guide_id", "project_id"}
    if not required_samples.issubset(samples.columns) or \
       not required_sequences.issubset(sequences.columns):
        return False

    sample_guides: set[str] = set()
    sample_projects: set[str] = set()
    for _, row in samples.iterrows():
        group = str(row["group"]).strip()
        guide = str(row["guide_id"]).strip()
        project = str(row["project"]).strip()
        if group == "0" and guide in {"", "0"}:
            if project:
                sample_projects.add(project)
            continue
        if guide:
            sample_guides.add(guide)
        if project:
            sample_projects.add(project)

    spacer = sequences["spacer"].astype(str).str.strip()
    amplicon = sequences["amplicon"].astype(str).str.strip()

    sequence_guides = set(
        sequences.loc[spacer.ne(""), "guide_id"].astype(str).str.strip()
    )
    sequence_guides.discard("")

    sequence_projects = set(
        sequences.loc[amplicon.ne(""), "project_id"].astype(str).str.strip()
    )
    sequence_projects.discard("")

    missing_guides = sorted(sample_guides - sequence_guides)
    missing_projects = sorted(sample_projects - sequence_projects)

    if missing_guides or missing_projects:
        print("ERROR: sequences.csv is missing definitions required by samples.csv.")
        if missing_guides:
            print(
                "guide_id(s) used in samples.csv but not defined in sequences.csv "
                f"(no row with a non-empty spacer): {', '.join(missing_guides)}"
            )
        if missing_projects:
            print(
                "project(s) used in samples.csv but not defined in sequences.csv "
                f"(no row with a non-empty amplicon): {', '.join(missing_projects)}"
            )
        print(
            "Please check that every guide_id used in samples.csv has a row "
            "with spacer in sequences.csv, and every project has a row "
            "with amplicon in sequences.csv."
        )
        return False

    print("Sample/sequence guide_id and project ID checks passed.")
    return True


def validate_inputs(
    fastq_dir: Optional[Path | str] = None,
    parameters_path: Path | str = "parameters.csv",
    samples_path: Path | str = "samples.csv",
    sequences_path: Path | str = "sequences.csv",
) -> bool:
    """Run all preflight checks. Return True only when it is safe to start the pipeline."""
    fastq_dir = Path(fastq_dir or os.environ.get("FASTQ_DIR", "fastq"))
    parameters_path = Path(parameters_path)
    samples_path = Path(samples_path)
    sequences_path = Path(sequences_path)

    print("\\n=== Input validation ===")
    checks = [
        validate_fastq_directory(fastq_dir),
        validate_parameters(parameters_path),
    ]

    samples_ok, samples = validate_samples(samples_path)
    sequences_ok, sequences = validate_sequences(sequences_path)
    checks.extend([samples_ok, sequences_ok])

    if samples is not None and sequences is not None:
        checks.append(compare_sample_and_sequence_ids(samples, sequences))

    if all(checks):
        print("All input checks passed. Starting pipeline.\\n")
        return True

    print("\\nInput validation failed. Please fix the issues above before running the pipeline.")
    return False


if __name__ == "__main__":
    raise SystemExit(0 if validate_inputs() else 1)
