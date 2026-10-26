# CRISPResso2-AntiNoise-Wrapper

Statistical multi-project wrapper for CRISPResso2 with background noise filtering and error handling for genome editing data.

## Overview

This pipeline processes CRISPR editing experiments from FASTQ files using project-specific guides, amplicons, and analysis parameters. It extends CRISPResso2 with:

- Sequential execution (allows resuming after failures)
- Automated fallback for prime editing alignment errors
- Control-based noise modeling (Bayesian prior for small control sets)
- Project-level organization (metadata-driven)

## Requirements

- Linux (or Google Colab)
- Python 3.10+
- [CRISPResso2](https://github.com/pinellolab/CRISPResso2) (v2.2.12 or newer)
- Python packages: `biopython`, `statsmodels`, `matplotlib-inline`

---

## Quick Start with Google Colab (Recommended)

The easiest way to run the pipeline is via Google Colab. No local installation is required.

### Step 1: Set up the environment

Open a new Colab notebook and run the following cell. Note: after this cell, Colab will automatically restart the runtime (this is normal).

```python
!pip install -q condacolab
import condacolab
condacolab.install()
```

### Step 2: Install CRISPResso2 and dependencies

After the runtime restarts, run this cell:

```python
!conda install -c bioconda -c conda-forge crispresso2 -y
!pip install biopython statsmodels matplotlib-inline
import os
os.environ['MPLBACKEND'] = 'Agg'
```

If you encounter a matplotlib backend error, run:
```python
import os
os.environ['MPLBACKEND'] = 'Agg'
```
```
!pip install matplotlib==3.8.4
```
and restart the runtime.

### Step 3: Upload your project
When you clone or download this repository, it already contains **example data**:
- Example `samples.csv`, `sequences.csv`, and `parameters.csv` (filled with demonstration values).
- Example FASTQ files in `fastq/`.
- All `.py` files.

**Before running the pipeline, replace the example data with your own:**
- Delete or overwrite the example FASTQ files in `fastq/` with your own sequencing files.
- Replace `samples.csv`, `sequences.csv`, and `parameters.csv` with your own tables (use the examples as templates for the required columns).

Then:
1. Zip your project folder (containing all the `.py` files from the repository, `fastq/`, `samples.csv`, `sequences.csv`, and `parameters.csv`).
2. Upload the ZIP file to Colab (drag and drop it into the file browser).
3. Unzip it:
   ```python
   !unzip your_project.zip

### Step 4: Run the pipeline

Navigate to the project folder. You can copy the directory path by right-clicking the folder in the file browser and choosing "Copy path", then run:

```python
%cd /content/your_project_folder
%run main.py
```

The pipeline will process your samples and save results to `outcome/`. To download them:

```python
!zip -r results.zip outcome/
from google.colab import files
files.download('results.zip')
```

## Local Installation

The pipeline is built on CRISPResso2. For the best experience and to prevent environment issues, install it within a separate Conda environment:

```bash
conda create -n crispresso2_env -c bioconda crispresso2
conda activate crispresso2_env
pip install biopython statsmodels
```

## Directory Structure

```text
project/
├── fastq/                       # Input FASTQ files
├── samples.csv                  # Sample annotation table
├── sequences.csv                # Guide and sequence information
├── parameters.csv               # Analysis parameters
├── outcome/
│   ├── alleles/                 # Allele frequency tables
│   ├── annotation/              # Annotated allele frequency tables
│   ├── base/                    # Base editing results (matrix)
│   ├── processed/               # Intermediate results
│   ├── final_summary_noise.csv  # Final summary after noise extraction
│   └── results_noise.csv        # Summary after CRISPResso2 analysis
└── report/                      # CRISPResso2 results
```

## Input Files

### 1. FASTQ files

Place all sequencing files into `fastq/`. Supported formats: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`.

Paired-end naming examples: `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`.

### 2. `samples.csv`

Describes all experimental samples. Sample names are generated automatically from FASTQ filenames; the user fills the remaining columns.

Example (from the provided dataset):

| sample_name | guide_id | group | edit_mode | project | noise |
|---|---:|---:|---|---:|---|
| spike_CRISPR_high | 3 | 1 | cas9 | 3 | yes |
| spike_CRISPR_low | 3 | 1 | cas9 | 3 | yes |
| spike_CRISPR_med | 3 | 1 | cas9 | 3 | yes |
| spike_base_high | 1 | 1 | base | 1 | yes |
| spike_base_low | 1 | 1 | base | 1 | yes |
| spike_base_med | 1 | 1 | base | 1 | yes |
| spike_prime_high | 2 | 1 | prime | 2 | yes |
| spike_prime_low | 2 | 1 | prime | 2 | yes |
| spike_prime_med | 2 | 1 | prime | 2 | yes |
| masked_control_1 | — | 0 | cas9 | 3 | yes |
| masked_control_1 | — | 0 | base | 1 | yes |
| masked_control_1 | — | 0 | prime | 2 | yes |

Column descriptions:

| Column | Description |
|---|---|
| `sample_name` | Sample identifier generated from FASTQ filenames. |
| `guide_id` | Guide RNA ID corresponding to `sequences.csv`. Controls may be left empty; the control will be processed against all guides within the project. |
| `group` | `0` = control sample; `1+` = experimental replicates. |
| `edit_mode` | Editing mode: `cas9`, `base`, or `prime`. |
| `project` | Project identifier. |
| `noise` | `yes` = run noise pipeline; `no` = standard mode. |

### 3. `sequences.csv`

This file contains two independent tables stored in the same CSV for convenience:

| Table | Columns | Purpose |
|---|---|---|
| Guide table | `guide_id`, `spacer`, `scaffold`, `extension` | Defines guide RNAs used in the experiment. |
| Project table | `project_id`, `amplicon`, `expected_seq` | Defines reference sequences for each project. |

These two tables are independent: a guide can exist without a project, and a project can exist without any specific guide row (if the amplicon is shared). A single row may contain both guide and project information, but this is not required. The association between guides and projects is defined only in `samples.csv`, not here.

Example (from the provided dataset, with guide 4 added as a standalone guide without a project):

| guide_id | spacer | scaffold | extension | project_id | amplicon | expected_seq |
|---:|---|---|---|---:|---|---|
| 1 | `CTTCGGCCCACCCCGCGTCG` | — | — | 1 | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGGAAGAGGGTTGCGACGACCGCGCACGCTCAC` | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGCAAGAGGGTTGCGACGACCGCGCACGCTCAC` |
| 2 | `CTTCGGCCCACCCCGCGTCG` | `gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc` | `ACCCTCTTTTTACGCGGGGTGGGCC` | 2 | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGGAAGAGGGTTGCGACGACCGCGCACGCTCAC` | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGAAAAAGAGGGTTGCGACGACCGCGCACGCTCAC` |
| 3 | `CTTCGGCCCACCCCGCGTCG` | — | — | 3 | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGGAAGAGGGTTGCGACGACCGCGCACGCTCAC` | `CCATCTCAGAACCACCGGGTCGGTGAGGTGTCAGGGGTTGAACGGCGATTGAATCTTCGGCCCACCCCGCGTCGAAAAAGAGGGTTGCGACGACCGCGCACGCTCAC` |
| 4 | `GCTAGCTAGCTAGCTAGCTA` | — | — | — | — | — |

Column descriptions:

| Column | Description |
|---|---|
| `guide_id` | Unique guide identifier referenced from `samples.csv`; part of the guide table. |
| `spacer` | Guide spacer sequence. Required for all projects; part of the guide table. |
| `scaffold` | pegRNA scaffold sequence (prime editing only); part of the guide table. |
| `extension` | pegRNA extension sequence (prime editing only); part of the guide table. |
| `project_id` | Project identifier; part of the project table. |
| `amplicon` | Reference amplicon sequence; part of the project table. |
| `expected_seq` | Expected edited sequence; can be empty. Part of the project table. |

### 4. `parameters.csv`

Contains CRISPResso2 parameters. Each parameter appears only once; project-specific overrides are defined using the `project` and `change` columns. The parameters table must include --suppress_plots. This flag prevents CRISPResso2 from generating plot files. In headless environments (e.g., Google Colab) or on systems without proper display libraries, plot generation can cause CRISPResso2 to crash, which stops the entire pipeline. If --suppress_plots is missing, the validation step will prompt you to add it automatically.

Example (from the provided dataset):

| parameters | values | project | change |
|---|---:|---:|---:|
| `--default_min_aln_score` | 60 | — | — |
| `-q` | 30 | — | — |
| `-s` | 0 | — | — |
| `--plot_window_size` | 20 | 2 | 10 |
| `--min_bp_quality_or_N` | 0 | — | — |
| `--exclude_bp_from_left` | 15 | 2 | 5 |
| `--exclude_bp_from_right` | 15 | 2 | 5 |
| `-w` | 6 | — | — |
| `-wc` | -3 | — | — |
| `--write_cleaned_report` | — | — | — |
| `--place_report_in_output_folder` | — | — | — |
| `--prime_editing_pegRNA_extension_quantification_window_size` | 4 | — | — |
| `--conversion_nuc_from` | G | — | — |
| `--conversion_nuc_to` | C | — | — |
| `--suppress_plots` | — | — | — |

#### How overrides work

| Situation | Behavior |
|---|---|
| `project` is empty | The parameter is applied to all projects. |
| `project` is filled | The parameter is overridden for the specified project using the value in `change`. |
| `change` is `false` | Removes the parameter for the specified project. |
| `values` is `additional` | Adds a parameter only for the specified project. If a value is needed, put it in `change`. |
| Multiple projects are specified | Separate project IDs with commas (e.g. `3,4`) and corresponding values in `change` with commas (e.g. `5,10`). Values are matched by index. |

In the example above:

| Parameter | Default value | Project-specific override |
|---|---:|---|
| `--plot_window_size` | 20 | Project 2 uses 10. |
| `--exclude_bp_from_left` | 15 | Project 2 uses 5. |
| `--exclude_bp_from_right` | 15 | Project 2 uses 5. |

## Running the Pipeline

### Step 1. Prepare input files

Place FASTQ files into `fastq/` and CSV tables into the project root.

### Step 2. Generate sample table (optional)

The pipeline can automatically generate `samples.csv` from FASTQ filenames:

```bash
python main.py
```

Choose `y` when prompted. After generation, manually fill in the remaining columns.

### Step 3. Run analysis

```bash
python main.py
```

## Output

Results are written into `outcome/`:

| File or directory | Description |
|---|---|
| `outcome/results_noise.csv` | Main analysis results. |
| `outcome/final_summary_noise.csv` | Noise/control analysis summary. |
| `outcome/alleles/` | Allele frequency tables. |
| `outcome/annotation/` | Annotated allele frequency tables. |
| `outcome/base/` | Base editing results (matrix). |
| `outcome/processed/` | Intermediate results. |
| `report/` | CRISPResso2 reports. |

## Known Limitations and Troubleshooting

### 1. CSV encoding and Excel issues

If you edit CSV files in Excel on Windows, they may be saved in Windows-1251 encoding or contain `#NAME?` errors (Excel converts text starting with `-` into formulas).

**Solutions:**

- Save files as **CSV UTF-8 (Comma delimited)** in Excel.
- If a value starts with `-` (e.g., `--suppress_plots`), prefix it with an apostrophe (`'--suppress_plots`). The apostrophe will disappear when saved as CSV.
- Alternatively, the pipeline attempts to auto-detect CP1251 encoding if UTF-8 fails.

### 2. Short quantification windows (fallback samples)

For prime editing fallback samples, CRISPResso2 determines the quantification window. If the window is short (4 nt or less), it may map to multiple locations in the amplicon, causing misannotation.

**Solutions:**

- Use `--prime_editing_pegRNA_extension_quantification_window_size` of at least 6 nt when possible.
- The pipeline includes an automatic realignment procedure: if empty rows are detected, it searches for alternative window positions where all modified reads have variants.
- If realignment fails, a warning is printed. Check the window size and extension sequence.

### 3. Missing --suppress_plots in parameters.csv
The pipeline requires --suppress_plots to be present in parameters.csv. Without this flag, CRISPResso2 may attempt to generate plots, which can crash the process in headless environments (e.g., Google Colab) or on systems lacking display libraries, stopping the entire pipeline.

**Solutions:**

Add a row with parameters = --suppress_plots and leave values, project, and change empty. The pipeline’s input validation will detect the missing flag and offer to add it automatically. If you prefer to add it manually, ensure it appears only once in the table.
