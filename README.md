# CRISPResso-AntiNoise-Wrapper
Statistical multi-project wrapper for CRISPResso2 with background noise filtering and error handling for genome editing data.

## Overview
This pipeline processes CRISPR editing experiments from FASTQ files using project-specific guides, amplicons, and analysis parameters.

#### The program requires:  
FASTQ files placed into the fastq/ directory  
Three CSV configuration tables:  
samples.csv  
sequences.csv  
parameters.csv  

## Directory structure
```
project/  
 │  
 ├── fastq/              	   # Input FASTQ files  
 │  
 ├── samples.csv         	   # Sample annotation table  
 ├── sequences.csv           # Guide and sequence information  
 ├── parameters.csv      	   # Analysis parameters  
 │  
 ├── outcome/  
 │   ├── alleles/                # Alleles frequency tables   
 │   ├── annotation/		   # Annotated alleles frequency tables  
 │   ├── base/			   # Base editing results (matrix)  
 │   ├── processed/		   # Intermediate result  
 │   ├── final_summary_noise.csv # Final summary after noise extraction  
 │   └── results_noise.csv       # Summary after CRISPResso2 analysis  
 │  
 └── report/				   # CRISPresso2 results  
```
## Input files
### 1. FASTQ files
#### All sequencing files must be placed into:  
fastq/  

#### Supported formats:  
.fastq  
.fq  
.fastq.gz  
.fq.gz  

#### Paired-end naming examples:  
sample1_R1.fastq.gz  
sample1_R2.fastq.gz  
sample1_r1.fastq.gz  
sample1_r2.fastq.gz  
sample1_R1_text.fastq.gz  
sample1_R2_text.fastq.gz  

#### Single-end naming example:  
sample1.fastq.gz  
The pipeline automatically detects paired and single-end reads.  

### 2. Parameters
This table contains CRISPResso or pipeline analysis parameters.

#### Format
```
parameters                     | values    | project | change  
--default_min_aln_score        | 60        | 3       | 50  
-q                             | 30        |         |  
-s                             | 0         |         |  
--plot_window_size             | 20        |         |  
--exclude_bp_from_left         | 15        |         |  
--exclude_bp_from_right        | additional| 2       | 15  
--write_cleaned_report         |           |         |  
--place_report_in_output_folder|           |         |  
--conversion_nuc_from          | A         |         |  
--conversion_nuc_to            | G         |         |  
```

#### Column description
#### parameters
CRISPResso2 parameters.

##### values
Default parameter value used for all projects.
If a parameter is project-specific and should not be applied globally, use:
additional
In this case column project should be filled with project number and in column change the parameter value should be provided.

##### project
Project ID where the parameter value differs from the default.

##### change
Alternative value used only for the specified project.
Example:
```
parameters              | values | project | change
--default_min_aln_score | 60     | 3       | 50

```
##### Interpretation:
all projects use 60
project 3 uses 50

### 3. Samples
The table “samples.csv” describes all experimental samples.
Sample names are generated automatically from FASTQ filenames, after which the user fills remaining columns manually.

#### Format
```
sample_name     | guide_id | group | edit_mode | project | noise
CFTR-LRG-802    | 2        | 1     | prime     | 2       | no
CFTR-LRG-803    | 2        | 1     | prime     | 2       | no
CFTR-LRG-835    |          | 0     | prime     | 2       | no
LRG-DMD-1201_L01| 4        | 1     | cas9      | 3       | no
mut-117_S10_L001| 1        | 1     | base      | 1       | no
```
#### Column description
##### sample_name
Sample identifier.
Automatically generated from FASTQ filenames. The program asks if user needs table to be generated before starting the pipeline.

##### guide_id
Guide RNA ID corresponding to the guide entry in sequences.csv.
Controls may be left empty.
If empty, the control sample will be processed against all guides within the project.

##### group
Biological or technical replicate group.
identical group numbers = replicates
0 = control sample

##### edit_mode
Editing mode used for the sample.  
Supported modes:
cas9
base
prime

##### project
Project identifier.
Each project contains its own guides, amplicons, and parameters.

##### project
Project identifier.
Each project contains its own guides, amplicons, and parameters.

### 4. Sequences
The table sequences.csv contains guide RNA and amplicon information.

#### Format
```
guide_id | spacer          | scaffold          | extension           | project_id | amplicon          | expected_seq
1        | spacer sequence |                   |                     | 1          | amplicon sequence | expected sequence
2        | spacer sequence | scaffold sequence | extension sequence  | 2          | amplicon sequence | expected sequence
3        | spacer sequence |                   |                     | 3          | amplicon sequence |
4        | spacer sequence |                   |                     |            |                   |

```

#### Column description
##### guide_id
Unique guide identifier.
Referenced from samples.csv.

##### spacer
Guide spacer sequence.
Required for all projects.

##### scaffold
pegRNA scaffold sequence.
Leave empty for non-prime editing projects.

##### extension
pegRNA extension sequence.
Leave empty for non-prime editing projects.

##### project_id
This and further columns are not associated with the guides. Project ids are needed for referencing amplicon and expected sequences.

##### amplicon
Reference amplicon sequence for the project.

##### expected_seq
Expected edited sequence for the project. Can be empty if no expected edited sequence is defined.

## Running the pipeline
The pipeline is built on CRISPResso2. For the best experience and to prevent environment issues, please install it within a separate Conda environment:  
Tool: CRISPResso2  
Env Setup: conda create -n crispresso2_env -c bioconda crispresso2  
Activation: conda activate crispresso2_env  

### Step 1. Prepare input files
Place:
FASTQ files into fastq/
CSV tables into the project root directory

### Step 2. Generate sample table (optional)
The pipeline can automatically generate samples.csv from FASTQ filenames.
python main.py
Choose:
Create samples.csv? (y/n)
After generation, manually fill remaining columns.

### Step 3. Run analysis
python main.py

## Output
Results are written into:  
outcome/  

Main outputs:  
results_noise.csv: Main analysis results  
final_summary_noise.csv:  Noise/control analysis

Reports after CRISPResso2 are written into:
report/
