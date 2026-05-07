from pathlib import Path
import re

READ_REGEX = re.compile(
    r"(?P<sample>.+?)(?:[ _\-.]r(?P<read>[12]))(?=.*\.(fastq|fq))|(?P<single>.+?)(?=\.(fastq|fq))",
    re.IGNORECASE
)

RESULTS_PATH = Path("outcome/results.csv")
RESULTS_PATH_NOISE = Path("outcome/results_noise.csv")
REPORT_DIR = "report"
ANNOTATION_DIR = Path("outcome/annotation")