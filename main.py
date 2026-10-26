from orchestrator import run_pipeline
from generate_tables import generate_samples_table
from config import RESULTS_PATH
from pathlib import Path
from validate_inputs import validate_inputs

def ask_create_samples():
    while True:
        answer = input("Create samples.csv? (y/n): ").strip().lower()
        if answer in ("y", "yes"):
            return True
        elif answer in ("n", "no"):
            return False
        else:
            print("Please enter 'y' or 'n'")


if __name__ == "__main__":
    if ask_create_samples():
        print("Generating samples table...")

        out_dir = Path.cwd()
        generate_samples_table(out_dir=out_dir)

        print("Done. Please edit the file before running pipeline.")
        exit(0)

    if not validate_inputs():
        raise SystemExit("Pipeline stopped. Please fix the issues above.")

    run_pipeline()
