import pandas as pd
from pathlib import Path
from config import READ_REGEX, RESULTS_PATH, RESULTS_PATH_NOISE, REPORT_DIR, ANNOTATION_DIR
def find_fastq_files(fastq_dir="fastq"):
    fastq_dir = Path(fastq_dir)
    samples = {}

    for f in fastq_dir.iterdir():
        if not f.is_file():
            continue

        m = READ_REGEX.search(f.name)
        if not m:
            continue

        if m.group("single") and not m.group("read"):
            sample = m.group("single").strip()
            samples.setdefault(sample, {"r1": None, "r2": None})
            samples[sample]["r1"] = str(f)
            continue

        sample = m.group("sample").strip()
        read = m.group("read")

        samples.setdefault(sample, {"r1": None, "r2": None})

        
        if read == "1":
            samples[sample]["r1"] = str(f)
        elif read == "2":
            samples[sample]["r2"] = str(f)

    return samples

#generate table
def generate_samples_table(fastq_dict=None, out_dir="", out_name="samples.csv"):
    if fastq_dict is None:
        fastq_dict = find_fastq_files()
    """Generate samples.csv"""
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    df = pd.DataFrame({
        "sample_name": sorted(fastq_dict.keys()),
        "guide_id": "",
        "group": "",
        "edit_mode": "",
        "project": "",
        "noise": "no"      })

    out_path = out_dir / out_name
    df.to_csv(out_path, index=False)
    print(f"Generated {out_path}")
    print("   Fill the table. For noise samples check noise=yes")
    return out_path

#load tables, make task_df
def load_tables():
    seq_df = pd.read_csv("sequences.csv", sep=None, engine="python",
                         dtype=str, keep_default_na=False)

    param_df = pd.read_csv("parameters.csv", sep=None, engine="python",
                           dtype=str, keep_default_na=False)

    samp_df = pd.read_csv("samples.csv", sep=None, engine="python",
                          dtype=str, keep_default_na=False)

    # --- базовый словарь (без additional) ---
    base_params = {
        row["parameters"].strip(): row["values"
            ].strip()
        for _, row in param_df.iterrows()
        if row["values"].strip().lower() != "additional"
    }

    return seq_df, samp_df, param_df, base_params
 
def expand_samples(samp_df, seq_df):
    tasks = []

    # dictionaries project → amplicon / expected
    amplicon_map = (
        seq_df.dropna(subset=["amplicon"])
              .drop_duplicates("project_id")
              .set_index("project_id")["amplicon"]
              .to_dict()
    )
        # checking if noise column exists
    has_noise_col = "noise" in samp_df.columns
    expected_map = (
        seq_df.dropna(subset=["expected_seq"])
              .drop_duplicates("project_id"
                              )
              .set_index("project_id")["expected_seq"]
              .to_dict()
    )

    for _, row in samp_df.iterrows():
        base = row.to_dict()
        project = row.get("project")

        base["amplicon"] = amplicon_map.get(project, "")
        base["expected_seq"] = expected_map.get(project, "")
    
        # определяем режим обработки 
        if has_noise_col:
            noise_flag = str(row.get("noise", "no")).lower()
        else:
            noise_flag = "yes"   

        base["noise"] = noise_flag

        if noise_flag == "yes":
            base["mode"] = "cas9"
        else:
            base["mode"] = "standard"

        # controls
        if int(row["group"]) == 0:
            # гайды берём ИЗ samp_df
            guide_ids = (
                samp_df[
                    (samp_df["project"] == project) &
                    # (samp_df["guide_id"].notna())
                    (samp_df["guide_id"]!="")
                ]["guide_id"]
                .unique()
                .tolist()
            )

            for gid in guide_ids:
                t = base.copy()
                t["guide_id"] = gid
                t["expanded_name"] = f"{row['sample_name']}__g{int(gid)}"
                tasks.append(t)

        # regular samples
        else:
            t = base.copy()
            t["expanded_name"] = row["sample_name"]
            tasks.append(t)

    return pd.DataFrame(tasks)
