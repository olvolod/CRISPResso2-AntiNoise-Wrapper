import pandas as pd
from pathlib import Path
import zipfile
import os
from Bio.Align import PairwiseAligner
from config import RESULTS_PATH
#samples before fallback (CRISPResso cmd)
def apply_project_overrides(cmd, project, param_df):

    for _, row in param_df.iterrows():
        param = row["parameters"].strip()
        value = row["values"].strip()
        project_list = row.get("project", "").strip()
        change = row.get("change", "").strip()

        if not project_list:
            continue

        projects = [p.strip() for p in project_list.split(",")]

        if project not in projects:
            continue

        # project index
        idx = projects.index(project)

        # additional 
        if value.lower() == "additional":
            change_val = row.get("change", "").strip()  
            if change_val:
                cmd += [param, change_val]
            else:
                cmd.append(param)
            continue

        # false
        if change.lower() == "false":
            if param in cmd:
                i = cmd.index(param)
                cmd = cmd[:i] + cmd[i+2:] if i+1 < len(cmd) else cmd[:i]
            continue

        # override 
        if change:
            changes = [c.strip() for c in change.split(",")]

            if idx < len(changes):
                new_value = changes[idx]

                # удалить старый
                if param in cmd:
                    i = cmd.index(param)
                    cmd = cmd[:i] + cmd[i+2:]

                # добавить новый
                cmd += [param, new_value]

    return cmd

def is_sample_done(sample_name, report_dir="report"):
    p = Path(report_dir) / f"CRISPResso_on_{sample_name}/Alleles_frequency_table.zip"
    return p.exists()

def build_command(task, seq_df, params, param_df, fastq_dict, report_dir="report"):

    sample_name = task["expanded_name"].strip()
    base_sample = task["sample_name"].strip()

    if base_sample not in fastq_dict:
        raise ValueError(f"No fastq found for sample {base_sample}")

    r1 = fastq_dict[base_sample].get("r1")
    r2 = fastq_dict[base_sample].get("r2")

    cmd = ["CRISPResso", "-r1", r1]
    if r2:
        cmd += ["-r2", r2]

    #  regime
    mode = task.get("edit_mode", "").strip().lower()

    if mode not in {"cas9", "base", "prime"}:
        raise ValueError(f"Unknown mode '{mode}' for sample {sample_name}")

    # forbidden params
    CAS9_FLAGS = {"-g"}
    BASE_FLAGS = {"--conversion_nuc_from", "--conversion_nuc_to", "--base_editor_output"}
    PRIME_FLAGS = {
        "--prime_editing_pegRNA_spacer_seq",
        "--prime_editing_pegRNA_extension_seq",
        "--prime_editing_pegRNA_scaffold_seq",
        "--prime_editing_pegRNA_extension_quantification_window_size",
        "--prime_editing_pegRNA_nicking_guide_seq",
        "--prime_editing_override_prime_edited_ref_seq",
        "--prime_editing_override_sequence_checks"
    }

    FORBIDDEN_PARAMS = {
        "cas9": BASE_FLAGS | PRIME_FLAGS,
        "base": PRIME_FLAGS,
        "prime": BASE_FLAGS | CAS9_FLAGS,
    }

    # Guide 
    guide_id = task.get("guide_id", "").strip()

    g = None
    if guide_id and guide_id in seq_df["guide_id"].values:
        g = seq_df.set_index("guide_id").loc[guide_id]

    if mode == "prime":
        if g is None:
            raise ValueError(f"No guide found for prime sample {sample_name}")

        peg_spacer = g.get("spacer", "").strip()
        peg_extension = g.get("extension", "").strip()
        peg_scaffold = g.get("scaffold", "").strip()

        if not peg_spacer or not peg_extension:
            print(f" Skipping {sample_name} — missing spacer or extension")
            return None

        cmd += ["--prime_editing_pegRNA_spacer_seq", peg_spacer]
        cmd += ["--prime_editing_pegRNA_extension_seq", peg_extension]

        if peg_scaffold:
            cmd += ["--prime_editing_pegRNA_scaffold_seq", peg_scaffold]

    else:
        # Cas9 / Base
        if g is not None:
            spacer = g.get("spacer", "").strip()
            if spacer:
                cmd += ["-g", spacer]

    # amplicon / expected 
    amplicon = task.get("amplicon", "").strip()
    expected_seq = task.get("expected_seq", "").strip()

    if amplicon != "":
        cmd += ["-a", amplicon]
    
    if expected_seq != "":
        cmd += ["-e", expected_seq]

    # name and output 
    cmd += ["-n", sample_name, "--output_folder", report_dir]

    # user parameters 
    for p, v in params.items():
        p = p.strip()
        v = v.strip()

        if not p:
            continue

        # not doubling
        if p in {"-a", "amplicon"} and "-a" in cmd:
            continue
        if p in {"-e", "expected_seq"} and "-e" in cmd:
            continue

        # forbidden for regime
        if p in FORBIDDEN_PARAMS.get(mode, set()):
            continue

        # flag no value
        if v == "":
            cmd.append(p)
        else:
            cmd += [p, v]
    project = task.get("project", "").strip()
    cmd = apply_project_overrides(cmd, project, param_df)
    print(f"\nDEBUG command for {sample_name}:")
    print(" ".join(cmd), "...")

    return cmd

def analyze_sample(sample_name):
    zip_path = f"report/CRISPResso_on_{sample_name}/Alleles_frequency_table.zip"
    sample_dir = f"report/CRISPResso_on_{sample_name}/"
    with zipfile.ZipFile(zip_path) as z:
        z.extract("Alleles_frequency_table.txt", "outcome/alleles")

    df = pd.read_csv("outcome/alleles/Alleles_frequency_table.txt", sep="\t")
    os.rename("outcome/alleles/Alleles_frequency_table.txt",
              f"outcome/alleles/{sample_name}_Alleles_frequency_table.txt")

    filtered_df = df[(df["Reference_Name"] == "Reference") & (df["Read_Status"] == "MODIFIED")]
    hdr_df = df[df["Reference_Name"] == "HDR"]
    hdr_perf = hdr_df[hdr_df["Read_Status"] == "UNMODIFIED"]
    hdr_imp = hdr_df[hdr_df["Read_Status"] == "MODIFIED"]

    deleted = filtered_df[filtered_df["n_deleted"] != 0]["#Reads"].sum()
    inserted = filtered_df[(filtered_df["n_deleted"] == 0) & (filtered_df["n_inserted"] != 0)]["#Reads"].sum()
    mut = filtered_df[(filtered_df["n_deleted"] == 0) & (filtered_df["n_inserted"] == 0)]["#Reads"].sum()

    total = df["#Reads"].sum()
    # Quantification window nucleotide frequency 

    possible_files = [
        "Reference.Quantification_window_nucleotide_frequency_table.txt",
        "Quantification_window_nucleotide_frequency_table.txt"
    ]

    qwin_file = None
    for fname in possible_files:
        candidate = os.path.join(sample_dir, fname)
        if os.path.exists(candidate):
            qwin_file = candidate
            break

    quant_window_seq = ""

    if qwin_file is not None:
        print(f" Using file: {os.path.basename(qwin_file)}")

        qwin_df = pd.read_csv(qwin_file, sep="\t", index_col=0)

        seq_list = []
        for col in qwin_df.columns:
            nuc = qwin_df[col].idxmax()
            seq_list.append(nuc)

        quant_window_seq = "".join(seq_list)
    else:
        print(f"⚠ Quantification window table not found for {sample_name}")
    return {
        "Sample": sample_name,
        "correction_%": float(f"{(hdr_df['#Reads'].sum() / total) * 100:.2f}"),
        "Perfect_correction_%": float(f"{(hdr_perf['#Reads'].sum() / total) * 100:.2f}"),
        "Imperfect_correction_%": float(f"{(hdr_imp['#Reads'].sum() / total) * 100:.2f}"),
        "undesired_%": float(f"{(filtered_df['#Reads'].sum() / total) * 100:.2f}"),
        "Indel_%": float(f"{((deleted + inserted) / total) * 100:.2f}"),
        "Mutation_%": float(f"{((filtered_df['#Reads'].sum() / total) * 100 - ((deleted + inserted) / total) * 100):.2f}"),
        "Reads": total,
        "QuantWindow": quant_window_seq
    }
def update_results_table(result):
    if RESULTS_PATH.exists():
        df = pd.read_csv(RESULTS_PATH)
        df = pd.concat([df, pd.DataFrame([result])], ignore_index=True)
    else:
        df = pd.DataFrame([result])

    df.to_csv(RESULTS_PATH, index=False)

#samples fallback
def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]
def get_aligned_strings(alignment, seqA, seqB): #seqA и seqB из таблицы, которую мы изначально добавляли
    alnA, alnB = [], []
    a_i = b_i = 0
    for (a_start, a_end), (b_start, b_end) in zip(alignment.aligned[0], alignment.aligned[1]):
        while a_i < a_start:
            alnA.append(seqA[a_i]); alnB.append("-"); a_i += 1
        while b_i < b_start:
            alnA.append("-"); alnB.append(seqB[b_i]); b_i += 1
        for _ in range(a_end - a_start):
            alnA.append(seqA[a_i]); alnB.append(seqB[b_i])
            a_i += 1; b_i += 1
    while a_i < len(seqA):
        alnA.append(seqA[a_i]); alnB.append("-"); a_i += 1
    while b_i < len(seqB):
        alnA.append("-"); alnB.append(seqB[b_i]); b_i += 1
    return "".join(alnA), "".join(alnB)

def align_expected_to_amplicon(amplicon, expected):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aln = aligner.align(amplicon, expected)[0]
    return get_aligned_strings(aln, amplicon, expected)


def project_extension_to_amplicon(amplicon, extension):

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    aln = aligner.align(amplicon, extension)[0]

    aligned_amp, aligned_ext = get_aligned_strings(aln, amplicon, extension)

    projected = []

    for a, b in zip(aligned_amp, aligned_ext):
        if b != "-":
            projected.append(a)

    return "".join(projected).replace("-", "")

def build_fallback_command(task, seq_df, params, param_df, fastq_dict, guide_seq, is_extension):
    
    #using existing build_command, but replacing the regime and guide.

  
    fallback_task = task.copy()


    fallback_task["edit_mode"] = "cas9"

    suffix = "extension" if is_extension else "spacer"
    fallback_task["expanded_name"] = f"{task['expanded_name']}__fallback_{suffix}"

    # replacing guide_id logically 
    # temporarily deliting guide_id to build command without it, then adding guide_seq manually
    fallback_task["guide_id"] = ""

    # if extension → deleting expected_seq
    # if is_extension:
        # fallback_task["expected_seq"] = ""

    # building command
    cmd = build_command(fallback_task, seq_df, params, param_df, fastq_dict)

    # deleting -g
    if "-g" in cmd:
        idx = cmd.index("-g")
        cmd = cmd[:idx] + cmd[idx+2:]
#    cmd.append("--suppress_plots")
    # adding guide
    insert_pos = cmd.index("-a")  # вставляем перед -a
    cmd = cmd[:insert_pos] + ["-g", guide_seq] + cmd[insert_pos:]
    # extra logic for extension: quantification window size 
    if is_extension:
        param_key = "--prime_editing_pegRNA_extension_quantification_window_size"
        if param_key in params and params[param_key].strip():
            w_value = params[param_key].strip()
            # deleting -w, if exists
            if "-w" in cmd:
                idx = cmd.index("-w")
                cmd = cmd[:idx] + cmd[idx+2:]
            # inserting new -w before -wc 
            if "-wc" in cmd:
                idx_wc = cmd.index("-wc")
                cmd = cmd[:idx_wc] + ["-w", w_value] + cmd[idx_wc:]
            else:
                cmd += ["-w", w_value]
    project = task.get("project", "").strip()
    cmd = apply_project_overrides(cmd, project, param_df)
    return cmd
# fallback analysis
def extract_variants_at_positions(aligned_read, aligned_ref, positions):
    positions = set(positions)
    variants = []

    ref_pos = 0  # 1-based по reference

    for r, a in zip(aligned_ref, aligned_read):
        if r != "-":
            ref_pos += 1

        if ref_pos not in positions:
            continue

        if r == a:
            continue

        # mutation
        if r != "-" and a != "-":
            variants.append({
                "change": f"{r}>{a}",
                "pos_type": f"mut{ref_pos}"
            })

        # deletion
        elif r != "-" and a == "-":
            variants.append({
                "change": f"del{r}",
                "pos_type": f"del{ref_pos}"
            })

        # insertion
        elif r == "-" and a != "-":
            variants.append({
                "change": f"ins{a}",
                "pos_type": f"ins{ref_pos}"
            })

    return variants
    
def annotate_window_variants(df, positions):
    df = df.copy()

    def process_row(row):
        variants = extract_variants_at_positions(
            row["Aligned_Sequence"],
            row["Reference_Sequence"],
            positions
        )

        if not variants:
            row["change"] = ""
            row["pos_type"] = ""
            return row

        row["change"] = ";".join([v["change"] for v in variants])
        row["pos_type"] = ";".join([v["pos_type"] for v in variants])
        return row

    df = df.apply(process_row, axis=1)
    return df
    
def get_quant_window(sample_name, mode):
    suffix = "__fallback_extension" if mode == "extension" else "__fallback_spacer"
    sample_dir = f"report/CRISPResso_on_{sample_name}{suffix}"

    # possible variants of file name
    possible_files = [
        "Reference.Quantification_window_nucleotide_frequency_table.txt",
        "Quantification_window_nucleotide_frequency_table.txt"
    ]

    qwin_file = None
    for fname in possible_files:
        candidate = os.path.join(sample_dir, fname)
        if os.path.exists(candidate):
            qwin_file = candidate
            break

    quant_window_seq = ""

    if qwin_file is not None:
        print(f"📂 Using file: {os.path.basename(qwin_file)}")

        qwin_df = pd.read_csv(qwin_file, sep="\t", index_col=0)

        seq_list = []
        for col in qwin_df.columns:
            nuc = qwin_df[col].idxmax()
            seq_list.append(nuc)

        quant_window_seq = "".join(seq_list)

        print(f"  {mode} quant window:")
        print(f"   sequence: {quant_window_seq}")
        print(f"   length:   {len(quant_window_seq)}")

    else:
        print(f" Quantification window table not found for {sample_name} ({mode})")

    return quant_window_seq

def get_windows_and_overlap(sample_name, task_df, seq_df):
    ext_qwin = get_quant_window(sample_name, "extension")
    spacer_qwin = get_quant_window(sample_name, "spacer")

    # chosing task
    row = task_df.query("expanded_name == @sample_name")
    if row.empty:
        raise ValueError(f"Sample {sample_name} not found in task_df")
    project_id = row.iloc[0]["project"]

    # chosing amplicon from seq_df
    seq_row = seq_df.query("project_id == @project_id")
    if seq_row.empty:
        raise ValueError(f"Project {project_id} not found in seq_df")
    amplicon = seq_row.iloc[0]["amplicon"].strip()

    ext_start = amplicon.find(ext_qwin)
    spacer_start = amplicon.find(spacer_qwin)

    if ext_start == -1 or spacer_start == -1:
        raise ValueError("Quant window not found in amplicon")

    ext_range = range(ext_start + 1, ext_start + 1 + len(ext_qwin))
    spacer_range = range(spacer_start + 1, spacer_start + 1 + len(spacer_qwin))

    overlap = set(ext_range) & set(spacer_range)

    print("Extension range:", ext_range.start, "-", ext_range.stop - 1)
    print("Spacer range:", spacer_range.start, "-", spacer_range.stop - 1)
    print("Overlap:", overlap)

    return overlap

def load_alleles(sample_name, mode):
    suffix = "__fallback_extension" if mode == "extension" else "__fallback_spacer"
    sample_dir = f"report/CRISPResso_on_{sample_name}{suffix}"
    zip_path = os.path.join(sample_dir, "Alleles_frequency_table.zip")

    with zipfile.ZipFile(zip_path) as z:
        with z.open("Alleles_frequency_table.txt") as f:
            df = pd.read_csv(f, sep="\t")

    print(f"\n {mode} alleles loaded:")
    print("rows:", len(df))
    print("sum reads:", df["#Reads"].sum())

    return df

def load_tables_fallback(sample_name):

    ext_df = load_alleles(sample_name, "extension")
    spacer_df = load_alleles(sample_name, "spacer")

    ext_total = ext_df["#Reads"].sum()
    spacer_total = spacer_df["#Reads"].sum()

    return ext_df, spacer_df, ext_total, spacer_total
import os
import pandas as pd


def update_spacer_from_extension(spacer_df, ext_df, sample_name):
    """
    Обновляет spacer по результатам extension.

    Логика:
    1. Берём только MODIFIED аллели из extension.
    2. Если последовательность есть в spacer:
       - переключаем Read_Status в MODIFIED;
       - переносим n_deleted/n_inserted/n_mutated.
    3. Последовательности, которых нет в spacer,
       сохраняем отдельным файлом.
    """

    spacer = spacer_df.copy()

    # словарь MODIFIED аллелей extension
    ext_lookup = {
        row["Aligned_Sequence"]: row
        for _, row in ext_df.iterrows()
        if row["Read_Status"] == "MODIFIED"
    }

    updated = 0

    for idx, row in spacer.iterrows():

        seq = row["Aligned_Sequence"]

        if seq not in ext_lookup:
            continue

        ext_row = ext_lookup[seq]

        # переключаем статус
        if row["Read_Status"] == "UNMODIFIED":
            spacer.at[idx, "Read_Status"] = "MODIFIED"
            updated += 1

        # переносим аннотацию из extension
        spacer.at[idx, "n_deleted"] = ext_row["n_deleted"]
        spacer.at[idx, "n_inserted"] = ext_row["n_inserted"]
        spacer.at[idx, "n_mutated"] = ext_row["n_mutated"]

        # последовательность обработана
        del ext_lookup[seq]

    print(f"Updated spacer statuses: {updated}")

    # ---------- сохраняем extension-only ----------
    extra = pd.DataFrame(ext_lookup.values())

    if not extra.empty:

        os.makedirs("outcome/extension_only_reads", exist_ok=True)

        save_path = (
            f"outcome/extension_only_reads/"
            f"{sample_name}_extension_only.tsv"
        )

        extra.to_csv(save_path, sep="\t", index=False)

        print(f"Saved {len(extra)} extension-only alleles")
        print(f"File: {save_path}")

        print("Extra reads:", extra["#Reads"].sum())

        wrong = extra[extra["n_deleted"] == 0]

        if not wrong.empty:
            print(
                f"[WARNING] {len(wrong)} extension-only alleles "
                f"without deletions (possible alignment artefacts)"
            )

    return spacer
# def update_spacer_from_extension(spacer_df, ext_df):
#     """
#     Обновляет Read_Status в spacer по результатам extension.

#     Логика:
#     1. Берём только MODIFIED аллели extension.
#     2. Если такая последовательность есть в spacer и там UNMODIFIED,
#        переключаем её в MODIFIED.
#     3. Если последовательности нет в spacer,
#        считаем её новой и добавляем в таблицу.
#     """

#     spacer = spacer_df.copy()

#     # словарь MODIFIED-последовательностей extension
#     ext_lookup = {
#         row["Aligned_Sequence"]: row
#         for _, row in ext_df.iterrows()
#         if row["Read_Status"] == "MODIFIED"
#     }

#     spacer_sequences = set(spacer["Aligned_Sequence"])

#     updated = 0

#     for idx, row in spacer.iterrows():

#         seq = row["Aligned_Sequence"]

#         if seq in ext_lookup:

#             if row["Read_Status"] == "UNMODIFIED":
#                 spacer.at[idx, "Read_Status"] = "MODIFIED"
#                 updated += 1

#             # последовательность уже обработали
#             del ext_lookup[seq]

#     print(f"Updated spacer statuses: {updated}")

#     # Остались только последовательности,
#     # которых вообще нет в spacer
#     extra = pd.DataFrame(ext_lookup.values())


#     if not extra.empty:

#         extra["Source"] = "Extension_only"

#         print(
#         f"Saved {len(extra)} extension-only alleles"
#         )

#         merged = pd.concat(
#         [spacer, extra],
#         ignore_index=True
#         )

#     else:
#         spacer = spacer
#     return spacer    
# def filter_extension_not_in_spacer(ext_df, spacer_df):

#     cols = ["Aligned_Sequence", "#Reads"]

#     check = ext_df.merge(
#         spacer_df[cols],
#         on=cols,
#         how="left",
#         indicator=True
#     )

#     extra = check[check["_merge"] == "left_only"]

#     print("Extra rows in extension:", len(extra))
#     print("Extra reads:", extra["#Reads"].sum())
    

#     if not extra.empty:
#         # Проверяем, есть ли строки НЕ являющиеся делециями
#         non_deletions = extra[extra["n_deleted"] == 0]
#         if not non_deletions.empty:
#             print(f"\n[WARNING] Найдено {len(non_deletions)} строк без делеций в extra!")
#             # Печатаем типы мутаций, которые вызывают проблему:
#             print(non_deletions[["Aligned_Sequence", "#Reads", "Read_Status", "n_deleted", "n_inserted"]].head(10))
#         else:
#             print("All mismatches are deletions → OK")
# #    if not extra.empty:
# #        if not (extra["n_deleted"] > 0).all():
# #            # display(extra)
# #            raise ValueError("Non-deletion rows found in extension mismatch")

# #        print("All mismatches are deletions → OK")

#     ext_filtered = check[check["_merge"] == "both"].drop(columns=["_merge"])

#     return ext_filtered

# def process_extension(ext_df, overlap):

#     ext_annot = annotate_window_variants(ext_df, overlap)

#     ext_no_overlap = ext_annot[
#         ext_annot["pos_type"] == ""
#     ].copy()

#     ext_modified = ext_no_overlap[
#         ext_no_overlap["Read_Status"] == "MODIFIED"
#     ].copy()

#     # удалить служебные колонки
#     for col in ["pos_type", "change"]:
#         if col in ext_modified.columns:
#             ext_modified.drop(columns=col, inplace=True)

#     print("Extension after annotation+filter:")
#     print("rows:", len(ext_modified))
#     print("reads:", ext_modified["#Reads"].sum())

#     return ext_modified

# def apply_hdr_logic(spacer_df, ext_df):

#     spacer_df = spacer_df.copy()
#     ext_df = ext_df.copy()

#     ext_sequences = set(ext_df["Aligned_Sequence"])

#     for idx, row in spacer_df.iterrows():
#         if row["Reference_Name"] == "HDR":
#             seq = row["Aligned_Sequence"]

#             if seq in ext_sequences:

#                 if row["Read_Status"] == "UNMODIFIED":
#                     spacer_df.loc[idx, "Read_Status"] = "MODIFIED"

#                 ext_df = ext_df[
#                     ext_df["Aligned_Sequence"] != seq
#                 ]

#     return spacer_df, ext_df

# def remove_extension_from_spacer(spacer_df, ext_df):

#     ext_sequences = set(ext_df["Aligned_Sequence"])

#     spacer_filtered = spacer_df[
#         ~spacer_df["Aligned_Sequence"].isin(ext_sequences)
#     ].copy()

#     return spacer_filtered

# def final_merge(spacer_df, ext_df, original_total):

#     merged = pd.concat([spacer_df, ext_df], ignore_index=True)

#     merged = merged.sort_values("#Reads", ascending=False)

#     final_reads = merged["#Reads"].sum()

#     if final_reads != original_total:
#         raise ValueError(
#             f"Read mismatch: expected {original_total}, got {final_reads}"
#         )

#     print("Final rows:", len(merged))
#     print("Final reads:", final_reads)

#     return merged

def process_all_samples(sample_list, task_df, seq_df):

    results = {}

    for sample in sample_list:
        try:
            merged = process_sample(sample, task_df, seq_df)
            results[sample] = merged
        except Exception as e:
            print(f"❌ Error in {sample}: {e}")

    return results

def is_fallback_sample_done(sample_name, report_dir="report"):
    p = Path(report_dir) / f"CRISPResso_on_{sample_name}__fallback_spacer/Alleles_frequency_table.zip"
    return p.exists()
def build_non_overlapping_quant_window(sample_name, task_df, seq_df):
    ext_q = get_quant_window(sample_name, "extension")
    spacer_q = get_quant_window(sample_name, "spacer")

    row = task_df.query("expanded_name == @sample_name")
    if row.empty:
        raise ValueError(f"Sample {sample_name} not found in task_df")
    project_id = row.iloc[0]["project"]

    seq_row = seq_df.query("project_id == @project_id")
    if seq_row.empty:
        raise ValueError(f"Project {project_id} not found in seq_df")
    amplicon = seq_row.iloc[0]["amplicon"]

    ext_start = amplicon.find(ext_q)
    spacer_start = amplicon.find(spacer_q)

    ext_positions = range(ext_start, ext_start + len(ext_q))
    spacer_positions = range(spacer_start, spacer_start + len(spacer_q))

    overlap = set(ext_positions) & set(spacer_positions)

    ext_filtered = "".join(
        base for i, base in enumerate(ext_q, start=ext_start)
        if i not in overlap
    )

    combined = spacer_q + ext_filtered

    return combined

def process_sample(sample_name, task_df, seq_df):

    print(f"\n=== Processing {sample_name} ===")

    overlap = get_windows_and_overlap(sample_name, task_df, seq_df)

    ext_df_raw, spacer_df_raw, ext_total, spacer_total = load_tables_fallback(sample_name)

    original_total = spacer_total 

    merged = update_spacer_from_extension(
        spacer_df_raw,
        ext_df_raw,
        sample_name
    )

    final_reads = merged["#Reads"].sum()

    print("Final rows:", len(merged))
    print("Final reads:", final_reads)

    if final_reads != original_total:
        print(
            f"[WARNING] Read count mismatch: "
            f"{final_reads} vs {original_total}"
        )

    return merged
    # ext_df_clean = filter_extension_not_in_spacer(
    #     ext_df_raw,
    #     spacer_df_raw
    # )

    # ext_df_processed = process_extension(ext_df_clean, overlap)

    # spacer_df_hdr, ext_df_hdr = apply_hdr_logic(
    #     spacer_df_raw,
    #     ext_df_processed
    # )

    # spacer_df_final = remove_extension_from_spacer(
    #     spacer_df_hdr,
    #     ext_df_hdr
    # )

    # merged = final_merge(
    #     spacer_df_final,
    #     ext_df_hdr,
    #     original_total
    # )

   # return merged

def analyze_fallback_sample(sample_name, task_df, seq_df):

    merged = process_sample(
        sample_name,
        task_df=task_df,
        seq_df=seq_df
    )

    os.makedirs("outcome", exist_ok=True)
    os.makedirs("outcome/alleles", exist_ok=True)
    save_path = f"outcome/alleles/{sample_name}_fallback_Alleles_frequency_table.txt"
    merged.to_csv(save_path, sep="\t", index=False)

    print(f"✅ Saved fallback table → {save_path}")

    df = merged

    filtered_df = df[
        (df["Reference_Name"] == "Reference") &
        (df["Read_Status"] == "MODIFIED")
    ]

    hdr_df = df[df["Reference_Name"] == "HDR"]
    hdr_perf = hdr_df[hdr_df["Read_Status"] == "UNMODIFIED"]
    hdr_imp = hdr_df[hdr_df["Read_Status"] == "MODIFIED"]

    deleted = filtered_df[filtered_df["n_deleted"] != 0]["#Reads"].sum()

    inserted = filtered_df[
        (filtered_df["n_deleted"] == 0) &
        (filtered_df["n_inserted"] != 0)
    ]["#Reads"].sum()

    total = df["#Reads"].sum()
    return {
        "Sample": f'{sample_name}_fallback',
        "correction_%": float(f"{(hdr_df['#Reads'].sum() / total) * 100:.2f}"),
        "Perfect_correction_%": float(f"{(hdr_perf['#Reads'].sum() / total) * 100:.2f}"),
        "Imperfect_correction_%": float(f"{(hdr_imp['#Reads'].sum() / total) * 100:.2f}"),
        "undesired_%": float(f"{(filtered_df['#Reads'].sum() / total) * 100:.2f}"),
        "Indel_%": float(f"{((deleted + inserted) / total) * 100:.2f}"),
        "Mutation_%": float(f"{((filtered_df['#Reads'].sum() / total) * 100 - ((deleted + inserted) / total) * 100):.2f}"),
        "Reads": total,
        "QuantWindow": build_non_overlapping_quant_window(
            sample_name,
            task_df,
            seq_df
        )
    }
