import pandas as pd
from pathlib import Path
import os
import re
from Bio.Align import PairwiseAligner
from scipy.stats import beta
from statsmodels.stats.proportion import proportion_confint
from config import ANNOTATION_DIR
import crispresso_functions as cf
# post-processing functions 
# annotation
def load_alleles_table(path): #для всех таблиц из папки с таблицами
    df = pd.read_csv(path, sep="\t")
    return df

def filter_alleles(df):
    df = df.copy()
    hdr_mask = df["Reference_Name"].isin(["HDR", "Prime-edited"]) & \
               (df["Read_Status"] == "MODIFIED") 
               # (df['%Reads'] > 0.01)
    ref_mask = ~hdr_mask
    mut_mask = (ref_mask & (df["n_mutated"] > 0))
    ins_mask = (
        ref_mask &
        (df["n_inserted"] > 0)
    )

    del_mask = (
        ref_mask &
        (df["n_deleted"] > 0)
    )

    return df[hdr_mask | mut_mask | ins_mask | del_mask].copy()

def build_hdr_amplicon_map(aln_amp, aln_exp):
    """
    Creates a mapping:
    hdr_pos (1-based) → amp_pos (1-based)

    Returns a dict
    """

    hdr_to_amp = {}

    amp_pos = 0
    hdr_pos = 0

    for a, e in zip(aln_amp, aln_exp):

        if a != "-":
            amp_pos += 1

        if e != "-":
            hdr_pos += 1

        # только если обе позиции реальные
        if a != "-" and e != "-":
            hdr_to_amp[hdr_pos] = amp_pos

    return hdr_to_amp

def extract_variants_with_mapping(
    aligned_read,
    aligned_ref,
    positions,
    hdr_to_amp_map=None
):
    """
    If hdr_to_amp_map is provided, the position recounts in amplicon coordinates.

    Returns:
    change
    pos_type
    """

    positions = set(positions)
    variants = []

    ref_pos = 0  # 1-based

    for r, a in zip(aligned_ref, aligned_read):

        if r != "-":
            ref_pos += 1

        if ref_pos not in positions:
            continue

        if r == a:
            continue

        # пересчет позиции
        final_pos = ref_pos
        if hdr_to_amp_map:
            final_pos = hdr_to_amp_map.get(ref_pos, None)
            if final_pos is None:
                continue

        # mutation
        if r != "-" and a != "-":
            variants.append({
                "change": f"{r}>{a}",
                "pos_type": f"mut{final_pos}"
            })

        # deletion
        elif r != "-" and a == "-":
            variants.append({
                "change": f"del{r}",
                "pos_type": f"del{final_pos}"
            })

        # insertion
        elif r == "-" and a != "-":
            variants.append({
                "change": f"ins{a}",
                "pos_type": f"ins{final_pos}"
            })

    return variants

def annotate_window_variants_stat(
    df,
    positions,
    hdr_map=None,
    expected_changes=None
):

    df = df.copy()

    def process_row(row):

        # HDR / Prime-edited
        if row["Reference_Name"] in ["HDR", "Prime-edited"] and hdr_map:
        
            variants = extract_variants_with_mapping(
                row["Aligned_Sequence"],
                row["Reference_Sequence"],
                positions,
                hdr_to_amp_map=hdr_map
            )
        
            variants = mark_expected_changes(
                variants,
                row,
                expected_changes,
                hdr_map
            )

        # Reference
        else:

            variants = extract_variants_with_mapping(
                row["Aligned_Sequence"],
                row["Reference_Sequence"],
                positions
            )

        if not variants:
            row["change"] = ""
            row["pos_type"] = ""
            return row

        row["change"] = ";".join(
            [v["change"] for v in variants]
        )

        row["pos_type"] = ";".join(
            [v["pos_type"] for v in variants]
        )

        return row

    return df.apply(process_row, axis=1)

def mark_expected_changes(
    variants,
    row,
    expected_changes,
    hdr_map
):

    if not expected_changes or not hdr_map:
        return variants

    # making several expected positions in amplicon coordinates
    expected_amp_positions = {}

    for ch in expected_changes:
        amp_pos = hdr_map.get(ch["pos"])
        if amp_pos:
            expected_amp_positions[amp_pos] = ch

    aligned_read = row["Aligned_Sequence"]
    aligned_ref = row["Reference_Sequence"]

    new_variants = []

    for v in variants:

        pos = int("".join(filter(str.isdigit, v["pos_type"])))

        if pos not in expected_amp_positions:
            new_variants.append(v)
            continue

        expected_ch = expected_amp_positions[pos]

        # checking in aligned
        ref_pos = 0

        for r, a in zip(aligned_ref, aligned_read):

            if r != "-":
                ref_pos += 1

            if ref_pos == pos:

                if expected_ch["type"] == "mut":
                    if a != expected_ch["alt"]:
                        v["pos_type"] = v["pos_type"] + "change"

                elif expected_ch["type"] == "del":
                    if a != "-":
                        v["pos_type"] = v["pos_type"] + "change"

                elif expected_ch["type"] == "ins":
                    if a != expected_ch["alt"]:
                        v["pos_type"] = v["pos_type"] + "change"

                break

        new_variants.append(v)

    return new_variants

def detect_expected_changes(aln_amp, aln_exp): 
    changes = []
    amp_pos = 0
    i = 0

    while i < len(aln_amp):
        a = aln_amp[i]
        e = aln_exp[i]

        # counting amplicon position
        if a != "-":
            amp_pos += 1

        # match
        if a == e:
            i += 1
            continue

        # INSERTION 
        if a == "-" and e != "-":
            ins_seq = []
            start_pos = amp_pos

            while i < len(aln_amp) and aln_amp[i] == "-" and aln_exp[i] != "-":
                ins_seq.append(aln_exp[i])
                i += 1

            changes.append({
                "pos": start_pos,
                "ref": "",
                "alt": "".join(ins_seq),
                "type": "ins"
            })
            continue

        # DELETION 
        if a != "-" and e == "-":
            changes.append({
                "pos": amp_pos,
                "ref": a,
                "alt": "",
                "type": "del"
            })
            i += 1
            continue

        # MUTATION 
        if a != "-" and e != "-":
            changes.append({
                "pos": amp_pos,
                "ref": a,
                "alt": e,
                "type": "mut"
            })
            i += 1
            continue

        i += 1

    return changes

def find_all_occurrences(seq, sub): 
    starts = []
    start = 0
    while True:
        idx = seq.find(sub, start)
        if idx == -1:
            break
        starts.append(idx)
        start = idx + 1
    return starts
def split_window_into_chunks(window, reference, min_len=6):
    chunks = []
    i = 0

    while i < len(window):
        found = False
        for j in range(len(window), i + min_len - 1, -1):
            sub = window[i:j]
            if sub in reference:
                chunks.append(sub)
                i = j
                found = True
                break
        if not found:
            i += 1

    return chunks
def normalize_window_chunks(chunks, window_len):
    if len(chunks) != 2:
        return chunks  # ничего не делаем

    c1, c2 = chunks

    if (len(c1) + len(c2)) != window_len:
        return chunks  # странный случай, не трогаем

    if len(c1) % 2 != 0:
        # перенос последнего nt из c1 в начало c2
        c2 = c1[-1] + c2
        c1 = c1[:-1]

    return [c1, c2]
def window_positions_in_reference(window, reference):
    chunks = split_window_into_chunks(window, reference)
    chunks = normalize_window_chunks(chunks, len(window))

    positions = []

    for chunk in chunks:
        start = reference.find(chunk)
        if start == -1:
            continue
        positions.extend(range(start + 1, start + 1 + len(chunk)))

    return sorted(positions)

def find_window_gaps(positions):
    """
    positions: list[int] - window positions (1-based)
    return: list[int] — all positions that are missing between windows
    """
    positions = sorted(positions)
    gaps = []

    for i in range(len(positions) - 1):
        cur_pos = positions[i]
        next_pos = positions[i + 1]

        if next_pos > cur_pos + 1:
            gaps.extend(range(cur_pos + 1, next_pos))

    return gaps

def annotate_gap_variants(df, gap_positions):
    """
    df: result of annotate_window_variants
    gap_positions: list of positions (1-based) that are missing between windows
    """
    df = df.copy()

    def process_row(row):
        # if already annotated, leave as is
        if row["change"] or row["pos_type"]:
            return row

        variants = cf.extract_variants_at_positions(
            row["Aligned_Sequence"],
            row["Reference_Sequence"],
            gap_positions
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
    
def compare_to_expected(variants, expected_changes):
    """
    Checking if the variant matches the expected change.
    If not — marks it as *change
    """

    expected_set = {
        (c["type"], c["pos"], c["alt"])
        for c in expected_changes
    }

    new_variants = []

    for v in variants:

        change_type = v["pos_type"][:3]
        pos = int(''.join(filter(str.isdigit, v["pos_type"])))
        alt = v["change"].split(">")[-1] if ">" in v["change"] else ""

        if (change_type, pos, alt) in expected_set:
            continue
        else:
            v["pos_type"] = v["pos_type"] + "2change"

        new_variants.append(v)

    return new_variants
    
def collapse_positions(pos_string):
    if pd.isna(pos_string):
        return pos_string

    parts = [p.strip() for p in str(pos_string).split(";") if p.strip()]
    if not parts:
        return pos_string

    result = []

    # DEL / INS
    for prefix in ["del", "ins"]:
        nums = []

        for p in parts:
            if p.startswith(prefix):
                num = re.findall(r"\d+", p)
                if num:
                    nums.append(int(num[0]))

        if not nums:
            continue

        nums = sorted(set(nums))

        # группируем подряд идущие
        ranges = []
        start = nums[0]
        prev = nums[0]

        for n in nums[1:]:
            if n == prev + 1:
                prev = n
            else:
                ranges.append((start, prev))
                start = n
                prev = n

        ranges.append((start, prev))

        for r in ranges:
            if r[0] == r[1]:
                result.append(f"{prefix}{r[0]}")
            else:
                result.append(f"{prefix}{r[0]}-{r[1]}")

    # MUT
    mut_nums = []
    for p in parts:
        if p.startswith("mut"):
            num = re.findall(r"\d+", p)
            if num:
                mut_nums.append(int(num[0]))

    mut_nums = sorted(set(mut_nums))

    for n in mut_nums:
        result.append(f"mut{n}")

    return ";".join(result)

def process_all_samples_annotation(task_df, seq_df):


    alleles_dir = Path("outcome/alleles")
    output_dir = Path("outcome/annotation")
    results_path = Path("outcome/results_noise.csv")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Загружаем results один раз
    
    results_df = pd.read_csv(results_path)

    for alleles_file in alleles_dir.glob("*_Alleles_frequency_table.txt"):

        
        # Names
        base_with_fallback = alleles_file.name.replace(
            "_Alleles_frequency_table.txt", ""
        )

        base_without_fallback = base_with_fallback.replace("_fallback", "")

        print(f"\n=== Processing {base_with_fallback} ===")

        # Retrieve project_id
        task_row = task_df[
            task_df["expanded_name"] == base_without_fallback
        ]

        if task_row.empty:
            print("⚠ No project found in task_df")
            continue

        project_id = task_row.iloc[0]["project"]

        # Retrieve amplicon + expected
        seq_row = seq_df[
            seq_df["project_id"] == project_id
        ]

        if seq_row.empty:
            print("⚠ No sequences found")
            continue

        amplicon = seq_row.iloc[0]["amplicon"]
        expected = seq_row.iloc[0]["expected_seq"]

        # 4. Retrieve window from results.csv
        res_row = results_df[
            results_df["Sample"] == base_with_fallback
        ]

        if res_row.empty:
            print("⚠ No window in results.csv")
            continue

        window_seq = res_row.iloc[0]["QuantWindow"]

      
        # 5. Loading alleles
        df = load_alleles_table(alleles_file)
        hdr_mask = df["Reference_Name"].isin(["HDR", "Prime-edited"]) & \
                   (df["Read_Status"] == "MODIFIED") 
                   # (df['%Reads'] > 0.01)
        
        ref_mask = ~hdr_mask
        df = df[(~ref_mask) | ((df["n_deleted"] <= 5) & (df["n_inserted"] <= 5))]
        df = filter_alleles(df)

        # 6. Preparing HDR mapping
     
        hdr_map = None
        expected_changes = []
        
        if isinstance(expected, str) and expected.strip() != "":
        
            print("→ HDR detected, building alignment map")
        
            aln_amp, aln_exp = cf.align_expected_to_amplicon(
                amplicon,
                expected
            )
        
            expected_changes = detect_expected_changes(
                aln_amp,
                aln_exp
            )
        
            hdr_map = build_hdr_amplicon_map(
                aln_amp,
                aln_exp
            )
        
        else:
            print("→ No expected_seq, skipping HDR logic")

        # 7. Determining window positions 
        window_positions = window_positions_in_reference(
            window_seq,
            amplicon
        )

        if not window_positions:
            print("⚠ Window does not map")
            continue

        # 8. Annotating variants
        df = annotate_window_variants_stat(
            df=df,
            positions=window_positions,
            hdr_map=hdr_map,
            expected_changes=expected_changes
        )

        # 9. Annotating gaps within the window
        gap_positions = find_window_gaps(window_positions)

        if gap_positions:
            df = annotate_gap_variants(
                df,
                gap_positions
            )
            
        df["pos_type"] = df["pos_type"].apply(collapse_positions)

        # 10. Sorting

        df = df.sort_values(
            by="%Reads",
            ascending=False
        )

        # 11. Saving
        out_path = output_dir / f"{base_with_fallback}_positions.txt"

        df.to_csv(
            out_path,
            sep="\t",
            index=False
        )

        print(f"✔ Saved {out_path}")

#логика reference
def load_annotation_file(path):
    """
    Downloads annotation file and returns:
    - df
    - sample_name_for_output
    - sample_name_for_task
    """

    print(f"\n Loading file: {path.name}")

    # имя без _positions.txt
    file_prefix = path.stem.replace("_positions", "")

    # имя для сохранения (убираем _fallback если есть)
    sample_name_for_task = file_prefix.replace("_fallback", "")

    # имя для output (оставляем fallback если есть)
    sample_name_for_output = file_prefix

    df = pd.read_csv(path, sep="\t")

    return df, sample_name_for_output, sample_name_for_task

def build_control_position_table(annotation_dir, task_df, results_df):

    rows = []

    controls = task_df[task_df["group"] == "0"]
    all_projects = task_df["project"].unique()
    projects_with_controls = controls["project"].unique()
    projects_without_controls = set(all_projects) - set(projects_with_controls)

    for _, row in controls.iterrows():

        control_name = row["expanded_name"]
        project = row["project"]
        guide = row["guide_id"]

        # depth
        depth_match = results_df[
            results_df["Sample"].str.replace("_fallback", "", regex=False)
            == control_name
        ]

        if depth_match.empty:
            continue

        depth = depth_match["Reads"].values[0]

        # positions file 
        matches = list(Path(annotation_dir).glob(f"{control_name}*_positions.txt"))
        if not matches:
            continue

        df = pd.read_csv(matches[0], sep="\t")

        grouped = (
            df.groupby("pos_type")
              .agg({
                  "#Reads": "sum",
                  "%Reads": "sum"
              })
              .reset_index()
        )

        for _, g in grouped.iterrows():
            rows.append({
                "control": control_name,
                "pos": g["pos_type"],
                "#reads": g["#Reads"],
                "%reads": g["%Reads"],
                "depth": depth,
                "project": project,
                "guide": guide
            })

    table = pd.DataFrame(rows)

  
    # explode and summary

    table["pos"] = table["pos"].str.split(";")
    table = table.explode("pos")

    # reconstruct
    table = (
        table.groupby(
            ["control", "pos", "depth", "project", "guide"],
            as_index=False
        )
        .agg({
            "#reads": "sum",
            "%reads": "sum"
        })
    )

    return table, projects_without_controls
def build_noise_model(control_table, ci_level=0.95):

    control_table = control_table.copy()

    control_table["#reads"] = pd.to_numeric(
        control_table["#reads"], errors="coerce"
    )

    control_table["depth"] = pd.to_numeric(
        control_table["depth"], errors="coerce"
    )

    control_table = control_table.dropna(subset=["#reads", "depth"])

    #  Global statistics (mut/del/ins)


    control_table["pos_type"] = control_table["pos"].str[:3]

    type_stats = {}

    for t in ["mut", "del", "ins"]:

        df_t = control_table[control_table["pos_type"] == t]

        if df_t.empty:
            continue

        pct = df_t["#reads"] / df_t["depth"]

        type_stats[t] = {
            "median": pct.median(),
            "weight": df_t["pos"].nunique()
        }

    thresholds = []

    grouped = control_table.groupby(["project", "guide", "pos"])

    for (project, guide, pos), df_pos in grouped:

        n_controls = df_pos["control"].nunique()

        pos_type = pos[:3]
        k_total = df_pos["#reads"].sum()
        n_total = df_pos["depth"].sum()

        #  n = 1 or 2 → Bayesian

        if n_controls <= 2:

            # k_total = df_pos["#reads"].sum()
            # n_total = df_pos["depth"].sum()

            if n_total == 0:
                threshold = 0,
                noise_mean = 0
            else:
                global_p = type_stats[pos_type]["median"]
                prior_weight = type_stats[pos_type]["weight"]

                alpha0 = global_p * prior_weight
                beta0 = (1 - global_p) * prior_weight

                alpha_post = alpha0 + k_total
                beta_post = beta0 + (n_total - k_total)
                # posterior mean
                noise_mean = (alpha_post / (alpha_post + beta_post))

                threshold = beta.ppf(
                    ci_level,
                    alpha_post,
                    beta_post
                ) * 100

        #  n ≥ 3 → Binomial CI (Clopper-Pearson)

        else:

            noise_mean = k_total / n_total
            if n_total == 0:
                threshold = 0,
                noise_mean = 0
            else:
                ci_low, ci_up = proportion_confint(
                    k_total,
                    n_total,
                    alpha=1-ci_level,
                    method="beta"
                )

                threshold = ci_up * 100
                
        thresholds.append({
            "project": project,
            "guide": guide,
            "pos": pos,
            "control_mean": noise_mean * 100,
            "threshold": threshold,
            "n_controls": n_controls
        })

    return pd.DataFrame(thresholds)
def check_samples_vs_threshold(task_df, noise_df):
    """
    For each sample:
    - explode pos_type
    - calculate sum of %Reads by position
    - pull control mean and threshold
    - compare with threshold
    """

    print("\n🟢 Checking samples vs threshold")

    results = []

    samples = task_df[task_df["group"] != "0"]

    print("Total samples:", len(samples))

    for _, sample_row in samples.iterrows():

        sample_name = sample_row["expanded_name"]
        project = sample_row["project"]
        guide = sample_row["guide_id"]
        group = sample_row["group"]

        print(f"\n➡ Sample: {sample_name}")
        print(f"   Project: {project} | Guide: {guide}")

        noise_subset = noise_df[
            (noise_df["project"] == project) &
            (noise_df["guide"] == guide)
        ]

        if noise_subset.empty:
            print("    No noise model found, skipping")
            continue

        matches = list(ANNOTATION_DIR.glob(f"{sample_name}*_positions.txt"))
        if not matches:
            print("   Positions file not found")
            continue

        df, _, _ = load_annotation_file(matches[0])

        # remove HDR / Prime
        df = df[~df["Reference_Name"].isin(["HDR", "Prime-edited"])].copy()
        #explode
        df["pos"] = df["pos_type"].str.split(";")
        df = df.explode("pos")
        df["pos"] = df["pos"].str.strip()
        
        # position sum
        pos_sums = (
            # df.groupby("pos_type")["%Reads"]
            df.groupby("pos")["%Reads"]
            .sum()
            .reset_index()
            .rename(columns={"%Reads": "sample_sum"})
        )

        # merging with noise (adding mean control and threshold)
        merged = pos_sums.merge(
            noise_subset[["pos", "control_mean", "threshold"]], #noise_subset[["pos", "sum", "threshold"]],
            # left_on="pos_type",
            # right_
            on="pos",
            how="left"
        )

        # merged = merged.rename(columns={
        #     "sum": "control_mean"
        # })

        # if pos not in control
        merged["control_mean"] = merged["control_mean"].fillna(0)
        merged["threshold"] = merged["threshold"].fillna(0)

        # comparison
        merged["over_threshold"] = merged["sample_sum"] > merged["threshold"]

        print("   Positions over threshold:",
              merged["over_threshold"].sum())

        # constructing result
        for _, row in merged.iterrows():
            results.append({
                "sample": sample_name,
                "guide": guide,
                "group": group,
                "project": project,
                "pos": row["pos"],
                "sample_sum": row["sample_sum"],
                "control_mean": row["control_mean"],
                "threshold": row["threshold"],
                # "control_mean": row ["control_mean"],
                "over_threshold": "yes" if row["over_threshold"] else "no"
            })

    final_df = pd.DataFrame(results)
    print("\n Done checking samples")
    print(final_df.head())

    return final_df
def filter_by_guide_frequency(variants_df, meta_df):
    """
    variants_df — первая таблица (с pos)
    meta_df — таблица с guide_id
    """

    # merging
    df = variants_df.merge(
        meta_df[["sample_name", "guide_id"]],
        left_on="sample",
        right_on="sample_name",
        how="left"
    )

    # counting samples per guide
    guide_sample_counts = (
        df.groupby("guide_id")["sample"]
        .nunique()
        .reset_index(name="n_samples")
    )

    df = df.merge(guide_sample_counts, on="guide_id", how="left")

    # counting how many samples have each pos within a guide
    pos_counts = (
        df.groupby(["guide_id", "pos"])["sample"]
        .nunique()
        .reset_index(name="pos_sample_count")
    )

    df = df.merge(pos_counts, on=["guide_id", "pos"], how="left")

    # 4 applying filtering rule
    filtered_df = df[
        (df["n_samples"] == 1) |        # if sample  → keep
        (df["pos_sample_count"] >= 2)   # if >=2 samples have the variant
    ].copy()

    # 5 deleting helper columns
    filtered_df = filtered_df.drop(
        columns=["sample_name", "n_samples", "pos_sample_count"]
    )

    return filtered_df

def sum_edit_types_per_sample(df):
    df = df.copy()

    # 1 определяем тип
    df["edit_type"] = df["pos"].str.extract(r"^(del|ins|mut)")

    # 2 суммируем
    result = (
        df.groupby(["sample", "edit_type"])["adjusted_sum"]
        .sum()
        .unstack(fill_value=0)
        .reset_index()
    )

    # гарантируем наличие всех колонок
    for col in ["del", "ins", "mut"]:
        if col not in result.columns:
            result[col] = 0

    result = result[["sample", "del", "ins", "mut"]]

    return result
def clean_multi_positions(df):
    """
    df — таблица после threshold фильтра
    """

    df = df.copy()

    # working on each sample separately
    samples = df["sample"].unique()

    cleaned_frames = []
    blacklist_dict = {}   # ← here goes no_list
    
    for sample in samples:

        sample_df = df[df["sample"] == sample].copy()

        # blacklist (no_list) 
        no_rows = sample_df[sample_df["over_threshold"] == "no"]

        no_list = set()

        for pos in no_rows["pos"]:
            parts = str(pos).split(";")
            for p in parts:
                if p.strip():
                    no_list.add(p.strip())


        # saving blacklist for sample
        blacklist_dict[sample] = no_list

        #  processing yes strings 
        yes_rows = sample_df[sample_df["over_threshold"] == "yes"].copy()

        updated_rows = []

        for _, row in yes_rows.iterrows():

            pos = str(row["pos"])

            # if no ";", do not touch
            if ";" not in pos:
                updated_rows.append(row)
                continue

            parts = [p.strip() for p in pos.split(";") if p.strip()]

            # deleting elements that are in no_list
            filtered = [p for p in parts if p not in no_list]

            # if after filtering nothing is left → skip
            if len(filtered) == 0:
                continue

            # if one element is left → sum to it
            if len(filtered) == 1:
                row["pos"] = filtered[0]
                updated_rows.append(row)
                continue

            # if several elements are left → join them back
            row["pos"] = ";".join(filtered)
            updated_rows.append(row)

        cleaned_sample_df = pd.DataFrame(updated_rows)

        cleaned_frames.append(cleaned_sample_df)

    final_df = pd.concat(cleaned_frames, ignore_index=True)

    return final_df, blacklist_dict

#HDR logic imperfect
def check_HDR_vs_threshold(task_df, noise_df):
    """
    for each sample:
    - explode pos_type
    - count the sum of %Reads by position
    - pull in the control mean and threshold
    - compare with the threshold
    """

    print("\nChecking samples vs threshold")

    results = []
    no_noise_samples = []
    samples = task_df[task_df["group"] != "0"]

    print("Total samples:", len(samples))

    for _, sample_row in samples.iterrows():

        sample_name = sample_row["expanded_name"]
        project = sample_row["project"]
        guide = sample_row["guide_id"]
        group = sample_row["group"]

        print(f"\n➡ Sample: {sample_name}")
        print(f"   Project: {project} | Guide: {guide}")

        noise_subset = noise_df[
            (noise_df["project"] == project) &
            (noise_df["guide"] == guide)
        ]

        if noise_subset.empty:
            print("No noise model found, skipping")
            continue

        matches = list(ANNOTATION_DIR.glob(f"{sample_name}*_positions.txt"))
        if not matches:
            print("Positions file not found")
            continue

        df, _, _ = load_annotation_file(matches[0])

        # delete HDR / Prime
        df = df[~df["Reference_Name"].isin(["Reference"])].copy()

        # sum per position
        pos_sums = (
            df.groupby("pos_type")["%Reads"]
            .sum()
            .reset_index()
            .rename(columns={"%Reads": "sample_sum"})
        )

        # merge with noise (adding control mean and threshold)
        merged = pos_sums.merge(
            noise_subset[["pos", "threshold", "control_mean"]],
            left_on="pos_type",
            right_on="pos",
            how="left"
        )



        # if position not in control
        merged["threshold"] = merged["threshold"].fillna(0)
        merged["control_mean"] = merged["control_mean"].fillna(0)

        # comparison
        merged["over_threshold"] = merged["sample_sum"] > merged["threshold"]

        print("   Positions over threshold:",
              merged["over_threshold"].sum())

        # collecting results
        for _, row in merged.iterrows():
            results.append({
                "sample": sample_name,
                "guide": guide,
                "group": group,
                "project": project,
                "pos": row["pos_type"],
                "sample_sum": row["sample_sum"],
                "threshold": row["threshold"],
                "control_mean": row["control_mean"],
                "over_threshold": "yes" if row["over_threshold"] else "no"
            })

    final_df = pd.DataFrame(results)
    print("\n Done checking samples")
    print(final_df.head())

    return final_df, no_noise_samples

def clean_hdr_positions_with_noise(df, blacklist_dict, noise_df, debug=False):

    df = df.copy()
    cleaned_rows = []

    # to faster find threshold
    noise_lookup = (
        noise_df
        .set_index(["project", "guide", "pos", "control_mean"])["threshold"]
        .to_dict()
    )

    for _, row in df.iterrows():

        sample = row["sample"]
        guide = row["guide"]
        project = row["project"]
        sample_sum = row["sample_sum"]

        blacklist = blacklist_dict.get(sample, set())

        parts = [p.strip() for p in str(row["pos"]).split(";") if p.strip()]
        kept_positions = []

        if debug:
            print("\n")
            print("Sample:", sample)
            print("Original:", parts)

        for p in parts:

            base_p = p.replace("change", "")

            # if  change → saving
            if "change" in p:
                kept_positions.append(p)
                if debug:
                    print(p, "→ KEEP (change)")
                continue

            # finding threshold в noise_df
            noise_threshold = noise_lookup.get((project, guide, base_p))

            if debug:
                print("Checking:", p,
                      "| in blacklist:", base_p in blacklist,
                      "| noise_threshold:", noise_threshold)

            # if position not in noise_df — can be kept 
            if noise_threshold is None:
                kept_positions.append(p)
                continue

            if base_p in blacklist:
                if sample_sum > noise_threshold:
                    kept_positions.append(p)
                    if debug:
                        print("→ KEEP (blacklist but > noise)")
                else:
                    if debug:
                        print("→ DROP (blacklist and <= noise)")
            else:
                kept_positions.append(p)
                if debug:
                    print("→ KEEP (not in blacklist)")

        if kept_positions:
            row["pos"] = ";".join(kept_positions)
            cleaned_rows.append(row)
        else:
            if debug:
                print("Row removed entirely")

    if not cleaned_rows:
        return pd.DataFrame(columns=df.columns)

    return pd.DataFrame(cleaned_rows).reset_index(drop=True)

def adjust_corrections(edit_df, results_df):

    # 1. counting imperf_fin
    edit_df = edit_df.copy()
    edit_df["imperf_fin"] = edit_df[["del", "ins", "mut"]].sum(axis=1)

    imperf_lookup = (
        edit_df
        .set_index("sample")["imperf_fin"]
        .to_dict()
    )

    results_df = results_df.copy()

    adjusted_rows = []

    for _, row in results_df.iterrows():

        sample_full = row["Sample"]
        sample = sample_full.replace("_fallback", "")

        #  imperf_fin (or 0)
        imperf_fin = imperf_lookup.get(sample, 0)

        #  % → в доли
        perfect = row["Perfect_correction_%"] 
        imperfect = row["Imperfect_correction_%"] 

        # counting a
        a = float(imperfect) - float(imperf_fin)

        perfect_adj = float(perfect) + float(a)
        imperfect_adj = imperf_fin

        adjusted_rows.append({
            "sample": sample,
            "perfect_correction_adj": perfect_adj,
            "imperfect_correction_adj": imperfect_adj
        })

    return pd.DataFrame(adjusted_rows)

#final summary
def build_final_summary(
    final_adjusted,
    results_df,
    ref_adjusted,
    task_df,
    info_projects_without_noise
):

    # preparation       
    results_df = results_df.copy()
    task_df = task_df.copy()
    final_adjusted = final_adjusted.copy()
    ref_adjusted = ref_adjusted.copy()

    # deleting _fallback
    results_df["sample_clean"] = results_df["Sample"].str.replace("_fallback", "", regex=False)

    # merge results + task
    merged = results_df.merge(
        task_df[["expanded_name", "project", "guide_id", "group"]],
        left_on="sample_clean",
        right_on="expanded_name",
        how="left"
    )

    # merge final_adjusted
    merged = merged.merge(
        final_adjusted,
        left_on="sample_clean",
        right_on="sample",
        how="left"
    )

    # merge ref_adjusted
    ref_adjusted["undesired_calc"] = (
        ref_adjusted[["del", "ins", "mut"]].sum(axis=1)
    )

    merged = merged.merge(
        ref_adjusted[["sample", "del", "ins", "mut", "undesired_calc"]],
        left_on="sample_clean",
        right_on="sample",
        how="left"
    )

    merged.fillna(0, inplace=True)
    print(merged)
    output_rows = []
    numeric_cols = [
        "correction_%",
        "Perfect_correction_%",
        "Imperfect_correction_%",
        "undesired_%",
        "Indel_%",
        "Mutation_%",
        "perfect_correction_adj",
        "imperfect_correction_adj",
        "del",
        "ins",
        "mut",
        "undesired_calc"
    ]
    
    for col in numeric_cols:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)
            
    for _, row in merged.iterrows():

        sample = row["sample_clean"]
        project = str(row["project"])
        group = str(row["group"])
        guide = row["guide_id"]

        is_control = group == "0"
        no_noise = project in info_projects_without_noise

        #  CONTROLS 
        if is_control:
            output_rows.append({
                "sample": sample,
                "correction": round(row["correction_%"], 2),
                "correction_adj": "N/A",
                "perfect": round(row["Perfect_correction_%"], 2),
                "perfect_adj": "N/A",
                "imperfect": round(row["Imperfect_correction_%"], 2),
                "imperfect_adj": "N/A",
                "undesired": "N/A",
                "indel": "N/A",
                "mutation": "N/A",
                "noise_model": "N/A",
                "project": project,
                "guide": guide,
                "group": group
            })
            continue

        # NO NOISE MODEL 
        if no_noise:
            output_rows.append({
                "sample": sample,
                "correction": round(row["correction_%"], 2),
                "correction_adj": "N/A",
                "perfect": round(row["Perfect_correction_%"], 2),
                "perfect_adj": "N/A",
                "imperfect": round(row["Imperfect_correction_%"], 2),
                "imperfect_adj": "N/A",
                "undesired": round(row["undesired_%"], 2),
                "indel": round(row["Indel_%"], 2),
                "mutation": round(row["Mutation_%"], 2),
                "noise_model": "no",
                "project": project,
                "guide": guide,
                "group": group
            })
            continue

        # NORMAL SAMPLE
        perfect_adj = row["perfect_correction_adj"]
        imperfect_adj = row["imperfect_correction_adj"]

        correction = perfect_adj + imperfect_adj

        output_rows.append({
            "sample": sample,
            "correction": round(correction, 2),
            "correction_adj": round(correction, 2),
            "perfect": round(perfect_adj, 2),
            "perfect_adj": round(perfect_adj, 2),
            "imperfect": round(imperfect_adj, 2),
            "imperfect_adj": round(imperfect_adj, 2),
            "undesired": round(row["undesired_calc"], 2),
            "indel": round(row["del"] + row["ins"], 2),
            "mutation": round(row["mut"], 2),
            "noise_model": "yes",
            "project": project,
            "guide": guide,
            "group": group
        })

    final_df = pd.DataFrame(output_rows)

    return final_df

#логика base editing
def load_annotation(sample, annotation_dir, debug=False):
    """
    Loads annotation file for sample (expanded_name)

    sample: expanded_name
    annotation_dir: folder with *_positions.txt
    """

    filename = f"{sample}_positions.txt"
    filepath = os.path.join(annotation_dir, filename)

    if debug:
        print(f"Loading annotation: {filepath}")

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Annotation file not found: {filepath}")

    df = pd.read_csv(filepath, sep="\t")

    if debug:
        print(f"Loaded {len(df)} rows")

    return df

def collect_control_stats(df, sample, pos, nuc_from, nuc_to, depth, debug=False):

    df_pos = df[df["pos_type"] == f"mut{pos}"]

    mut_reads = df_pos["#Reads"].sum()

    target_change = f"{nuc_from}>{nuc_to}"

    specific_reads = df_pos[df_pos["change"] == target_change]["#Reads"].sum()

    if debug:
        print(f"\nSample: {sample}")
        print(f"Position: {pos}")
        print(f"mut_reads: {mut_reads}")
        print(f"specific_reads: {specific_reads}")
        print(f"depth: {depth}")

    return {
        "ctrl": sample,
        "pos": pos,
        "mut": mut_reads,
        "specific": specific_reads,
        "depth": depth
    }

def build_control_table(
        task_df,
        results_df,
        params,
        load_annotation,
        debug=False):

    nuc_from = params["--conversion_nuc_from"]
    nuc_to = params["--conversion_nuc_to"]

    if debug:
        print(f"\nConversion: {nuc_from}>{nuc_to}")

    # Only base controls
    controls = task_df[
        (task_df["group"] == "0") &
        (task_df["edit_mode"] == "base")
    ]

    rows = []

    for _, ctrl_row in controls.iterrows():

        sample = ctrl_row["expanded_name"]
        project = ctrl_row["project"]
        guide = ctrl_row["guide_id"]
        amplicon = ctrl_row["amplicon"]

        if debug:
            print(f"\nProcessing control: {sample}")

        # depth + window из results
        depth_row = results_df[results_df["Sample"] == sample]

        if depth_row.empty:
            if debug:
                print(f"No depth info for {sample}")
            continue

        depth = float(depth_row["Reads"].iloc[0])
        window = depth_row["QuantWindow"].iloc[0]

        # Automatically find editable positions in the window
        positions = find_editable_positions(
            amplicon=amplicon,
            window=window,
            nuc_from=nuc_from
        )

        if debug:
            print("Editable positions:", positions)

        # annotation loading (only once)
        df_annotation = load_annotation(
            sample,
            "outcome/annotation",
            debug
        )

        # counting stats for each position
        for pos in positions:

            stats = collect_control_stats(
                df_annotation,
                sample,
                pos,
                nuc_from,
                nuc_to,
                depth,
                debug
            )

            row = {
                "ctrl": sample,
                "project": project,
                "guide": guide,
                "pos": pos,
                "mut": stats["mut"],
                "specific": stats["specific"],
                "depth": stats["depth"],
                "spec_mut_type": f"{nuc_from}>{nuc_to}",
                "QuantWindow": window
            }

            rows.append(row)

    return pd.DataFrame(rows)

def build_base_noise_model(control_table, debug=False):

    if control_table.empty:
        if debug:
            print("Control table is empty")
        return pd.DataFrame()

    # to numeric
    for col in ["mut", "specific", "depth"]:
        control_table[col] = pd.to_numeric(control_table[col], errors="coerce")

    grouped = (
        control_table
        .groupby([
            "project",
            "guide",
            "QuantWindow",
            "pos",
            "spec_mut_type"
        ])
        .agg({
            "mut": "sum",
            "specific": "sum",
            "depth": "sum"
        })
        .reset_index()
    )

    # counting percentages
    def calc_row(row):

        depth = row["depth"]

        if depth == 0:
            return pd.Series({
                "percent_mut": 0,
                "percent_specific": 0,
                "percent_specific_from_mut": 0
            })

        percent_mut = (row["mut"] / depth) * 100
        percent_specific = (row["specific"] / depth) * 100

        if percent_mut == 0:
            percent_specific_from_mut = 0
        else:
            percent_specific_from_mut = (percent_specific / percent_mut) * 100

        return pd.Series({
            "percent_mut": percent_mut,
            "percent_specific": percent_specific,
            "percent_specific_from_mut": percent_specific_from_mut
        })

    stats = grouped.join(grouped.apply(calc_row, axis=1))

    if debug:
        print("\nNoise model built:")
        print(stats.head())

    return stats
    
def find_editable_positions(amplicon, window, nuc_from):
    """
    Automatically finds positions of nuc_from in the window and returns coordinates
    relative to the amplicon.
    """

    start = amplicon.find(window)

    if start == -1:
        raise ValueError("Window not found in amplicon")

    positions = []

    for i, base in enumerate(window):
        if base == nuc_from:
            positions.append(start + i + 1)  # 1-based

    return positions

def get_pos_stats(annotation_df, pos, nuc_from, nuc_to):

    pos_label = f"mut{pos}"
    annotation_df["%Reads"] = pd.to_numeric(annotation_df["%Reads"], errors="coerce")
    df_pos = annotation_df[annotation_df["pos_type"] == pos_label]
    
    # all mutations at this position
    total_mut = df_pos["%Reads"].sum()

    # needed mutation at this position
    specific_mut = df_pos[
        df_pos["change"] == f"{nuc_from}>{nuc_to}"
    ]["%Reads"].sum()

    return total_mut, specific_mut

def compute_control_percentages(ctrl_df):

    grouped = ctrl_df.groupby("pos").sum()

    grouped["percent_mut"] = grouped["mut"] / grouped["depth"]
    grouped["percent_specific"] = grouped["specific"] / grouped["depth"]

    grouped["percent_specific_from_mut"] = (
        grouped["specific"] / grouped["mut"]
    ) * 100

    grouped = grouped.fillna(0)

    return grouped
def build_sample_base_table(
    task_df,
    results_df,
    params,
    load_annotation,
    debug=False
):

    nuc_from = params["--conversion_nuc_from"]
    nuc_to = params["--conversion_nuc_to"]

    rows = []

    # not controls, only base editing
    samples = task_df[
        (task_df["edit_mode"] == "base") &
        (task_df["group"] != "0")
    ]

    if debug:
        print(f"Samples to process: {len(samples)}")

    for _, sample_row in samples.iterrows():

        sample = sample_row["expanded_name"]
        project = sample_row["project"]
        guide = sample_row["guide_id"]
        amplicon = sample_row["amplicon"]

        if debug:
            print(f"\nProcessing {sample}")

        # depth + window
        depth_row = results_df[results_df["Sample"] == sample]

        if depth_row.empty:
            if debug:
                print(f"No results for {sample}")
            continue

        window = depth_row["QuantWindow"].iloc[0]

        # finding positions in the window automatically
        try:
            positions = find_editable_positions(
                amplicon,
                window,
                nuc_from
            )
        except ValueError as e:
            if debug:
                print(f"Window error for {sample}: {e}")
            continue

        # loading annotation (only once)
        annotation_df = load_annotation(
            sample,
            "outcome/annotation",
            debug
        )

        # counting stats for each position
        
        pos_results = []

        for pos in positions:

            total_mut, specific_mut = get_pos_stats(
                annotation_df,
                pos,
                nuc_from,
                nuc_to
            )

            percent_mut = 0
            percent_specific = 0
            percent_specific_from_mut = 0

            if total_mut > 0:
                percent_mut = total_mut
                percent_specific = specific_mut

                if specific_mut > 0:
                    percent_specific_from_mut = (
                        specific_mut / total_mut * 100
                    )

            pos_results.append({
                "sample": sample,
                "project": project,
                "guide": guide,
                "QuantWindow": window,
                "pos": pos,
                "percent_mut": percent_mut,
                "percent_specific": percent_specific,
                "percent_specific_from_mut": percent_specific_from_mut,
                "amplicon": amplicon
            })

        rows.extend(pos_results)

    return pd.DataFrame(rows)

def get_window_range(amplicon, window):
    """
    Returns the start and end positions of the window in the amplicon (1-based)
    """

    start = amplicon.find(window)

    if start == -1:
        raise ValueError("Window not found in amplicon")

    start += 1
    end = start + len(window) - 1

    return start, end

def get_over_threshold(sample, pos, comparison_df):

    df_sample = comparison_df[comparison_df["sample"] == sample]

    if df_sample.empty:
        return "N/A"

    row = df_sample[df_sample["pos"] == f"mut{pos}"]

    if row.empty:
        return "N/A"

    return row["over_threshold"].iloc[0]
def build_sample_matrix(base_table, control_df, amplicon, window, comparison_df):

    sample = base_table["sample"].iloc[0]

    start, end = get_window_range(amplicon, window)

    positions = list(range(start, end + 1))
    bases = list(window)

    # dictionaries for fast access
    sample_map = base_table.set_index("pos").to_dict("index")

    if control_df is not None and not control_df.empty:
        control_map = control_df.set_index("pos").to_dict("index")
    else:
        control_map = {}

    rows = []

    # string 1 — window
    rows.append(["window"] + bases)

    # string 2 — positions
    rows.append(["pos"] + positions)

    # placeholders
    percent_sample = []
    percent_sample_pos = []
    percent_sample_from_mut = []

    percent_ctrl = []
    percent_ctrl_pos = []
    percent_ctrl_from_mut = []
    
    over_threshold = []
    
    for pos in positions:
    
        value = get_over_threshold(sample, pos, comparison_df)
        over_threshold.append(value)

    for pos in positions:

        if pos in sample_map:
            percent_sample.append(fmt(sample_map[pos]["percent_specific"]))
            percent_sample_pos.append(fmt(sample_map[pos]["percent_mut"]))
            percent_sample_from_mut.append(
                fmt(sample_map[pos]["percent_specific_from_mut"])
            )

        else:
            percent_sample.append("N/A")
            percent_sample_pos.append("N/A")
            percent_sample_from_mut.append("N/A")

        if pos in control_map:
            percent_ctrl.append(fmt(control_map[pos]["percent_specific"]))
            percent_ctrl_pos.append(fmt(control_map[pos]["percent_mut"]))
            percent_ctrl_from_mut.append(
                fmt(control_map[pos]["percent_specific_from_mut"])
            )

        else:
            percent_ctrl.append("N/A")
            percent_ctrl_pos.append("N/A")
            percent_ctrl_from_mut.append("N/A")

    rows.append(["percent_sample"] + percent_sample)
    rows.append(["percent_sample_pos"] + percent_sample_pos)
    rows.append(["percent_sample_from_mut"] + percent_sample_from_mut)

    rows.append(["percent_control"] + percent_ctrl)
    rows.append(["percent_control_pos"] + percent_ctrl_pos)
    rows.append(["percent_control_from_mut"] + percent_ctrl_from_mut)

    rows.append(["over_threshold"] + over_threshold)

    return pd.DataFrame(rows)
    

def build_base_reports(base_table, noise_model, task_df, comparison_df):

    os.makedirs("outcome/base", exist_ok=True)

    for sample, sample_df in base_table.groupby("sample"):

        task_row = task_df[task_df["sample_name"] == sample].iloc[0]

        amplicon = task_row["amplicon"]
        window = sample_df["QuantWindow"].iloc[0]

        project = sample_df["project"].iloc[0]
        guide = sample_df["guide"].iloc[0]

        if noise_model is None or noise_model.empty:
            control_df = pd.DataFrame()
        else:
            control_df = noise_model[
                (noise_model["project"] == project) &
                (noise_model["guide"] == guide) &
                (noise_model["QuantWindow"] == window)
            ]

        matrix = build_sample_matrix(
            sample_df,
            control_df,
            amplicon,
            window,
            comparison_df
        )

        outfile = f"outcome/base/{sample}_base.csv"

        matrix.to_csv(outfile, index=False, header=False)

        print(f"Saved {outfile}")

def project_extension_to_amplicon(amplicon, extension):

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    aln = aligner.align(amplicon, extension)[0]

    aligned_amp, aligned_ext = cf.get_aligned_strings(aln, amplicon, extension)

    projected = []

    for a, b in zip(aligned_amp, aligned_ext):
        if b != "-":
            projected.append(a)

    return "".join(projected).replace("-", "")

def fmt(value):
    if value == "N/A" or pd.isna(value):
        return "N/A"
    return round(float(value), 2)