import pandas as pd
import subprocess
import os
from pathlib import Path
import noise_evaluation as ne
import crispresso_functions as cf
import generate_tables as gt
from config import RESULTS_PATH, RESULTS_PATH_NOISE

def initialize_pipeline():
    
    RESULTS_PATH.parent.mkdir(exist_ok=True)
    RESULTS_PATH_NOISE.parent.mkdir(exist_ok=True)
    
    # finding fastq 
    fastq_dict = gt.find_fastq_files("fastq")
    

    # Loading tables
    seq_df, samp_df, param_df, base_params = gt.load_tables()
    
    # Checking for the presence of the noise column
    if "noise" not in samp_df.columns:
        print("No noise column in samples.csv")
        print("Adding noise column with default 'no' value")
        samp_df["noise"] = "no"
        samp_df.to_csv("samples.csv", index=False)
    
    # Expansion of tasks
    task_df = gt.expand_samples(samp_df, seq_df)
    
    print(f" Initialization completed. Found tasks: {len(task_df)}")
    return fastq_dict, seq_df, samp_df, param_df, base_params, task_df




def run_noise_mode(fastq_dict, seq_df, param_df, params, tasks):
    """
    Noise mode for tasks with noise=yes
    """
    print("\n" + "="*50)
    print(" STARTING NOISE MODE")
    print("="*50)
    
    results = []
    
    # Separating by edit_mode
    cas9_base_tasks = tasks[tasks['edit_mode'].str.lower().isin(['cas9', 'base'])]
    prime_tasks = tasks[tasks['edit_mode'].str.lower() == 'prime']
    
    print(f"\n Cas9/Base: {len(cas9_base_tasks)}")
    print(f"   Prime: {len(prime_tasks)}")

    # Processing Cas9/Base 
    if not cas9_base_tasks.empty:
        cas9_base_results = run_cas9_base_noise(
            cas9_base_tasks, fastq_dict, seq_df, params, param_df
        )
        results.extend(cas9_base_results)


    # Processing PRIME 
    if not prime_tasks.empty:
        prime_results = run_prime_fallback_noise(
            prime_tasks, fastq_dict, seq_df, params, param_df
        )
        results.extend(prime_results)
    
    # Saving all results once
    if results:
        save_noise_results(results)
        print(f"\n Saved {len(results)} results")    
    
    # Post-processing
    if results:
        print(f"\n Processed samples: {len(results)}")
        print(f"   Final results: {[r.get('Sample', 'Unknown') for r in results]}")
        
        print("\n" + "="*50)
        print(" LAUNCHING POST-PROCESSING")
        print("="*50)
        
        try:
            run_noise_post_processing(tasks, seq_df, params)
        except Exception as e:
            print(f" Error in post-processing: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("\n No results to process")
    
    return results






def run_fallback_extension(sample_name, task_dict, guide_seq, is_extension, 
                          fastq_dict, seq_df, params, param_df):
    """Launching extension for fallback"""
    
    cmd = cf.build_fallback_command(
        task=task_dict,
        seq_df=seq_df,
        params=params,
        param_df=param_df,
        fastq_dict=fastq_dict,
        # report_dir="report",
        guide_seq=guide_seq,
        is_extension=is_extension
    )
    
    print(f"   Launching extension...")
    
    try:
        subprocess.run(cmd, check=True)
        
        # Checking the result
        report_dir = f"report/CRISPResso_on_{sample_name}_fallback_extension"
        alleles_file = os.path.join(report_dir, "Alleles_frequency_table.zip")
        
        if os.path.exists(alleles_file):
            print(f"   Extension OK")
            return True
        else:
            print(f"   Extension: no result file found")
            return False
    except subprocess.CalledProcessError as e:
        print(f"  Extension error: {e}")
        # checking file
        report_dir = f"report/CRISPResso_on_{sample_name}_fallback_extension"
        alleles_file = os.path.join(report_dir, "Alleles_frequency_table.zip")
        if os.path.exists(alleles_file):
            print(f"   Extension failed, but data is available")
            return True
        return False


def run_fallback_spacer(sample_name, task_dict, guide_seq, is_extension,
                       fastq_dict, seq_df, params, param_df):
    """Launching spacer for fallback"""
    
    cmd = cf.build_fallback_command(
        task=task_dict,
        seq_df=seq_df,
        params=params,
        param_df=param_df,
        fastq_dict=fastq_dict,
        # report_dir="report",
        guide_seq=guide_seq,
        is_extension=is_extension
    )

    print(f"   Launching spacer...")
    
    try:
        subprocess.run(cmd, check=True)
        print(f"   Spacer OK")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  Spacer error: {e}")
        return False


def run_fallback_for_samples(failed_samples, fastq_dict, seq_df,
                             param_df, params, task_df):
    results = []

    for task_dict in failed_samples:
        sample_name = task_dict.get("expanded_name")
        guide_id = task_dict.get("guide_id")
        fallback_name = f"{sample_name}_fallback"

        print(f"\n Fallback for {sample_name}")

        if cf.is_fallback_sample_done(sample_name):
            print(f" {sample_name} already processed through fallback")

            if RESULTS_PATH.exists():
                existing_df = pd.read_csv(RESULTS_PATH)
                if fallback_name in existing_df["Sample"].values:
                    print(" Already in results, skipping")
                    continue

            try:
                res = cf.analyze_fallback_sample(sample_name, task_df, seq_df)
                if res:
                    results.append(res)
                    cf.update_results_table(res)
                    print(f"Fallback results updated for {sample_name}")
            except Exception as e:
                print(f" Error analyzing fallback: {e}")
            continue

        # Getting sequences
        try:
            g = seq_df.set_index("guide_id").loc[guide_id]
        except KeyError:
            print(f" Guide {guide_id} not found")
            continue

        spacer = g.get("spacer", "").strip()
        extension = g.get("extension", "").strip()
        amplicon = task_dict.get("amplicon", "")

        rev_comp = cf.revcomp(extension)
        processed_extension = cf.project_extension_to_amplicon(amplicon, rev_comp)

        success_ext = run_fallback_extension(
            sample_name, task_dict, guide_seq=processed_extension,
            is_extension=True, fastq_dict=fastq_dict,
            seq_df=seq_df, params=params, param_df=param_df
        )
        if not success_ext:
            continue

        run_fallback_spacer(
            sample_name, task_dict, guide_seq=spacer,
            is_extension=False, fastq_dict=fastq_dict,
            seq_df=seq_df, params=params, param_df=param_df
        )

        try:
            res = cf.analyze_fallback_sample(sample_name, task_df, seq_df)
            if res:
                results.append(res)
                cf.update_results_table(res)
                print(f"Fallback results updated for {sample_name}")
        except Exception as e:
            print(f"Error in fallback analysis: {e}")

    return results


def save_standard_results(results):
    """Saving results of the standard mode"""
    if not results:
        return
    
    if RESULTS_PATH.exists():
        df = pd.read_csv(RESULTS_PATH)
        df = pd.concat([df, pd.DataFrame(results)], ignore_index=True)
    else:
        df = pd.DataFrame(results)
    
    df.to_csv(RESULTS_PATH, index=False)
    print(f" Saved {len(results)} results to {RESULTS_PATH}")

def save_noise_results(results):
    """Saving results of the noise mode"""
    print(f"\n START SAVING: {len(results)} results")
    for r in results:
        print(f"   - {r.get('Sample', 'Unknown')}")
    
    if not results:
        return
    
    if RESULTS_PATH_NOISE.exists():
        existing_df = pd.read_csv(RESULTS_PATH_NOISE)
        print(f"   In the file there are already: {existing_df['Sample'].tolist()}")
    else:
        existing_df = pd.DataFrame()
        print(f"   File not found")
    
    # Saving
    new_df = pd.DataFrame(results)
    combined_df = pd.concat([existing_df, new_df], ignore_index=True)
    combined_df.to_csv(RESULTS_PATH_NOISE, index=False)
    
    print(f"   SAVED. Now in the file: {combined_df['Sample'].tolist()}")
    print(f" END SAVING\n")

def run_pipeline(samples_file='samples.csv'):
    """
    Running the pipeline directly in the notebook
    Mode (standard/noise) is determined by the 'noise' column in samples.csv
    
    Parameters:
    -----------
    samples_file : str
        path to samples.csv (if already filled)
    """
    
    # Initialization
    print("\n Step 1: Initialization")
    fastq_dict, seq_df, samp_df, param_df, params, task_df = initialize_pipeline()
    
    # Check filling of samples.csv
    print("\n Step 2: Check samples.csv")
    print("\nFirst 5 rows of samples.csv:")
    # display(samp_df.head())
    
    # Check filling
    empty_cols = []
    for col in ['guide_id', 'group', 'edit_mode', 'project']:
        if samp_df[col].isna().all() or (samp_df[col] == '').all():
            empty_cols.append(col)
    
    if empty_cols:
        print(f"\n Attention! Columns {empty_cols} are not filled!")
        print("Please fill in the samples.csv file before running the pipeline")
        print("Currently in samples.csv:")
        # display(samp_df)
        
        response = input("\nContinue with current data? (y/n): ")
        if response.lower() != 'y':
            return
    
    # Information about tasks
    print("\n Step 3: Distribution of tasks")
    print(f"Total tasks: {len(task_df)}")
    print(f"noise = yes: {len(task_df[task_df['noise'] == 'yes'])} (will go to noise pipeline)")
    print(f"noise = no: {len(task_df[task_df['noise'] == 'no'])} (will go to standard pipeline)")
    
    print("\n By editing modes:")
    print(task_df['edit_mode'].value_counts())
    
    # Running the processing depending on noise
    print("\n Step 4: Running the processing")
    
    start_time = pd.Timestamp.now()
    
    # Separate tasks by noise
    standard_tasks = task_df[task_df['noise'] == 'no']
    noise_tasks = task_df[task_df['noise'] == 'yes']
    
    all_results = []
    
    # Processing standard tasks (noise = no)
    if not standard_tasks.empty:
        print(f"\n Processing standard tasks (noise=no): {len(standard_tasks)}")
        std_results = run_standard_mode(fastq_dict, seq_df, param_df, params, standard_tasks, task_df)
        all_results.extend(std_results)
        results_file_std = RESULTS_PATH
    
    # Processing noise tasks (noise = yes)
    if not noise_tasks.empty:
        print(f"\n🧬 Processing noise tasks (noise=yes): {len(noise_tasks)}")
        noise_results = run_noise_mode(fastq_dict, seq_df, param_df, params, noise_tasks)
        all_results.extend(noise_results)
        results_file_noise = RESULTS_PATH_NOISE
    
    end_time = pd.Timestamp.now()
    duration = end_time - start_time
    
    # Results
    print("\n" + "="*70)
    print(" PIPELINE COMPLETED")
    print("="*70)
    
    print(f"\n Execution time: {duration}")
    print(f" Total samples processed: {len(all_results)}")
    
    # Show standard results
    if RESULTS_PATH.exists():
        print(f"\n Standard results saved to: {RESULTS_PATH}")
        results_df = pd.read_csv(RESULTS_PATH)
        print(f"   Records: {len(results_df)}")
        if not results_df.empty:
            print(f"\n First 5 standard results:")
            # display(results_df.head())
    
    # Show noise results
    if RESULTS_PATH_NOISE.exists():
        print(f"\n Noise results saved to: {RESULTS_PATH_NOISE}")
        results_df_noise = pd.read_csv(RESULTS_PATH_NOISE)
        print(f"   Rows: {len(results_df_noise)}")
        # if not results_df_noise.empty:

        #     # display(results_df_noise.head())
    
    return all_results


def run_standard_mode(fastq_dict, seq_df, param_df, params, tasks, task_df):
    """
    Standard mode for tasks with noise=no
    """
    print("\n" + "="*50)
    print(" RUNNING STANDARD MODE")
    print("="*50)
    
    results = []
    failed_samples = []
    
    for _, task in tasks.iterrows():
        task_dict = task.to_dict()
        sample_name = task_dict["expanded_name"]
        print(f"\n Processing: {sample_name}")
        
        # Check for existing results
        if cf.is_sample_done(sample_name):
            print(f" {sample_name} CRISPResso report found")
            
            if RESULTS_PATH.exists():
                existing_df = pd.read_csv(RESULTS_PATH)
                if sample_name in existing_df["Sample"].values:
                    print(f" {sample_name} already analyzed, skipping completely")
                    continue
            
            try:
                res = cf.analyze_sample(sample_name)
                if res:
                    results.append(res)
                    cf.update_results_table(res)
                    print(f" Analysis done and results updated for {sample_name}")
            except Exception as e:
                print(f"⚠ Failed to analyze {sample_name}: {e}")
            
            continue
        
        # Forming the CRISPResso command
        try:
            cmd = cf.build_command(
                task=task_dict,
                seq_df=seq_df,
                params=params,
                param_df=param_df,
                fastq_dict=fastq_dict
            )
        except Exception as e:
            print(f" Failed to build CRISPResso command for {sample_name}: {e}")
            continue
        
        # Running CRISPResso
        try:
            subprocess.run(cmd, check=True)
            print(f" CRISPResso finished for {sample_name}")
        except subprocess.CalledProcessError as e:
            print(f" CRISPResso failed for {sample_name} (exit {e.returncode})")
            # Save only errors 6 for fallback
            if task_dict["edit_mode"].lower() == "prime" and e.returncode == 6:
                failed_samples.append(task_dict)
            continue
        
        # Analyzing results
        try:
            res = cf.analyze_sample(sample_name)
            if res:
                results.append(res)
                cf.update_results_table(res)
                print(f" Analysis done and results updated for {sample_name}")
        except Exception as e:
            print(f" Failed to analyze sample {sample_name}: {e}")
        
        # Save intermediate results
        if results:
            try:
                pd.DataFrame(results).to_csv(RESULTS_PATH, index=False)
            except Exception as e:
                print(f" Failed to save intermediate summary: {e}")
    
    # Fallback for failed Prime samples
    if failed_samples:
        print(f"\n {len(failed_samples)} samples failed, starting fallback run")
        fallback_results = run_fallback_for_samples(
            failed_samples, fastq_dict, seq_df, 
            param_df = param_df,
            params = params, 
            task_df = task_df, 
        )
        results.extend(fallback_results)
    
    return results



def run_noise_post_processing(task_df, seq_df, params):
    """
    Full post-processing for noise samples
    """
    
    # Making necessary directories
    ANNOTATION_DIR = Path("outcome/annotation")
    ALLELES_NOISE_DIR = Path("outcome/alleles")
    BASE_DIR = Path("outcome/base")
    PROCESSED_DIR = Path("outcome/processed")
    for d in [ANNOTATION_DIR, ALLELES_NOISE_DIR, BASE_DIR, PROCESSED_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Annotation of all samples
    print("\n Step 1: Annotation of variants...")
    try:
        # Temporarily redefine global variables
        import sys
        module = sys.modules[__name__]
        module.task_df = task_df
        module.seq_df = seq_df
        
        
        ne.process_all_samples_annotation(task_df, seq_df)
        print("    Annotation completed")
    except Exception as e:
        print(f"    Annotation error: {e}")
        import traceback
        traceback.print_exc()
        
    #  Building noise model for reference 
    print("\n Step 2: Building noise model...")
    try:
        results_df = pd.read_csv(RESULTS_PATH_NOISE)
        
        # Building control table
        control_table, projects_without_controls = ne.build_control_position_table(
            ANNOTATION_DIR, 
            task_df, 
            results_df
        )
        
        # Building noise model
        noise_df = ne.build_noise_model(control_table)
        noise_df.to_csv("outcome/processed/noise_model.csv", index=False)
        print(f"    Noise model saved: {len(noise_df)} positions")
        
        # Saving projects without controls for final summary
        global projects_without_noise
        projects_without_noise = projects_without_controls
    except Exception as e:
        print(f"    Error building noise model: {e}")
        projects_without_noise = set()
    
    # Checking samples against thresholds (reference) 
    print("\n Step 3: Checking reference samples...")
    comparison_df = pd.DataFrame()  # <- default value
    try:
        comparison_df = ne.check_samples_vs_threshold(task_df, noise_df)
#        comparison_df = comparison_df[comparison_df["sample_sum"] > 0.1]
        
        # Cleaning multiple positions
        clean_df = ne.clean_multi_positions(comparison_df)
        clean_df.to_csv("outcome/processed/sample_vs_threshold.csv", index=False)
                # Annotate each sample using the threshold table
        print("\n Step 4: Annotating samples with thresholds...")

        annotation_files = list(ANNOTATION_DIR.glob("*_positions.txt"))

        for annotation_path in annotation_files:
            try:
                # getting sample name
                df_annotation, _, sample_name = ne.load_annotation_file(
                    annotation_path
                )

                # adding statuses based on clean_df
                df_annotated = ne.annotate_with_thresholds(
                    df_annotation=df_annotation,
                    df_thresholds=clean_df,
                    sample_name=sample_name
                )

                # rewrite annotation
                df_annotated.to_csv(
                    annotation_path,
                    sep="\t",
                    index=False
                )

                print(f"    Annotated: {sample_name}")

            except Exception as e:
                print(
                    f"    Annotation error for "
                    f"{annotation_path.name}: {e}"
                )

        print("    Threshold annotation completed")

    except Exception as e:
        print(f"    Error checking samples against thresholds: {e}")
        import traceback
        traceback.print_exc()

    # Build the final noise summary
    print("\n Step 5: Building noise summary...")

    results_df = pd.read_csv(RESULTS_PATH_NOISE)

    # === step 5 in run_noise_post_processing ===
    print("\n Step 5: Building noise summary...")
    results_df = pd.read_csv(RESULTS_PATH_NOISE)
    
    # projects without controls 
    final_summary_df = ne.build_noise_summary(
        annotation_dir=ANNOTATION_DIR,
        results_df=results_df,
        task_df=task_df,
        info_projects_without_noise=projects_without_noise
    )
    
    summary_path = "outcome/noise_summary.csv"
    final_summary_df.to_csv(
        summary_path,
        index=False
    )
    print(
        f"    Noise summary saved: {summary_path} "
        f"({len(final_summary_df)} samples)"
    )

    # --- Step 7: Processing Base editors ---
    base_tasks = task_df[task_df['edit_mode'] == 'base']
    if not base_tasks.empty:
        control_table = ne.build_control_table(
            task_df=task_df,
            results_df=results_df,
            params=params,
            load_annotation=ne.load_annotation,
            debug=True
        )
        
        print(control_table)
        
        noise_model = ne.build_base_noise_model(control_table, debug=True)
        
        print(noise_model)
        
        base_table = ne.build_sample_base_table(
            task_df=task_df,
            results_df=results_df,
            params=params,
            load_annotation=ne.load_annotation,
            debug=True
        )
        
        print(base_table.head())
        print("\n DEBUG: Checking comparison_df before build_base_reports")
        print(f"   comparison_df type: {type(comparison_df)}")
        print(f"   comparison_df shape: {comparison_df.shape if hasattr(comparison_df, 'shape') else 'N/A'}")
        if hasattr(comparison_df, 'columns'):
            print(f"   comparison_df columns: {comparison_df.columns.tolist()}")
            if not comparison_df.empty:
                print(f"   comparison_df head:\n{comparison_df.head()}")
        ne.build_base_reports(
            base_table=base_table,
            noise_model=noise_model,
            task_df=task_df,
            comparison_df=comparison_df
        )
    else:
        print("\n Step 7: Base editors not found")

def process_single_annotation(alleles_file, task_df, seq_df, results_df):
    """
    Processing a single annotation file
    """
    #  Names
    base_with_fallback = alleles_file.name.replace("_Alleles_frequency_table.txt", "")
    base_without_fallback = base_with_fallback.replace("_fallback", "")
    
    print(f"\n=== Processing {base_with_fallback} ===")
    
    # Retrieve project_id 
    task_row = task_df[task_df["expanded_name"] == base_without_fallback]
    if task_row.empty:
        print(" No project found in task_df")
        return
    
    project_id = task_row.iloc[0]["project"]
    
    # Retrieve amplicon + expected 
    seq_row = seq_df[seq_df["project_id"] == project_id]
    if seq_row.empty:
        print("No sequences found")
        return
    
    amplicon = seq_row.iloc[0]["amplicon"]
    expected = seq_row.iloc[0]["expected_seq"]
  
    # Retrieve window from results.csv
    res_row = results_df[results_df["Sample"] == base_with_fallback]
    if res_row.empty:
        print(" No window in results.csv")
        return
    
    window_seq = res_row.iloc[0]["QuantWindow"]
    
    # Load alleles 
    df = ne.load_alleles_table(alleles_file)
    
    # Filtering 
    hdr_mask = df["Reference_Name"].isin(["HDR", "Prime-edited"]) & \
               (df["Read_Status"] == "MODIFIED")
    ref_mask = ~hdr_mask
    df = df[(~ref_mask) | ((df["n_deleted"] <= 5) & (df["n_inserted"] <= 5))]
    df = ne.filter_alleles(df)
    
    # Preparation of HDR mapping
    hdr_map = None
    expected_changes = []
    
    if isinstance(expected, str) and expected.strip() != "":
        aln_amp, aln_exp = cf.align_expected_to_amplicon(amplicon, expected)
        expected_changes = ne.detect_expected_changes(aln_amp, aln_exp)
        hdr_map = ne.build_hdr_amplicon_map(aln_amp, aln_exp)
    
    # Retrieve window positions
    window_positions = ne.window_positions_in_reference(window_seq, amplicon)
    if not window_positions:
        print("⚠ Window does not map")
        return
    
    # Variants annotation
    df = ne.annotate_window_variants_stat(
        df=df,
        positions=window_positions,
        hdr_map=hdr_map,
        expected_changes=expected_changes
    )
    
    # Gap annotation within the window
    gap_positions = ne.find_window_gaps(window_positions)
    if gap_positions:
        df = ne.annotate_gap_variants(df, gap_positions)
    
    df["pos_type"] = df["pos_type"].apply(ne.collapse_positions)
    
    # Sorting and saving
    df = df.sort_values(by="%Reads", ascending=False)
    df = ne.annotate_hdr_window(
        df=df,
        amplicon=amplicon,
        expected=expected,
    )
    out_path = Path("outcome/annotation") / f"{base_with_fallback}_positions.txt"
    df.to_csv(out_path, sep="\t", index=False)
    print(f" Saved {out_path}")


def run_cas9_base_noise(tasks, fastq_dict, seq_df, params, param_df):
    """Processing Cas9/Base samples in noise mode"""
    results = []
    
    for _, task in tasks.iterrows():
        task_dict = task.to_dict()
        sample_name = task_dict["expanded_name"]
        print(f"\n Cas9/Base: {sample_name}")
        
        # Check if the sample is already processed
        if cf.is_sample_done(sample_name):
            if RESULTS_PATH_NOISE.exists():
                existing_df = pd.read_csv(RESULTS_PATH_NOISE)
                if sample_name in existing_df["Sample"].values:
                    print(f"{sample_name} already processed, skipping")
                    continue
            
            try:
                res = cf.analyze_sample(sample_name)
                if res:
                    results.append(res)
                continue
            except Exception as e:
                print(f"Analysis error: {e}")
        
        # cas9 mode for base 
        if task_dict["edit_mode"].lower() == "base":
            task_dict["edit_mode"] = "cas9"
        
        # Formation of the command
        try:
            cmd = cf.build_command(
                task=task_dict,
                seq_df=seq_df,
                params=params,
                param_df=param_df,  
                fastq_dict=fastq_dict
            )
        except Exception as e:
            print(f"Error building command: {e}")
            continue
        
        # Launch CRISPResso
        try:
            subprocess.run(cmd, check=True)
            print(f"CRISPResso executed")
        except Exception as e:
            print(f"Error launching: {e}")
            continue
        
        # Analysis
        try:
            res = cf.analyze_sample(sample_name)
            if res:
                results.append(res)
        except Exception as e:
            print(f"Error analyzing: {e}")
    
    return results


def run_prime_fallback_noise(tasks, fastq_dict, seq_df, params, param_df):
    """Processing Prime samples in noise mode (immediate fallback)"""
    results = []
    
    for _, task in tasks.iterrows():
        task_dict = task.to_dict()
        sample_name = task_dict["expanded_name"]
        guide_id = task_dict.get("guide_id")
        
        print(f"\n Prime fallback: {sample_name}")

        # Check if the fallback sample is already processed
        if cf.is_fallback_sample_done(sample_name):
            print(f" {sample_name} already processed through fallback")
            if RESULTS_PATH_NOISE.exists():
                existing_df = pd.read_csv(RESULTS_PATH_NOISE)
                if f"{sample_name}_fallback" in existing_df["Sample"].values:
                    print(f" Already in results, skipping")
                    continue
            
            try:
                res = cf.analyze_fallback_sample(sample_name, tasks, seq_df)
                if res:
                    results.append(res)
                continue
            except Exception as e:
                print(f" Error analyzing fallback: {e}")
        
        # Retrieve sequences
        try:
            g = seq_df.set_index("guide_id").loc[guide_id]
        except KeyError:
            print(f" Guide {guide_id} not found")
            continue
        
        spacer = g.get("spacer", "").strip()
        extension = g.get("extension", "").strip()
        amplicon = task_dict.get("amplicon", "")
        
        # Processing extension
        rev_comp = cf.revcomp(extension)
        processed_extension = cf.project_extension_to_amplicon(amplicon, rev_comp)
        
        # Launch extension
        success_ext = run_fallback_extension(
            sample_name, task_dict, guide_seq=processed_extension,
            is_extension=True, fastq_dict=fastq_dict,
            seq_df=seq_df, params=params, param_df=param_df
        )
        
        if not success_ext:
            continue
        
        # Launch spacer
        run_fallback_spacer(
            sample_name, task_dict, guide_seq=spacer,
            is_extension=False, fastq_dict=fastq_dict,
            seq_df=seq_df, params=params, param_df=param_df
        )
        
        # Analysis
        try:
            res = cf.analyze_fallback_sample(sample_name, tasks, seq_df)
            if res:
                results.append(res)
        except Exception as e:
            print(f" Error analyzing fallback: {e}")
    
    return results



# FUNCTIONS FOR GENERATING THE TEMPLATE

def generate_samples_template():
    """Separate function only for creating the samples.csv template"""
    print(" Generating samples.csv template...")
    
    fastq_dict = gt.find_fastq_files("fastq")
    
    df = pd.DataFrame({
        "sample_name": sorted(fastq_dict.keys()),
        "guide_id": "",
        "group": "",
        "edit_mode": "",
        "project": "",
        "noise": "no"  
    })
    
    df.to_csv("samples.csv", index=False)
    print("Template created: samples.csv")
    print("   Fill in the table and save the file")
    print("\n Columns:")
    print("   - sample_name: sample name (from fastq)")
    print("   - guide_id: guide ID from sequences.csv")
    print("   - group: 0 for controls, 1/2/3 for samples")
    print("   - edit_mode: cas9 / base / prime")
    print("   - project: ID from sequences.csv")
    print("   - noise: no - standard, yes - noise")
    
    return df


def quick_check_samples(samples_file='samples.csv'):
    """Quick check of samples.csv"""
    if not Path(samples_file).exists():
        print(f" {samples_file} not found!")
        print(" Run generate_samples_template() to create the template")
        return
    
    df = pd.read_csv(samples_file, sep=None, engine="python", dtype=str, keep_default_na=False)
    print(f" {samples_file}:")
    # display(df)
    
    print("\n Checking:")
    for col in df.columns:
        empty = df[col].isna().sum() + (df[col] == '').sum()
        filled = len(df) - empty

   
    print(f"\n   noise = yes: {len(df[df['noise'] == 'yes'])}")
    print(f"   noise = no: {len(df[df['noise'] == 'no'])}")
