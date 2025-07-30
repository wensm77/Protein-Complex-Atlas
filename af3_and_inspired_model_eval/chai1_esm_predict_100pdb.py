import logging
import shutil
from pathlib import Path
import sys
import traceback
import time


import numpy as np
# Ensure chai_lab is installed or in PYTHONPATH
try:
    from chai_lab.chai1 import run_inference
except ImportError:
    print("Error: Cannot import 'chai_lab'. Please ensure it is properly installed or you are running in the correct environment.")
    sys.exit(1)

# --- Configuration ---

# 1. Directory containing input .fasta files
INPUT_DIR = "fasta"

# 2. Base directory for storing all prediction results
BASE_OUTPUT_DIR = "chai1_esm"


# --- Main script logic ---

def main():
    """
    Main execution function
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


    input_path = Path(INPUT_DIR)
    base_output_path = Path(BASE_OUTPUT_DIR)

    if not input_path.is_dir():
        logging.error(f"Input directory does not exist -> {input_path}")
        sys.exit(1)
        
    base_output_path.mkdir(parents=True, exist_ok=True)
    fasta_files = list(input_path.glob('*.fasta'))

    if not fasta_files:
        logging.warning(f"No .fasta files found in directory '{input_path}'.")
        return

    logging.info(f"Found {len(fasta_files)} .fasta files, starting prediction...")

    print("Speed test, 200 steps...")
    
    
    # Process each file
    cnt = 0
    time_dict = {}
    for fasta_file in fasta_files:
        # Get filename (e.g., "my_protein")
        cnt += 1
        if cnt >= 21: break
        base_name = fasta_file.stem
        # Create separate output subdirectory
        output_subdir = base_output_path / base_name
        
        print("\n" + "="*60)
        # start_time = time.time()
        logging.info(f"Starting processing: {fasta_file.name}")
        print("="*60)
        
        # Key step: ensure output directory is empty
        if output_subdir.exists():
            logging.warning(f"Removing existing old directory: {output_subdir}")
            shutil.rmtree(output_subdir)
        output_subdir.mkdir(exist_ok=True)
        
        try:
            logging.info(f"Running inference for {fasta_file.name}...")
            # Call run_inference function
            start_time = time.time()
            candidates = run_inference(
                fasta_file=fasta_file,
                output_dir=output_subdir,
                # Use parameters from your example
                num_trunk_recycles=3,
                num_diffn_timesteps=1,
                seed=42,
                device="cuda:2", # Can be modified to "cuda:1" etc. if needed
                use_esm_embeddings=True,
            )
            end_time = time.time()
            # Append elapsed time for this file to file
            with open(base_output_path / "time_1step.txt", "a") as f:
                f.write(f"{base_name}: {end_time - start_time:.2f}\n")

            # Simple processing of results
            if candidates and candidates.cif_paths:
                agg_scores = [rd.aggregate_score.item() for rd in candidates.ranking_data]
                best_score = max(agg_scores)
                logging.info(f"✅ Success: {fasta_file.name} prediction completed. Found {len(candidates.cif_paths)} candidate structures, best score: {best_score:.4f}")

                # You can also load more detailed score files here, e.g.
                # scores_npz_path = output_subdir / "scores.model_idx_0.npz"
                # if scores_npz_path.exists():
                #     scores = np.load(scores_npz_path)
                #     logging.info(f"  - pTM score for first model: {scores['ptm'].item():.4f}")
            else:
                logging.warning(f"Warning: {fasta_file.name} prediction ran but returned no valid candidate structures.")

        except Exception as e:
            logging.error(f"❌ Failed: {fasta_file.name} encountered a serious error during processing.")
            # Log detailed error information to this task's output directory
            error_details = traceback.format_exc()
            logging.error(error_details)
            error_file = output_subdir / "error.log"
            error_file.write_text(error_details)

    # Write time dictionary to file
    

    print("\n" + "="*60)
    logging.info("🎉 All processing completed!")
    print("="*60)

if __name__ == "__main__":
    main()