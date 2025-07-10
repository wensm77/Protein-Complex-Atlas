import pandas as pd
import argparse
import sys
from tqdm import tqdm
import os
import glob

def load_binding_sites(path):
    """
    Loads binding site data, represented as '01' strings.

    Args:
        path (str): The path to the CSV file containing 'chain_id' and 'binding_site' columns.

    Returns:
        dict: A dictionary mapping chain_id to its binding site string (e.g., '0101...').
    """
    print(f"Loading binding site information from {path}...")
    try:
        df = pd.read_csv(path)
        if 'chain_id' not in df.columns or 'binding_site' not in df.columns:
            print(f"Error: Binding site file {path} must contain 'chain_id' and 'binding_site' columns.", file=sys.stderr)
            sys.exit(1)
        
        df['binding_site'] = df['binding_site'].fillna('').astype(str)
        binding_sites_map = df.set_index('chain_id')['binding_site'].to_dict()
        print(f"Loaded {len(binding_sites_map)} binding site records.")
        return binding_sites_map
    except FileNotFoundError:
        print(f"Error: Binding site file not found at {path}", file=sys.stderr)
        sys.exit(1)

_interface_indices_cache = {}
def get_interface_indices(binding_site_string):
    """
    Converts a '01' string to a set containing the 1-based indices of the '1's.
    """
    if binding_site_string in _interface_indices_cache:
        return _interface_indices_cache[binding_site_string]
    
    indices = {i + 1 for i, char in enumerate(binding_site_string) if char == '1'}
    _interface_indices_cache[binding_site_string] = indices
    return indices

def check_similarity(q_id, t_id, q_st, q_ed, t_st, t_ed, binding_sites_map):
    """
    Performs the core similarity check logic based on interface overlap.
    """
    q_site_str = binding_sites_map.get(q_id)
    t_site_str = binding_sites_map.get(t_id)
    
    if q_site_str is None or t_site_str is None:
        return False

    region_q = set(range(q_st, q_ed + 1))
    region_t = set(range(t_st, t_ed + 1))
    interface_q = get_interface_indices(q_site_str)
    interface_t = get_interface_indices(t_site_str)

    intersection_qt = region_q.intersection(interface_t)
    # Using max() in the denominator makes the similarity score stricter.
    max_len_qt = max(len(region_q), len(interface_t))
    
    score_qt = 0.0
    if max_len_qt > 0:
        score_qt = len(intersection_qt) / max_len_qt
    if score_qt > 0.5:
        return True

    intersection_tq = region_t.intersection(interface_q)
    # Using max() in the denominator makes the similarity score stricter.
    max_len_tq = max(len(region_t), len(interface_q))

    score_tq = 0.0
    if max_len_tq > 0:
        score_tq = len(intersection_tq) / max_len_tq
    if score_tq > 0.5:
        return True

    return False

def file_has_similarity(foldseek_path, binding_sites_map):
    """
    Checks a single Foldseek result file for any similar protein pairs.
    Returns True immediately upon finding the first match for efficiency.
    """
    try:
        if os.path.getsize(foldseek_path) == 0:
            return False

        chunk_iterator = pd.read_csv(
            foldseek_path, 
            sep='\t', 
            header=None, 
            chunksize=10000,
            usecols=[0, 1, 4, 5, 6, 7],
            names=['q_id', 't_id', 'q_st', 'q_ed', 't_st', 't_ed'],
            on_bad_lines='skip'
        )
        
        for chunk in chunk_iterator:
            for _, row in chunk.iterrows():
                if check_similarity(
                    row['q_id'], row['t_id'],
                    row['q_st'], row['q_ed'],
                    row['t_st'], row['t_ed'],
                    binding_sites_map
                ):
                    return True  # Found a match, return early.
    
    except (pd.errors.EmptyDataError, ValueError):
        # Catch errors from empty files or mismatched columns.
        return False
    except Exception as e:
        print(f"\nWarning: An error occurred while processing file {os.path.basename(foldseek_path)}: {e}", file=sys.stderr)
        return False
        
    return False

def main(m8_directory, binding_sites_path, output_csv_path):
    """
    Main function: Iterates through all .m8 files in a directory, filters for file IDs
    containing similar protein pairs, and saves them to a CSV file.
    """
    binding_sites_map = load_binding_sites(binding_sites_path)
    
    m8_files = glob.glob(os.path.join(m8_directory, '*.m8'))
    
    if not m8_files:
        print(f"Error: No .m8 files found in the directory '{m8_directory}'.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Found {len(m8_files)} .m8 files in {m8_directory}. Starting processing...")
    
    positive_ids = []
    
    for m8_file in tqdm(m8_files, desc="Processing .m8 files"):
        if file_has_similarity(m8_file, binding_sites_map):
            # Extract the ID from the filename, e.g., '.../ID_2717.m8' -> 'ID_2717'
            file_id = os.path.basename(m8_file).split('.')[0]
            positive_ids.append(file_id)
            
    print(f"\nFiltering complete. Found {len(positive_ids)} files containing similar protein pairs.")
    
    if positive_ids:
        print(f"Writing the list of IDs to {output_csv_path}...")
        df_results = pd.DataFrame(positive_ids, columns=['id'])
        df_results.to_csv(output_csv_path, index=False)
        print("Write complete.")
    else:
        print("No files with similar protein pairs were found. No output file will be created.")

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description="Scans all .m8 files in a directory to find files containing protein pairs that are deemed similar based on interface overlap, then outputs their IDs to a CSV file."
    )
    parser.add_argument(
        "--m8_dir",
        required=True, 
        help="Directory containing the Foldseek search result files (in .m8 format)."
    )
    parser.add_argument(
        "--binding_sites_file",
        required=True, 
        help="CSV file containing the binding site information for all chains."
    )
    parser.add_argument(
        "--output_file", 
        required=True,
        help="Output CSV file path to save the IDs of .m8 files that contain similar proteins."
    )
    
    args = parser.parse_args()
    
    main(args.m8_dir, args.binding_sites_file, args.output_file)