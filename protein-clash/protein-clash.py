import os  
import numpy as np  
import pandas as pd  
import ray  
from Bio.SVDSuperimposer import SVDSuperimposer  
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection  
from Bio.SeqUtils import seq1  
import pymol2  
from tqdm import tqdm  

# Initialize Ray  
ray.init(num_cpus=10)  

def calculate_interface_overlap_with_pymol(pdb_file1, pdb_file2, pdb1_chain, pdb2_chain, distance_threshold=5.0):  
    """Calculate the interface overlap of common chains in two PDB files using all atoms."""  
    with pymol2.PyMOL() as pymol:  
        pymol.cmd.load(pdb_file1, "structure1")  
        pymol.cmd.load(pdb_file2, "structure2")  

        # Select all atoms from specified chains  
        chain1_atoms = pymol.cmd.get_model(f"structure1 and chain {pdb1_chain}").atom  
        chain2_atoms = pymol.cmd.get_model(f"structure2 and chain {pdb2_chain}").atom  

        # Prepare other chains assuming only two chains are present  
        other_chain1 = 'B' if pdb1_chain == 'A' else 'A'  
        other_chain2 = 'B' if pdb2_chain == 'A' else 'A'  

        # Get interface atoms within the distance threshold  
        interface1_atoms = pymol.cmd.get_model(f"structure1 and chain {other_chain1} within {distance_threshold} of (structure1 and chain {pdb1_chain})").atom  
        interface2_atoms = pymol.cmd.get_model(f"structure2 and chain {other_chain2} within {distance_threshold} of (structure2 and chain {pdb2_chain})").atom  

        # Extract coordinates  
        interface1_coords = np.array([atom.coord for atom in interface1_atoms])  
        interface2_coords = np.array([atom.coord for atom in interface2_atoms])  

        # Handle empty cases  
        if interface1_coords.size == 0 or interface2_coords.size == 0:  
            return 0.0  # No overlap possible if either interface is empty  

        # Calculate distances using broadcasting  
        distances = np.linalg.norm(interface1_coords[:, np.newaxis, :] - interface2_coords[np.newaxis, :, :], axis=-1)  

        # Identify overlapping atoms  
        overlap_atoms_1 = np.any(distances <= distance_threshold, axis=1)  # Atoms in interface1 that overlap  
        overlap_atoms_2 = np.any(distances <= distance_threshold, axis=0)  # Atoms in interface2 that overlap  

        # Count unique overlapping atoms  
        overlap_atom_count = min(np.sum(overlap_atoms_1), np.sum(overlap_atoms_2))  

        # Calculate overlap ratio  
        min_interface_atom_count = min(len(interface1_atoms), len(interface2_atoms))  
        overlap_ratio = overlap_atom_count / min_interface_atom_count if min_interface_atom_count > 0 else 0  # Divide by 2 to avoid overcounting  

        # Clean up selections  
        pymol.cmd.delete("structure1")  
        pymol.cmd.delete("structure2")  

    return overlap_ratio  


def analyze_interface_overlap_with_pymol(pdb_file1, pdb_file2, chain1_id="A", chain2_id="B", distance_threshold=5.0):  
    """  
    Analyze the overlap of interaction interfaces between specified chains in two PDB files using PyMOL.  
    
    :param pdb_file1: Path to the first PDB file  
    :param pdb_file2: Path to the second PDB file  
    :param chain1_id: ID of the first chain (default is "A")  
    :param chain2_id: ID of the second chain (default is "B")  
    :param distance_threshold: Distance threshold for defining interface atoms (default is 5.0)  
    :return: Returns the overlap ratio of the two chains (chain1, chain2)  
    """  
    with pymol2.PyMOL() as pymol:  
        pymol.cmd.load(pdb_file1, "protein1")  
        pymol.cmd.load(pdb_file2, "protein2")  
        
        # Select interface atoms of the two chains  
        pymol.cmd.select("interface1_chain1", f"protein1 and chain {chain1_id} within {distance_threshold} of protein2 and chain {chain2_id}")  
        pymol.cmd.select("interface1_chain2", f"protein2 and chain {chain2_id} within {distance_threshold} of protein1 and chain {chain1_id}")  
        
        # Get atom counts  
        overlap_count1 = pymol.cmd.count_atoms("interface1_chain1")  
        overlap_count2 = pymol.cmd.count_atoms("interface1_chain2")  
        
        total_atoms_chain1 = pymol.cmd.count_atoms(f"protein1 and chain {chain1_id}")  
        total_atoms_chain2 = pymol.cmd.count_atoms(f"protein2 and chain {chain2_id}")  
        
        # Calculate overlap ratio  
        overlap_ratio_chain1 = overlap_count1 / total_atoms_chain1 if total_atoms_chain1 > 0 else 0  
        overlap_ratio_chain2 = overlap_count2 / total_atoms_chain2 if total_atoms_chain2 > 0 else 0  

    return (overlap_ratio_chain1 + overlap_ratio_chain2)/2.0  

def load_and_filter_csv(csv_file):  
    df = pd.read_csv(csv_file)  
    print('before filter:', df.shape)  
    # Filter based on plddt and iptm criteria  
    after_df = df[(df['done_tag'] == 'done') & (df['plddt'] >= 70) & (df['iptm'] >= 0.6)]  
    print('after filter:', after_df.shape)  
    return after_df  

def extract_chain_sequence(chain):  
    return ''.join([seq1(residue.resname) for residue in chain.get_residues() if residue.id[0] == ' '])  

def find_matching_row(human_virus_sequence, human_human_df):  
    for _, row in human_human_df.iterrows():  
        if human_virus_sequence in row['sequence']:  
            return row  
    return None  

def find_matching_chain_id(human_virus_sequence, human_human_structure):  
    for chain in human_human_structure.get_chains():  
        human_human_sequence = extract_chain_sequence(chain)  
        if human_virus_sequence in human_human_sequence:  
            return chain.id  
    return None  

def align_and_save(human_human_pdb, human_virus_pdb, human_virus_chain_id, human_human_chain_id, output_dir):  
    parser = PDBParser(QUIET=True)  
    human_human_structure = parser.get_structure('human_human', human_human_pdb)  
    human_virus_structure = parser.get_structure('human_virus', human_virus_pdb)  

    human_human_chain = human_human_structure[0][human_human_chain_id]  
    human_virus_chain = human_virus_structure[0][human_virus_chain_id]  

    human_human_coords = np.array([atom.coord for atom in human_human_chain.get_atoms() if atom.get_name() == 'CA'])  
    human_virus_coords = np.array([atom.coord for atom in human_virus_chain.get_atoms() if atom.get_name() == 'CA'])  

    sup = SVDSuperimposer()  
    sup.set(human_virus_coords, human_human_coords)  
    sup.run()  
    rot, tran = sup.get_rotran()  

    for atom in human_human_structure.get_atoms():  
        atom.transform(rot, tran)  

    io = PDBIO()  
    io.set_structure(human_human_structure)  
    
    # Ensure the output directory exists  
    os.makedirs(output_dir, exist_ok=True)  
    
    output_path = os.path.join(output_dir, os.path.basename(human_virus_pdb)[:15] + os.path.basename(human_human_pdb))  
    print(f"Saving aligned file to {output_path}")  
    io.save(output_path)  
    return output_path  

def check_overlap_ca(coords1, coords2, threshold=5):  
    mat = np.append(coords1, coords2, axis=0)  
    a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]  
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T  
    contact_dists = dists[:len(coords1), len(coords1):]  

    overlap = np.argwhere(contact_dists <= threshold)  

    if len(overlap) == 0:  
        return 0.0  

    unique_overlap_atoms = np.unique(overlap[:, 0] if len(coords1) <= len(coords2) else overlap[:, 1])  
    overlap_amount = len(unique_overlap_atoms) / min(len(coords1), len(coords2))  
    
    return overlap_amount  

@ray.remote  
class ProgressActor:  
    def __init__(self, total):  
        self.total = total  
        self.current = 0  

    def update(self, n):  
        self.current += n  
        print(f"Progress: {self.current}/{self.total} ({(self.current / self.total) * 100:.2f}%)")  

@ray.remote  
def process_single_file(virus_row, human_human_df, aligned_dir, progress_actor):  
    try:  
        parser = PDBParser(QUIET=True)  
        virus_pdb_file = os.path.join(virus_row['batch_result_dir'], virus_row['target_pdb_name'])  
        
        human_virus_structure = parser.get_structure('human_virus', virus_pdb_file)  
        
        human_virus_chain = human_virus_structure[0].child_list[0]  
        human_virus_sequence = extract_chain_sequence(human_virus_chain)  
        
        matching_human_row = find_matching_row(human_virus_sequence, human_human_df)  
        
        result = None  
        if matching_human_row is not None:  
            human_human_pdb_file = os.path.join(matching_human_row['batch_result_dir'], matching_human_row['target_pdb_name'])  
            human_human_structure = parser.get_structure('human_human', human_human_pdb_file)  
            
            human_human_chain_id = find_matching_chain_id(human_virus_sequence, human_human_structure)  
            
            if human_human_chain_id:  
                aligned_human_human_file = align_and_save(human_human_pdb_file, virus_pdb_file, human_virus_chain.id, human_human_chain_id, aligned_dir)  
                
                aligned_human_human_structure = parser.get_structure('human_human_aligned', aligned_human_human_file)  

                human_virus_ca_coords = {  
                    chain.id: np.array([atom.coord for atom in chain.get_atoms() if atom.get_name() == 'CA'])  
                    for chain in human_virus_structure[0] if chain.id != human_virus_chain.id  
                }  
                human_human_ca_coords = {  
                    chain.id: np.array([atom.coord for atom in chain.get_atoms() if atom.get_name() == 'CA'])  
                    for chain in aligned_human_human_structure[0] if chain.id != human_human_chain_id  
                }  

                combined_human_human_coords = np.concatenate([coords for coords in human_human_ca_coords.values()])  
                combined_human_virus_coords = np.concatenate([coords for coords in human_virus_ca_coords.values()])  

                overlap_amount = check_overlap_ca(combined_human_human_coords, combined_human_virus_coords)  
                
                # interface_overlap = analyze_interface_overlap_with_pymol(aligned_human_human_file,  
                #                                                          virus_pdb_file)  
                interface_overlap = calculate_interface_overlap_with_pymol(aligned_human_human_file,  
                                                                           virus_pdb_file,   
                                                                           human_human_chain_id,   
                                                                           human_virus_chain.id)  
                result = {  
                    'human_human_sequence': matching_human_row['sequence'],  
                    'human_virus_sequence': virus_row['sequence'],  
                    'aligned_human_human_file': aligned_human_human_file,  
                    'virus_pdb_file': virus_pdb_file,  
                    'overlap_amount': overlap_amount,  
                    'interface_overlap': interface_overlap,  
                    'plddt': virus_row['plddt'],  
                    'iptm': virus_row['iptm']  
                }  
        
        # Update progress  
        ray.get(progress_actor.update.remote(1))  
        print(result)  
        
        return result  
    except Exception as e:  
        print(f"Error processing {virus_row['target_pdb_name']}: {e}")  
        return None  

def process_files_ray(human_virus_df, human_human_df, output_csv, aligned_dir):  
    total = len(human_virus_df)  
    progress_actor = ProgressActor.remote(total)  
    
    futures = [  
        process_single_file.remote(virus_row, human_human_df, aligned_dir, progress_actor)  
        for _, virus_row in human_virus_df.iterrows()  
    ]  
    
    results = []  
    for i, future in tqdm(enumerate(futures), total=len(futures), desc="Processing files"):  
        result = ray.get(future)  
        if result:  
            results.append(result)  
    
    # Save results to CSV  
    if results:  
        pd.DataFrame(results).to_csv(output_csv, mode='a', header=not os.path.exists(output_csv), index=False)  

    print("All files processed!")  

def merge_csv_files(file_list, output_file):  
    merged_df = pd.DataFrame()  

    for i, file in enumerate(file_list):  
        # Read CSV file  
        df = pd.read_csv(file)  

        # Add suffix to the ID column, suffix is the index of the file in the list (starting from 1)  
        suffix = f'_f{i+1}'  
        df['id'] = df['id'].astype(str) + suffix  

        # Merge dataframes  
        merged_df = pd.concat([merged_df, df], ignore_index=True)  

    df_unique = merged_df.drop_duplicates(subset='sequence')   
    # Save the merged dataframe to a new CSV file  
    df_unique.to_csv(output_file, index=False)  

# Example usage  
# Process human-human dataset  
file_list = ['../data/human_human.csv']  
human_human_complexes_file='merged_human_human_complexes_output.csv'  
merge_csv_files(file_list, human_human_complexes_file)  

# Process human-virus dataset  
file_list = ['../data/human_virus.csv']  
human_virus_complexes_file='merged_human_virus_complexes_output.csv'  
merge_csv_files(file_list, human_virus_complexes_file)  

print("Loading human-human complexes...")  
human_human_df = load_and_filter_csv(human_human_complexes_file)  
print("Loading human-virus complexes...")  
human_virus_df = load_and_filter_csv(human_virus_complexes_file)  

output_csv = 'collision_results.csv'  
aligned_dir = './data/aligned/'  

process_files_ray(human_virus_df, human_human_df, output_csv, aligned_dir)  

# Shut down Ray  
ray.shutdown()