import os
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor

error_log_file = "error_pdb_files.txt"

def pdb_to_fasta(pdb_file):
    parser = PDBParser()
    sequences = []
    
    try:
        structure = parser.get_structure('', pdb_file)

        for model in structure:
            for chain in model:
                
                sequence = ''.join([residue.resname for residue in chain if residue.id[0] == ' '])
                sequence = sequence.replace('ALA', 'A').replace('CYS', 'C').replace('ASP', 'D').replace('GLU', 'E')
                sequence = sequence.replace('PHE', 'F').replace('GLY', 'G').replace('HIS', 'H').replace('ILE', 'I')
                sequence = sequence.replace('LYS', 'K').replace('LEU', 'L').replace('MET', 'M').replace('ASN', 'N')
                sequence = sequence.replace('PRO', 'P').replace('GLN', 'Q').replace('ARG', 'R').replace('SER', 'S')
                sequence = sequence.replace('THR', 'T').replace('VAL', 'V').replace('TRP', 'W').replace('TYR', 'Y')

            
                fasta_sequence = Seq(sequence)
                record = SeqRecord(fasta_sequence, id=f"{os.path.basename(pdb_file)}_Chain_{chain.id}", description=f"Chain {chain.id} from {os.path.basename(pdb_file)}")
                sequences.append(record)

    except Exception:
        with open(error_log_file, "a") as log_file:
            log_file.write(f"{pdb_file}\n")  
        return []  

    return sequences


def process_pdb_files_parallel(pdb_files):
    with ProcessPoolExecutor(max_workers=65) as executor:  
        all_sequences = list(executor.map(pdb_to_fasta, pdb_files))
    return all_sequences


def collect_pdb_files(root_folder):
    pdb_files = []
    for dirpath, _, filenames in os.walk(root_folder):
        
        for file_name in filenames:
            if file_name.endswith(".pdb"):
                pdb_files.append(os.path.join(dirpath, file_name))
                
    return pdb_files


def convert_pdbs_in_folders_to_fasta(root_folder, output_folder, max_sequences_per_file=10000, pdb_batch_size=10000):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    pdb_files = collect_pdb_files(root_folder)

    file_count = 1
    sequence_count = 0
    all_sequences = []

    for i in range(0, len(pdb_files), pdb_batch_size):
        batch_files = pdb_files[i:i + pdb_batch_size]
        
        sequences_list = process_pdb_files_parallel(batch_files)

        print("{:.2%}".format(i/len(pdb_files)))
        for sequences in sequences_list:
            all_sequences.extend(sequences)
            sequence_count += len(sequences)
            
            #if sequence_count >= max_sequences_per_file:
    output_fasta_file = os.path.join(output_folder, f"all_sequences.fasta")
    SeqIO.write(all_sequences, output_fasta_file, "fasta")
    print(f"Saved {output_fasta_file} with {len(all_sequences)} sequences.")  
        #all_sequences = []
        #sequence_count = 0
        #file_count += 1

#    if all_sequences:
#        output_fasta_file = os.path.join(output_folder, f"all_sequences_part_{file_count}.fasta")
#        SeqIO.write(all_sequences, output_fasta_file, "fasta")
#        print(f"Saved {output_fasta_file} with {len(all_sequences)} sequences.")


if __name__ == "__main__":
    root_folder = ".../MPNN_data/pdb"  
    output_folder = ".../MPNN_data"

    if os.path.exists(error_log_file):
        os.remove(error_log_file)

    convert_pdbs_in_folders_to_fasta(root_folder, output_folder)
