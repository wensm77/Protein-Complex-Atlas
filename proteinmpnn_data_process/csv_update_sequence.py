import pandas as pd
from Bio import SeqIO
import argparse

def parse_fasta(fasta_file):
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def get_sequence(chainid, sequence_dict):
    if chainid in sequence_dict:
        return sequence_dict[chainid]
    
    if not chainid.startswith('uc'):
        uc_id = 'uc' + chainid
        if uc_id in sequence_dict:
            return sequence_dict[uc_id]
    
    if chainid.startswith('uc'):
        no_uc_id = chainid[2:]
        if no_uc_id in sequence_dict:
            return sequence_dict[no_uc_id]
    
    return ''

def main(csv_file, fasta_file, output_csv_file, use_biopython=False):
    df = pd.read_csv(csv_file)
    
    if use_biopython:
        sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        df['SEQUENCE'] = df['CHAINID'].apply(lambda x: str(sequence_dict.get(x, {}).seq) if x in sequence_dict else '')
    else:
        sequence_dict = parse_fasta(fasta_file)
        print(f"Read {len(sequence_dict)} sequences")
        
        df['SEQUENCE'] = df['CHAINID'].apply(lambda x: get_sequence(x, sequence_dict))
        
        missing_count = df['SEQUENCE'].apply(lambda x: len(x) == 0).sum()
        print(f"Missing sequences: {missing_count}")
    
    df.to_csv(output_csv_file, index=False)
    print(f"Updated CSV file saved as {output_csv_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update CSV file with sequences from FASTA file")
    parser.add_argument("--csv_file", required=True, help="Path to input CSV file")
    parser.add_argument("--fasta_file", required=True, help="Path to FASTA file containing sequences")
    parser.add_argument("--output_csv", required=True, help="Path to output CSV file")
    parser.add_argument("--use_biopython", action="store_true", help="Use BioPython for FASTA parsing (default: manual parsing)")
    
    args = parser.parse_args()
    
    main(args.csv_file, args.fasta_file, args.output_csv, args.use_biopython)