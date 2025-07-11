import pandas as pd
import argparse

def main(input_csv, output_csv):
    df = pd.read_csv(input_csv)

    unique_sequences = df['CHAINID'].unique()
    sequence_to_hash = {seq: idx + 1 for idx, seq in enumerate(unique_sequences)}

    df['HASH'] = df['CHAINID'].map(sequence_to_hash)

    representative_hashes = {seq: sequence_to_hash[seq] for seq in df['CLUSTER'].unique()}

    df['CLUSTER'] = df['CLUSTER'].map(representative_hashes)

    df.to_csv(output_csv, index=False)
    print(f"Updated CSV file with Hash column and replaced Cluster_Representative_ID written to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CSV file to add hash values and update cluster IDs")
    parser.add_argument("--input_csv", required=True, help="Path to input CSV file")
    parser.add_argument("--output_csv", required=True, help="Path to output CSV file")
    
    args = parser.parse_args()
    
    main(args.input_csv, args.output_csv)