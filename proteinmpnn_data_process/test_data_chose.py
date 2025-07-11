import pandas as pd
import random

def main():
  
    input_csv = 'list_new.csv'
    df = pd.read_csv(input_csv)
    
  
    all_clusters = df['CLUSTER'].unique().tolist()
    
 
    if len(all_clusters) < 500:
        raise ValueError("Not enough unique clusters to sample 250 for both validation and test sets.")
    
   
    valid_clusters = random.sample(all_clusters, 250)
    
    
    remaining_clusters = list(set(all_clusters) - set(valid_clusters))
    test_clusters = random.sample(remaining_clusters, 250)
    
     
    with open('valid_clusters.txt', 'w') as f:
        for cluster in valid_clusters:
            f.write(f"{cluster}\n")
    
    with open('test_clusters.txt', 'w') as f:
        for cluster in test_clusters:
            f.write(f"{cluster}\n")
    
    print("Clusters have been written to valid_clusters.txt and test_clusters.txt")

if __name__ == "__main__":
    main()
