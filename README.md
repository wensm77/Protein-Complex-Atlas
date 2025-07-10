# Protein-Complex-Atlas
## introduction
Protein-Complex-Atlas is the code repository related to the paper "Atlas of Predicted Protein Complex Structures Across Kingdoms".

## High-quality protein complex dataset Predicted by Colabfold
we constructed an atlas comprising **1.15 million** protein-protein interaction structures using ColabFold (https://github.com/YoshitakaMo/localcolabfold) , encompassing proteome-level interactions from bacteria, archaea, human, mice, plant, and human-virus pairs. Overall, we identified **170,001** high-confidence protein complex structures, especially over **33,226** in the human interactome.

We will be opening all high-quality protein complex structure data, as listed in the table below! 

| Type        | Source               | number of high-quality protein complex |
| ----------- | -------------------- | :---------------------------------: |
| Heterodimer | Bacteria and Archaea | 70393                                |
| Heterodimer | Human                | 33226                                |
| Heterodimer | Mouse                | 11529                                |
| Heterodimer | Arabidopsis          | 6281                                 |
| Heterodimer | Human-Virus          | 4646                                 |
| Homodimer   | Bacteria             | 35585                                |
| Homodimer   | Virulence factor     | 8082                                 |
| Homodimer   | Virus                | 259                                  |
|             | Total                | 170001                               |


## Co-Localization User Guide

The co-localization scripts are used to calculate the co-localization score matrices of 36 pathogen bacterias.

### Installation

Download the scripts to local host. 

The runtime environment must contain **MMseqs2 Release 15-6f452** (https://github.com/soedinglab/MMseqs2) and the other necessary packages. We recommend the scripts are processed in `Linux OS` and compiled on `Python 3.9` . Hardware requirements are as follows: 

- 200GB memory is required at least;
- 1TB local storage is required for faster I/O flows and storing databases.



### Getting Started

It is not difficult to use the scripts to calculate the matrices. The steps to run such part are listed:

1. Download the 36 **pathogens** and **NCBI** fasta-files to your local host ( https://huggingface.co/datasets/Protein-Complex-Lab/NCBI_Proteins )；
2. Run `setup_database.sh` to prepare all necessary  databases;
3. Run `src/co_localization.sh` to calculate co-localization score matrices and find all potential protein-protein interactions (PPIs) for downstream predicting workflow.

Fasta-files consist of two parts: 36 pathogens as query databases and 13 split NCBI databases as target databases.

`setup_database.sh` needs two required input parameters, which are query fasta directory and target fasta directory. The third parameter is optional, and it represents the directory of `MMseqs2` databases. If the third not provided, it will be set as `pwd/MMSEQS_DBS`.  

```bash
bash setup_database.sh /path/of/pathogens_fastas /path/of/ncbi_fastas /path/of/MMSEQS_DBS
```

Based on the above preparation, `src/co_localization.sh` can be processed. It will accept 3 required and 1 optional parameters:

- The query fasta-files directory (the same as the script `setup_database.sh`)
- The output directory
- The mmseqs_databases directory created by `setup_database.sh`
- (optional) The numbers of threads are used [Default = 15].

```bash
bash src/co_localization.sh /path/of/pathogens_fastas /path/of/output /path/of/MMSEQS_DBS 15
```

**Notice the first and the third parameters of 2 bash scripts must be identical!**

`setup_database.sh` and `src/co_localization.sh` will both run for two days approximately. We will improve them later.



### Output Directory Structure

If the scripts are processed successfully, the output directory will be made of 15 directories, including

- 13 directories named as `colocalization_x` to save the co-localzation intermediate results;
- 1 directory named as `m8_files` to save intermediate m8-files;
- 1 directory named as **`Final_Results`** to save all final TSV-format tables.

In the `Final_Results`, 36 pathogen PPIs TSV-format tables are listed, whose header is made of 6 items, such as

```
Pair Index	Sequence 1	Sequence 2	SequenceName 1	SequenceName 2	Min Distance
```

which represent PPI pair index, the first acid sequence, the second acid sequence, the first protein name, the second protein name and the minimum distance of the pairwise proteins, respectively.



## Analysis of Potential Spatial Clashes in Predicted Protein Complexes
This document provides an overview of analyzing potential spatial clashes involving predicted human protein-protein complexes and human protein-virus protein complexes.

### Code Location:
File Path: protein-clash/protein_clash.py

### Example Code Execution:

cd ./protein-clash
python protein_analysis_ray_v2.py

### Main Workflow of the Code:

The main workflow of the code involves:
1. load and filtering structural data of human protein-protein and human protein-virus protein complexes from CSV files, using pLDDT and iptm scores to select high-quality models. 
2. align sequences and superimposes structures by searching for matching human protein sequences in the human protein-protein complex database. If a match is found, the viral protein structure is superimposed onto the human complex to align the human proteins closely. The code calculates spatial distances between the viral protein and other human proteins in the complex, identifying potential spatial clashes when distances fall below a predefined threshold. 
3. Finally, the results, including clash information and structural details, are saved to a CSV file.

## ProteinMPNN Data Processing

You should clone ProteinMPNN repository first.
```bash
git clone https://github.com/dauparas/ProteinMPNN.git
```

###  Data Processing Pipeline

### 1. PDB Preprocessing

```bash
# Rename PDB files
python mapping_new_id.py

# Convert PDB to CIF format (multi-threaded)
python convert_pdb_to_cif_multiThread.py
```

### 2. Sequence Clustering

```bash
# Extract sequences from PDB files
python pdb2fasta.py

# Cluster sequences using Foldseek
foldseek easy-multimercluster pdb/ clu tmp --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65

# Convert to CSV format
python tsv2csv.py
```

### 3. Data Processing

```bash
# Update sequence information and metadata
python csv_update_sequence.py
python csv_fix.py
python final.csv.py

# Generate train/validation/test splits
python test_data_chose.py
```

### Training

Use `training.py` to fine-tune the model.

```bash
python ProteinMPNN/training/training.py \
    --path_for_training_data <path_to_your_data_folder> \
    --path_for_outputs <path_to_your_experiment_output> \
    --previous_checkpoint <path_to_pretrained_model.pt>
```

### Inference

To generate sequences for a new PDB structure, use `protein_mpnn_run.py`.

```bash
python ProteinMPNN/protein_mpnn_run.py \
    --pdb_path <path_to_input.pdb> \
    --pdb_path_chains "A" \
    --out_folder <path_to_output_folder> \
    --num_seq_per_target 8
```

## MPBind

### Installation

please ensure that your runtime environment meets the requirements.

### Environment Setup

1. **Clone the repository:**
```bash
git clone https://github.com/jianlin-cheng/MPBind.git
cd MPBind
```

2. **Create and activate conda environment:**
```bash
conda env create -f MPBind.yml
conda activate MPBind
```

### Training

MPBind supports training custom models on your own datasets. The training process is based on MPBind.


### Prediction

MPBind predicts multiple types of binding sites on protein structures from PDB files. We only focus on protein-protein binding sites.

### Usage

```bash
# Activate environment
conda activate MPBind

# Navigate to experiment folder
cd experiment

# Basic prediction command
python inference.py --input [pdb_folder] --output [prediction_folder] --version 2 --binding_type [0]
```

### Output Format

The binding site predictions are stored in the **B-factor column** of the output PDB files. Prediction scores range from 0.0 to 1.0, where higher values indicate higher binding probability.


### Evaluation

To evaluate prediction performance, we use a threshold of **0.5** on the output B-factor scores to classify each residue as an interface (1) or non-interface (0) residue. Standard metrics such as Precision, Recall, AUC and F1-score are then calculated by comparing these binary predictions against the ground truth labels using the `sklearn.metrics` library.


## Data Deduplication Scripts

Includes two main scripts for ensuring data quality and preventing data leakage between training and test sets by removing structurally similar proteins.

### 1. Chain-level Deduplication

The `chain_deduplicate.py` script filters a dataset to remove any protein chains that belong to the same structural cluster as chains in another dataset.

**Usage:**

```bash
python deduplicate/chain_deduplicate.py \
    --cluster_file <path_to_cluster.tsv> \
    --train_file <path_to_train_data.csv> \
    --test_file <path_to_test_data.csv> \
    --train_id_column "CHAINID" \
    --test_id_column "CHAINID" \
    --output_filtered_test "filtered_test_set.csv" \
    --output_overlap_report "overlap_report.csv"
```

### 2. Surface-level (Interface) Deduplication

The `surface_deduplicate.py` script scans Foldseek search results (in `.m8` format) to identify pairs of proteins that have significant overlap in their binding interface regions. This helps filter out proteins that may be functionally redundant even if their overall structure is different.

prepare the format like this:
```
ID,Label
ID_123, 0011...
```

**Usage:**

```bash
python deduplicate/surface_deduplicate.py \
    --m8_dir <path_to_folder_with_m8_files> \
    --binding_sites_file <path_to_binding_sites.csv> \
    --output_file "similar_surface_ids.csv"
```