# Protein-Complex-Atlas
## introduction
Protein-Complex-Atlas is the code repository related to the paper "Atlas of Predicted Protein Complex Structures Across Kingdoms".

## High-quality protein complex dataset Predicted by Colabfold
Protein complexes are fundamental to all biological processes. Public databases have been enriched with millions of potential protein-protein interaction entries for human and model organisms. However, accurately depicting large-scale of protein complex structures across kingdoms, including both the database-curated potential hits and a vast number of unidentified interactions in the prokaryotic genomes, has significantly lagged behind. Here, we constructed an atlas comprising 1.1 million protein-protein interaction structures using AlphaFold2-based ColabFold, encompassing proteome-level interactions from bacteria, archaea, human, mouse, plant, and human-virus pairs. Overall, we identified 181,671 high-confidence protein complex structures, including 37,855 in the human interactome and 108,879 among 224 types of prokaryotic species. Through structural clustering, protein complexes that exhibit conserved structures across multiple kingdoms have been detected, which aids in elucidating their previously unknown functions. Further analysis validated by Co-Immunoprecipitation experiments indicated multiple potential new viral receptors for Human mastadenovirus A and Papiine alphaherpesvirus 2. Additionally, we uncovered widespread gene fusion or fission events throughout evolution by integrating our protein complex structures with the AlphaFold-predicted monomer structure database. Finally, we extended the applicability of our dataset for enhanced protein binding surface prediction and protein sequence design by deep learning models, displayed the integrated value beyond each structure. We anticipate that our large-scale atlas of predicted protein complexes will prove valuable for further research and applications.



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

prepare binding_sites.csv format like this:
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

## af3 and inspired model evaluation


### chai-1 Prediction (ESM version)
```bash
python af3_and_inspiredmodel_eval/chai1_esm_predict.py
```
### chai-1 Prediction (MSA version)
```bash
python af3_and_inspiredmodel_eval/chai1_msa_predict.py
```

### Boltz Prediction
```bash
boltz predict input_path --recycling_steps 3 --sampling_steps 200 --diffusion_samples 5 --use_msa_server
```

### Protenix Prediction (ESM version)
```bash
protenix predict --input input.json --out_dir ./output_no_msa --seeds 101 --use_esm
```
### Protenix Prediction (MSA version)
```bash
bash af3_and_inspiredmodel_eval/inference_batch_msa.sh
```

## License

This project is licensed under the MIT License.