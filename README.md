# Protein-Complex-Atlas
## introduction
Protein Complex Atlas is the code repository related to the paper "Atlas of Predicted Protein Complex Structures Across Kingdoms".


## Analysis of Potential Spatial Clashes in Predicted Protein Complexes
This document provides an overview of analyzing potential spatial clashes involving predicted human protein-protein complexes and human protein-virus protein complexes.

### Code Location:
File Path: protein-clash/protein_clash.py

### Example Code Execution:

cd ./protein-clash
python protein_analysis_ray_v2.py

### Main Workflow of the Code:

The main workflow of the code involves loading and filtering structural data of human protein-protein and human protein-virus protein complexes from CSV files, using pLDDT and iptm scores to select high-quality models. It then aligns sequences and superimposes structures by searching for matching human protein sequences in the human protein-protein complex database. If a match is found, the viral protein structure is superimposed onto the human complex to align the human proteins closely. The code calculates spatial distances between the viral protein and other human proteins in the complex, identifying potential spatial clashes when distances fall below a predefined threshold. Finally, the results, including clash information and structural details, are saved to a CSV file.