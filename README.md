# Protein-Complex-Atlas
## introduction
Protein-Complex-Atlas is the code repository related to the paper "Atlas of Predicted Protein Complex Structures Across Kingdoms".

## High-quality protein complex dataset Predicted by Colabfold
we constructed an atlas comprising **1.15 million** protein-protein interaction structures using ColabFold(https://github.com/YoshitakaMo/localcolabfold),encompassing proteome-level interactions from bacteria, archaea, human, mice, plant, and human-virus pairs. Overall, we identified **170,001** high-confidence protein complex structures, especially over 33,226 in the human interactome.

We will be opening all high-quality protein complex structure data, as listed in the table below! 
**Data avaliable**: [Bio view]:http://www.biostructurehub-zhejianglab.com/

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


## Analysis of Potential Spatial Clashes in Predicted Protein Complexes
This document provides an overview of analyzing potential spatial clashes involving predicted human protein-protein complexes and human protein-virus protein complexes.

### Code Location:
File Path: protein-clash/protein_clash.py

### Example Code Execution:

cd ./protein-clash
python protein_analysis_ray_v2.py

### Main Workflow of the Code:

The main workflow of the code involves loading and filtering structural data of human protein-protein and human protein-virus protein complexes from CSV files, using pLDDT and iptm scores to select high-quality models. It then aligns sequences and superimposes structures by searching for matching human protein sequences in the human protein-protein complex database. If a match is found, the viral protein structure is superimposed onto the human complex to align the human proteins closely. The code calculates spatial distances between the viral protein and other human proteins in the complex, identifying potential spatial clashes when distances fall below a predefined threshold. Finally, the results, including clash information and structural details, are saved to a CSV file.