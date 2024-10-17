# Co-Localization User Guide

The co-localization scripts are used to calculate the co-localization score matrices of 36 pathogen bacterias.

## Installation

Download the scripts to local host. 

The runtime environment must contain **MMseqs2 Release 15-6f452** (https://github.com/soedinglab/MMseqs2) and the other necessary packages. We recommend the scripts are processed in `Linux OS` and compiled on `Python 3.9` . Hardware requirements are as follows: 

- 200GB memory is required at least;
- 1TB local storage is required for faster I/O flows and storing databases.



## Getting Started

It is not difficult to use the scripts to calculate the matrices. The steps to run such part are listed:

1. Download the 36 **pathogens** and **NCBI** fasta-files to your local host;
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



## Output Directory Structure

If the scripts are processed successfully, the output directory will be made of 15 directories, including

- 13 directories named as `colocalization_x` to save the co-localzation intermediate results;
- 1 directory named as `m8_files` to save intermediate m8-files;
- 1 directory named as **`Final_Results`** to save all final TSV-format tables.

In the `Final_Results`, 36 pathogen PPIs TSV-format tables are listed, whose header is made of 6 items, such as

```
Pair Index	Sequence 1	Sequence 2	SequenceName 1	SequenceName 2	Min Distance
```

which represent PPI pair index, the first acid sequence, the second acid sequence, the first protein name, the second protein name and the minimum distance of the pairwise proteins, respectively.

