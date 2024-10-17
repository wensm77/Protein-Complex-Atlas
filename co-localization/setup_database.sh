#!/bin/bash

# SUCH SCRIPT IS FOR SETTING UP THE DATABASE: QUERY AND TARGET
# QUERY: FASTA FILES SHOULD BE IN THE SAME FOLDER, CONTAINING 36 BACTERIA PATHOGENS FASTA-FILES
# TARGET: FASTA FILES SHOULD BE IN THE SAME FOLDER(NCBI DATABASE HAS BEEN SPLIT INTO 13 PARTS FOR MEMORY CONSIDERATION)
# MMSEQS_DBS: THE FOLDER WHERE THE DATABASE WILL BE STORED
# PLEASE MAKE SURE THAT THE BASENAMES WITHOUT EXTENSION OF FASTA FILES ARE UNIQUE

if [ $# -lt 2 ]; then
  echo "INPUT PARAMETER NUMBERS ERROR! Usage: bash setup_database.sh <query_fasta_dir> <target_fasta_dir> [mmseqs_database_path]"
  exit 1
fi

# GET INPUT FASTA FOLDERS AND DATABASE TARGET PATH
DEFAULT_DBS="$(pwd)/MMSEQS_DBS"
QUERY_FAA=$1
TARGET_FAA=$2
MMSEQS_DBS=${3:-$DEFAULT_DBS} 

# CHECK IF MMSEQS_DBS EXISTS. IF NOT, CREATE IT
if [ ! -d "$MMSEQS_DBS" ]; then
    mkdir -p "$MMSEQS_DBS"
fi

# CREATE QUERY DATABASE
echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Start Creating Query Database!"

cd $MMSEQS_DBS

for file in $QUERY_FAA/*.fa $QUERY_FAA/*.faa $QUERY_FAA/*.fasta; do
    if [[ -f "$file" ]]; then
        file_name=$(basename $file)
        file_extension="${file_name##*.}"
        db_name="query_$(basename $file ".${file_extension}")"
        mmseqs createdb "$file" $db_name
        echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Query $db_name Database Created!"        
    fi
done

echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Finish Creating Query Database!"

# CREATE TARGET DATABASE
echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Start Creating Target Database!"

# WE SPLIT NCBI DATABASE INTO 13 PARTS IN ADVANCE, AND THEY SHARE THE SAME PREFIX "merged_ncbi_", BUT END WITH NUMBERS 1~13

for file in $TARGET_FAA/*.fa $TARGET_FAA/*.faa $TARGET_FAA/*.fasta; do
    if [[ -f "$file" ]]; then
        file_name=$(basename $file)
        file_extension="${file_name##*.}"
        db_name="$(basename $file ".${file_extension}")"
        mmseqs createdb "$file" $db_name
        mmseqs createindex $db_name "tmp_${db_name}"
        rm -rf "tmp_${db_name}"
        echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Target $db_name Database Created!"
    fi
done

echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Finish Creating Target Database!"

