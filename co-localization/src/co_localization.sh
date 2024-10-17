#!/bin/bash

if [ $# -lt 3 ]; then
    echo "INPUT PARAMETER NUMBERS ERROR! Usage: bash co_localization.sh <query_fasta_dir> <output_dir> <mmseqs_database_path> [threads_number]"
    exit 1
fi

if [[ "$4" =~ ^-?[0-9]+$ ]]; then
    echo "TYPEERROR! Usage: threads_number MUST BE integer!"
    exit 1
fi


DEFAULT_THREADS=15

QUERY_FAA=$1   # The query fasta files directory (the same as the script "setup_database.sh")
OUTPUT_DIR=$2  # The output directory
MMSEQS_DBS=$3  # The mmseqs_databases directory created by setup_database.sh
THREADS=${4:-$DEFAULT_THREADS}  # The number of threads


sorted_files=$(ls "$QUERY_FAA"/*.faa | sort -t_ -k1,1n)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
M8_FILE_DIR="${OUTPUT_DIR}/m8_files"

cd $MMSEQS_DBS

for ((i=1; i<=13; i++)); do
    targetDB="merged_ncbi_$i"
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]: PART ${i}/13 IS STARTING!"
    mmseqs touchdb $targetDB

    for file in $sorted_files; do
        base_name="$(basename $file .faa)"
        db_name="query_${base_name}"
        out_suffix="${base_name}"
        aln_file="aln_${i}_${base_name}"

        mmseqs search "${db_name}" $targetDB "${aln_file}" "tmp$i" --db-load-mode 2 --max-seqs 10000
        mmseqs convertalis "${db_name}" $targetDB "${aln_file}" "${aln_file}.m8" --db-load-mode 2
        echo "[$(date +"%Y-%m-%d %H:%M:%S")]: ${base_name} in ${targetDB} search finished!"

        python "${SCRIPT_DIR}/extract_identity_from_result.py" -p "${aln_file}.m8" -i 0.3 -a 1.0 -o ${out_suffix} -s $M8_FILE_DIR
        co_folder="${OUTPUT_DIR}/colocalization_$i"
        curr_folder="${co_folder}/Job_${out_suffix}"
        if [ -d "${curr_folder}" ]; then
            echo "[$(date +"%Y-%m-%d %H:%M:%S")]: the '$curr_folder' has been exist, skip this step..."
        else
            mkdir -p "${curr_folder}"
        fi

        m8_file="${M8_FILE_DIR}/${out_suffix}_0.3.m8"
        python "${SCRIPT_DIR}/process_adpoints_mmseqs_mp.py" -c $m8_file -l 0 -p $THREADS -r 0.9 -o "${curr_folder}"
        wait

        python "${SCRIPT_DIR}/get_co_ad_sparse_rdd.py" -w "${curr_folder}"
        python "${SCRIPT_DIR}/figures_plot_latest1.py" -f ${file} -w "${curr_folder}" -o "${out_suffix}" -c 2000 -m 0.4 -s 0

        mmseqs rmdb $aln_file
        echo "[$(date +"%Y-%m-%d %H:%M:%S")]: ${base_name} in ${targetDB} co-localization finished!"
    done

    rm -rf "tmp$i"
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]: PART ${i}/13 IS FINISHED!"
done

python "${SCRIPT_DIR}/union_tsv.py" -o "${OUTPUT_DIR}/Final_Results"
echo "[$(date +"%Y-%m-%d %H:%M:%S")]: Colocalization finished!"
