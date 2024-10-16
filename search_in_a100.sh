#!/bin/bash

# 检查参数数量是否正确
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 nodeIndex type num1"
    exit 1
fi

# 获取命令行参数,四个参数分别为：节点编号，homo/heter, 输入编号1，输入编号2
node_index=$1
type=$2
num1=$3

# 定义通用参数
DATABASE_PATH="/mnt/tianwen-course/tianwen/home/David/mmseqs_db/mmseqs_dbs"

# 定义一个函数来运行colabfold_search
# 注意，需要注意文件路径
run_search () {
    local num=$1
    local INPUTFILE="/mnt/nfs-mdc/ai4s/qixianzhi/workstation_search/data/${type}/${type}_search_input_${num}.csv"
    local OUTPUTDIR="/mnt/nfs-mdc/ai4s/qixianzhi/workstation_search/search_work_${node_index}/output_${num}"
    local LOGFILE="/mnt/nfs-mdc/ai4s/qixianzhi/workstation_search/search_work_${node_index}/log_csv_${type}_${num}.log"
    
    start_time=$(date +%s)
    
    colabfold_search \
      --use-env 1 \
      --use-templates 1 \
      --db-load-mode 2 \
      --db2 pdb100_230517 \
      "${INPUTFILE}" \
      "${DATABASE_PATH}" \
      "${OUTPUTDIR}" &> "${LOGFILE}"
    
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    echo "Execution time for chunk ${num}: $elapsed_time seconds" | tee -a "${LOGFILE}"
}

# 启动两个独立的搜索 机器共192CPU，而单个search进程占用64线程
run_search $num1 &

# 等待所有后台进程完成
wait
echo "All searches are complete."
