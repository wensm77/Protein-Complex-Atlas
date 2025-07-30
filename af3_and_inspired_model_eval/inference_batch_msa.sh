# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

export CUTLASS_PATH="/mnt/sda/wensm/Protenix/opt/cutlass"
export LAYERNORM_TYPE=fast_layernorm
export USE_DEEPSPEED_EVO_ATTENTION=true

N_sample=5
N_step=200
N_cycle=10
seed=101

JSON_FILES_DIR="json"

DUMP_BASE_DIR="msa_result"

# 3. 检查JSON目录是否存在
if [ ! -d "$JSON_FILES_DIR" ]; then
    echo "error: JSON directory '$JSON_FILES_DIR' does not exist."
    exit 1
fi

# 4. 遍历目录下的所有 .json 文件
for json_file in ${JSON_FILES_DIR}/*.json; do
    

    filename=$(basename -- "$json_file")
    dirname="${filename%.*}"
    
    current_dump_dir="${DUMP_BASE_DIR}/${dirname}"
    
    echo "=========================================================="
    echo "==> start processing: ${filename}"
    echo "==> output directory: ${current_dump_dir}"
    echo "=========================================================="
    
    # 执行预测命令
    CUDA_VISIBLE_DEVICES=1 python3 runner/inference.py \
    --seeds ${seed} \
    --dump_dir "${current_dump_dir}" \
    --input_json_path "${json_file}" \
    --model.N_cycle ${N_cycle} \
    --sample_diffusion.N_sample ${N_sample} \
    --sample_diffusion.N_step ${N_step}

done

echo "finished!"


# The following is a demo to use DDP for inference
# torchrun \
#     --nproc_per_node $NPROC \
#     --master_addr $WORKER_0_HOST \
#     --master_port $WORKER_0_PORT \
#     --node_rank=$ID \
#     --nnodes=$WORKER_NUM \
#     runner/inference.py \
#     --seeds ${seed} \
#     --dump_dir ${dump_dir} \
#     --input_json_path ${input_json_path} \
#     --model.N_cycle ${N_cycle} \
#     --sample_diffusion.N_sample ${N_sample} \
#     --sample_diffusion.N_step ${N_step}