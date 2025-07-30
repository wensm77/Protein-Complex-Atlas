#!/usr/bin/env python3
"""
批处理脚本：分别处理 fig4-human 和 fig4-virus 文件夹中的PDB文件
识别界面残基并生成对应的CSV文件
"""

import os
import sys
import subprocess
import logging
from pathlib import Path

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_interface_analysis(input_dir: str, output_csv: str, distance_cutoff: float = 5.0):
    """
    运行界面残基分析
    """
    script_path = "/mnt/hdd/yecheng/eval_prot_dimers/identify_interface_residues.py"
    
    if not os.path.exists(script_path):
        logger.error(f"分析脚本不存在: {script_path}")
        return False
    
    if not os.path.exists(input_dir):
        logger.error(f"输入目录不存在: {input_dir}")
        return False
    
    try:
        cmd = ["python3", script_path, input_dir, output_csv, "--distance", str(distance_cutoff)]
        logger.info(f"执行命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"成功完成分析: {output_csv}")
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"分析失败: {e}")
        logger.error(f"错误输出: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"执行分析时出错: {e}")
        return False

def main():
    base_dir = "/mnt/hdd/yecheng/eval_prot_dimers"
    
    # 定义输入和输出路径
    tasks = [
        {
            "name": "fig4-human",
            "input_dir": os.path.join(base_dir, "fig4-human"),
            "output_csv": os.path.join(base_dir, "fig4_human_interface_residues.csv")
        },
        {
            "name": "fig4-virus", 
            "input_dir": os.path.join(base_dir, "fig4-virus"),
            "output_csv": os.path.join(base_dir, "fig4_virus_interface_residues.csv")
        }
    ]
    
    distance_cutoff = 5.0  # 5 Å 距离阈值
    
    logger.info(f"开始批处理界面残基分析")
    logger.info(f"距离阈值: {distance_cutoff} Å")
    
    success_count = 0
    
    for task in tasks:
        logger.info(f"\n{'='*50}")
        logger.info(f"处理任务: {task['name']}")
        logger.info(f"输入目录: {task['input_dir']}")
        logger.info(f"输出文件: {task['output_csv']}")
        logger.info(f"{'='*50}")
        
        success = run_interface_analysis(
            task['input_dir'], 
            task['output_csv'], 
            distance_cutoff
        )
        
        if success:
            success_count += 1
            logger.info(f"✓ 任务 {task['name']} 完成")
        else:
            logger.error(f"✗ 任务 {task['name']} 失败")
    
    logger.info(f"\n{'='*50}")
    logger.info(f"批处理完成")
    logger.info(f"成功: {success_count}/{len(tasks)} 个任务")
    logger.info(f"{'='*50}")
    
    if success_count == len(tasks):
        logger.info("所有任务都成功完成！")
        return 0
    else:
        logger.error(f"有 {len(tasks) - success_count} 个任务失败")
        return 1

if __name__ == "__main__":
    sys.exit(main())