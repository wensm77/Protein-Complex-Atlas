#!/usr/bin/env python3
"""
高性能批处理脚本：使用30个核心处理fig4-human和fig4-virus目录中的PDB文件
"""

import os
import sys
import subprocess
import logging
from datetime import datetime
import time
import shutil
from pathlib import Path

# 设置日志
log_filename = f"high_performance_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def get_pdb_files(directory):
    """获取目录中所有PDB文件"""
    pdb_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.pdb'):
                pdb_files.append(os.path.join(root, file))
    return pdb_files

def create_batch_directory(batch_files, batch_dir):
    """创建批次目录并复制文件"""
    os.makedirs(batch_dir, exist_ok=True)
    
    for i, pdb_file in enumerate(batch_files):
        dest_file = os.path.join(batch_dir, f"batch_{i:04d}_{os.path.basename(pdb_file)}")
        shutil.copy2(pdb_file, dest_file)
    
    return len(batch_files)

def process_batch(batch_dir, output_csv, threshold=5.0, processes=30):
    """处理一个批次，使用30个进程"""
    try:
        cmd = [
            'python3', 'identify_interface_residues_sasa_interact.py',
            batch_dir,
            output_csv,
            '--processes', str(processes)
        ]
        
        logger.info(f"处理批次: {batch_dir} -> {output_csv} (使用{processes}个进程)")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1小时超时
        )
        
        if result.returncode == 0:
            if os.path.exists(output_csv):
                size = os.path.getsize(output_csv)
                with open(output_csv, 'r') as f:
                    lines = sum(1 for _ in f)
                logger.info(f"✓ 批次完成: {lines} 行, {size} 字节")
                return True, lines - 1  # 减去表头
            else:
                logger.warning(f"输出文件不存在: {output_csv}")
                return False, 0
        else:
            logger.error(f"✗ 批次失败: 返回码 {result.returncode}")
            if result.stderr:
                logger.error(f"错误: {result.stderr[-500:]}")
            return False, 0
            
    except subprocess.TimeoutExpired:
        logger.error(f"✗ 批次超时: {batch_dir}")
        return False, 0
    except Exception as e:
        logger.error(f"✗ 批次异常: {e}")
        return False, 0

def merge_csv_files(csv_files, final_output):
    """合并多个CSV文件"""
    try:
        with open(final_output, 'w') as outfile:
            header_written = False
            total_rows = 0
            
            for csv_file in csv_files:
                if os.path.exists(csv_file):
                    with open(csv_file, 'r') as infile:
                        lines = infile.readlines()
                        if lines:
                            if not header_written:
                                outfile.write(lines[0])  # 写入表头
                                header_written = True
                            
                            # 写入数据行（跳过表头）
                            for line in lines[1:]:
                                outfile.write(line)
                                total_rows += 1
        
        logger.info(f"合并完成: {final_output} ({total_rows} 行数据)")
        return total_rows
        
    except Exception as e:
        logger.error(f"合并文件时发生错误: {e}")
        return 0

def process_directory_high_performance(input_dir, final_output, batch_size=1000, threshold=5.0, processes=30):
    """高性能分批处理目录"""
    logger.info(f"开始高性能处理: {input_dir}")
    logger.info(f"使用 {processes} 个进程")
    
    # 获取所有PDB文件
    pdb_files = get_pdb_files(input_dir)
    total_files = len(pdb_files)
    
    if total_files == 0:
        logger.error(f"目录 {input_dir} 中没有找到PDB文件")
        return False
    
    logger.info(f"找到 {total_files} 个PDB文件")
    logger.info(f"批次大小: {batch_size}")
    
    # 计算批次数量
    num_batches = (total_files + batch_size - 1) // batch_size
    logger.info(f"将分成 {num_batches} 个批次处理")
    
    # 创建工作目录
    work_dir = f"hp_batch_work_{os.path.basename(input_dir)}"
    os.makedirs(work_dir, exist_ok=True)
    
    csv_files = []
    total_processed = 0
    total_interface_residues = 0
    start_time = time.time()
    
    try:
        for batch_idx in range(num_batches):
            batch_start_time = time.time()
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, total_files)
            batch_files = pdb_files[start_idx:end_idx]
            
            logger.info(f"\n=== 批次 {batch_idx + 1}/{num_batches} ===")
            logger.info(f"处理文件 {start_idx + 1}-{end_idx}/{total_files}")
            
            # 创建批次目录
            batch_dir = os.path.join(work_dir, f"batch_{batch_idx:03d}")
            batch_csv = os.path.join(work_dir, f"batch_{batch_idx:03d}.csv")
            
            # 复制文件到批次目录
            file_count = create_batch_directory(batch_files, batch_dir)
            logger.info(f"创建批次目录: {batch_dir} ({file_count} 个文件)")
            
            # 处理批次
            success, residue_count = process_batch(batch_dir, batch_csv, threshold, processes)
            
            batch_elapsed = time.time() - batch_start_time
            
            if success:
                csv_files.append(batch_csv)
                total_processed += file_count
                total_interface_residues += residue_count
                files_per_sec = file_count / batch_elapsed if batch_elapsed > 0 else 0
                logger.info(f"✓ 批次 {batch_idx + 1} 成功: {residue_count} 个界面残基")
                logger.info(f"  批次用时: {batch_elapsed:.1f} 秒 ({files_per_sec:.1f} 文件/秒)")
            else:
                logger.error(f"✗ 批次 {batch_idx + 1} 失败")
            
            # 清理批次目录
            if os.path.exists(batch_dir):
                shutil.rmtree(batch_dir)
            
            # 显示总体进度和预估时间
            progress = (batch_idx + 1) / num_batches * 100
            elapsed = time.time() - start_time
            if batch_idx > 0:  # 避免除零错误
                avg_time_per_batch = elapsed / (batch_idx + 1)
                remaining_batches = num_batches - (batch_idx + 1)
                eta = remaining_batches * avg_time_per_batch
                logger.info(f"总体进度: {progress:.1f}% ({total_processed}/{total_files} 文件)")
                logger.info(f"预计剩余时间: {eta:.1f} 秒 ({eta/60:.1f} 分钟)")
        
        # 合并所有CSV文件
        if csv_files:
            logger.info(f"\n=== 合并结果 ===")
            final_rows = merge_csv_files(csv_files, final_output)
            
            # 清理临时CSV文件
            for csv_file in csv_files:
                if os.path.exists(csv_file):
                    os.remove(csv_file)
            
            total_elapsed = time.time() - start_time
            avg_files_per_sec = total_processed / total_elapsed if total_elapsed > 0 else 0
            
            logger.info(f"✓ 高性能处理完成: {final_output}")
            logger.info(f"  处理文件: {total_processed}/{total_files}")
            logger.info(f"  界面残基: {final_rows}")
            logger.info(f"  总用时: {total_elapsed:.1f} 秒 ({total_elapsed/60:.1f} 分钟)")
            logger.info(f"  平均速度: {avg_files_per_sec:.1f} 文件/秒")
            
            return True
        else:
            logger.error("没有成功的批次，无法生成最终结果")
            return False
            
    finally:
        # 清理工作目录
        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)

def main():
    """主函数"""
    logger.info("开始高性能分批处理 (30核心)")
    
    # 定义任务
    # 添加时间戳到输出文件名
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    tasks = [
        {
            'input_dir': 'fig4-human',
            'output_csv': f'fig4_human_interface_residues_hp_{timestamp}.csv',
            'description': 'Human蛋白质复合物'
        },
        {
            'input_dir': 'fig4-virus', 
            'output_csv': f'fig4_virus_interface_residues_hp_{timestamp}.csv',
            'description': 'Virus蛋白质复合物'
        }
    ]
    
    # 高性能参数设置
    batch_size = 1000  # 增大批次大小
    threshold = 5.0
    processes = 30     # 使用30个进程
    
    logger.info(f"批次大小: {batch_size}")
    logger.info(f"ΔASA阈值: {threshold}%")
    logger.info(f"并行进程数: {processes}")
    logger.info(f"日志文件: {log_filename}")
    
    success_count = 0
    start_time = time.time()
    
    for i, task in enumerate(tasks, 1):
        logger.info(f"\n{'='*60}")
        logger.info(f"任务 {i}/{len(tasks)}: {task['description']}")
        logger.info(f"{'='*60}")
        
        if not os.path.isdir(task['input_dir']):
            logger.error(f"输入目录不存在: {task['input_dir']}")
            continue
        
        if process_directory_high_performance(
            task['input_dir'], 
            task['output_csv'], 
            batch_size, 
            threshold, 
            processes
        ):
            success_count += 1
    
    # 总结
    elapsed = time.time() - start_time
    logger.info(f"\n{'='*60}")
    logger.info(f"高性能分批处理完成")
    logger.info(f"{'='*60}")
    logger.info(f"成功完成: {success_count}/{len(tasks)} 个任务")
    logger.info(f"总用时: {elapsed:.1f} 秒 ({elapsed/60:.1f} 分钟)")
    logger.info(f"日志文件: {log_filename}")
    
    # 显示最终结果
    for task in tasks:
        if os.path.exists(task['output_csv']):
            size = os.path.getsize(task['output_csv'])
            try:
                with open(task['output_csv'], 'r') as f:
                    lines = sum(1 for _ in f)
                logger.info(f"输出: {task['output_csv']} ({lines} 行, {size} 字节)")
            except:
                logger.info(f"输出: {task['output_csv']} ({size} 字节)")
    
    return 0 if success_count == len(tasks) else 1

if __name__ == '__main__':
    sys.exit(main())