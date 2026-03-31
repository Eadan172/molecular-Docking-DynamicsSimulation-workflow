#!/usr/bin/env python3
"""
RMSD分析脚本
计算蛋白骨架和配体的RMSD随时间的变化

使用方法:
    python rmsd.py -f trajectory.xtc -s system.tpr -o rmsd_analysis.png

输出:
    - RMSD随时间变化图
    - RMSD分布直方图
    - 统计信息 (均值、标准差)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description='RMSD分析')
    parser.add_argument('-f', '--trajectory', required=True, help='轨迹文件 (xtc)')
    parser.add_argument('-s', '--topology', required=True, help='拓扑文件 (tpr)')
    parser.add_argument('-o', '--output', default='rmsd_analysis.png', help='输出图片')
    parser.add_argument('--protein_group', default='Backbone', help='蛋白选择组')
    parser.add_argument('--ligand_group', default='LIG', help='配体选择组')
    return parser.parse_args()


def calculate_rmsd(trajectory, topology, selection, output_xvg):
    """使用gmx rms计算RMSD"""
    cmd = f'gmx rms -s {topology} -f {trajectory} -o {output_xvg} -xvg none'
    
    # 创建输入文件
    input_str = f'{selection}\n{selection}\n'
    
    subprocess.run(cmd, shell=True, input=input_str, 
                  capture_output=True, text=True)
    
    # 读取结果
    data = np.loadtxt(output_xvg)
    return data[:, 0], data[:, 1]  # 时间, RMSD


def plot_rmsd(time, rmsd_protein, rmsd_ligand=None, output='rmsd.png'):
    """绘制RMSD图"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. 蛋白RMSD随时间变化
    ax1 = axes[0, 0]
    ax1.plot(time, rmsd_protein, 'b-', linewidth=0.8, alpha=0.7)
    ax1.axhline(y=np.mean(rmsd_protein), color='r', linestyle='--', 
                label=f'Mean: {np.mean(rmsd_protein):.2f} nm')
    ax1.fill_between(time, 
                     np.mean(rmsd_protein) - np.std(rmsd_protein),
                     np.mean(rmsd_protein) + np.std(rmsd_protein),
                     alpha=0.2, color='red')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('RMSD (nm)')
    ax1.set_title('Protein Backbone RMSD')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 蛋白RMSD分布
    ax2 = axes[0, 1]
    sns.histplot(rmsd_protein, kde=True, ax=ax2, color='blue')
    ax2.axvline(x=np.mean(rmsd_protein), color='r', linestyle='--',
                label=f'Mean: {np.mean(rmsd_protein):.2f} nm')
    ax2.set_xlabel('RMSD (nm)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Protein RMSD Distribution')
    ax2.legend()
    
    # 3. 配体RMSD (如果有)
    if rmsd_ligand is not None:
        ax3 = axes[1, 0]
        ax3.plot(time, rmsd_ligand, 'g-', linewidth=0.8, alpha=0.7)
        ax3.axhline(y=np.mean(rmsd_ligand), color='r', linestyle='--',
                    label=f'Mean: {np.mean(rmsd_ligand):.2f} nm')
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('RMSD (nm)')
        ax3.set_title('Ligand RMSD')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. 对比图
        ax4 = axes[1, 1]
        ax4.plot(time, rmsd_protein, 'b-', linewidth=0.8, alpha=0.6, label='Protein')
        ax4.plot(time, rmsd_ligand, 'g-', linewidth=0.8, alpha=0.6, label='Ligand')
        ax4.set_xlabel('Time (ps)')
        ax4.set_ylabel('RMSD (nm)')
        ax4.set_title('Protein vs Ligand RMSD')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    else:
        # 如果没有配体,显示滚动平均
        ax3 = axes[1, 0]
        window = min(100, len(rmsd_protein)//10)
        rolling_mean = np.convolve(rmsd_protein, np.ones(window)/window, mode='valid')
        ax3.plot(time[window-1:], rolling_mean, 'r-', linewidth=1.5, label=f'Rolling Mean ({window} frames)')
        ax3.plot(time, rmsd_protein, 'b-', linewidth=0.5, alpha=0.3)
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('RMSD (nm)')
        ax3.set_title('RMSD with Rolling Average')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. 统计信息
        ax4 = axes[1, 1]
        ax4.axis('off')
        stats_text = f"""
        RMSD Statistics:
        
        Protein Backbone:
          Mean: {np.mean(rmsd_protein):.3f} nm
          Std:  {np.std(rmsd_protein):.3f} nm
          Min:  {np.min(rmsd_protein):.3f} nm
          Max:  {np.max(rmsd_protein):.3f} nm
          
        Simulation Time:
          Total: {time[-1]/1000:.1f} ns
          Frames: {len(time)}
          
        Equilibration:
          First 10ns Mean: {np.mean(rmsd_protein[:len(rmsd_protein)//10]):.3f} nm
          Last 10ns Mean:  {np.mean(rmsd_protein[-len(rmsd_protein)//10:]):.3f} nm
        """
        ax4.text(0.1, 0.5, stats_text, fontsize=11, family='monospace',
                verticalalignment='center')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"RMSD分析图已保存: {output}")
    
    # 打印统计信息
    print("\nRMSD统计信息:")
    print(f"蛋白骨架:")
    print(f"  均值: {np.mean(rmsd_protein):.3f} ± {np.std(rmsd_protein):.3f} nm")
    print(f"  范围: {np.min(rmsd_protein):.3f} - {np.max(rmsd_protein):.3f} nm")
    if rmsd_ligand is not None:
        print(f"配体:")
        print(f"  均值: {np.mean(rmsd_ligand):.3f} ± {np.std(rmsd_ligand):.3f} nm")
        print(f"  范围: {np.min(rmsd_ligand):.3f} - {np.max(rmsd_ligand):.3f} nm")


def main():
    args = parse_args()
    
    print("RMSD分析")
    print("="*60)
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        # 计算蛋白RMSD
        print("\n计算蛋白骨架RMSD...")
        protein_xvg = os.path.join(tmpdir, 'protein_rmsd.xvg')
        time, rmsd_protein = calculate_rmsd(
            args.trajectory, args.topology, 
            args.protein_group, protein_xvg
        )
        
        # 计算配体RMSD (如果存在)
        rmsd_ligand = None
        try:
            print("计算配体RMSD...")
            ligand_xvg = os.path.join(tmpdir, 'ligand_rmsd.xvg')
            _, rmsd_ligand = calculate_rmsd(
                args.trajectory, args.topology,
                args.ligand_group, ligand_xvg
            )
        except:
            print("警告: 无法计算配体RMSD (配体可能不存在)")
        
        # 绘图
        print("\n生成图表...")
        plot_rmsd(time, rmsd_protein, rmsd_ligand, args.output)
    
    print("\nRMSD分析完成!")


if __name__ == '__main__':
    main()
