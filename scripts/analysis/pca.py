#!/usr/bin/env python3
"""
PCA (主成分分析) 脚本
分析MD轨迹的主要运动模式

使用方法:
    python pca.py -f trajectory.xtc -s system.tpr -o pca_analysis.png

输出:
    - 投影到PC1-PC2的散点图 (颜色表示时间)
    - 特征值分布图 (解释方差比例)
    - 自由能形貌图 (基于PC1-PC2)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import subprocess
import os
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description='PCA分析')
    parser.add_argument('-f', '--trajectory', required=True, help='轨迹文件')
    parser.add_argument('-s', '--topology', required=True, help='拓扑文件')
    parser.add_argument('-o', '--output', default='pca_analysis.png', help='输出图片')
    parser.add_argument('--selection', default='Backbone', help='原子选择')
    parser.add_argument('--n_components', type=int, default=10, help='主成分数')
    return parser.parse_args()


def extract_coordinates(trajectory, topology, selection, output_gro):
    """提取选定原子的坐标"""
    # 使用gmx trjconv提取坐标
    cmd = f'gmx trjconv -s {topology} -f {trajectory} -o {output_gro} -fit rot+trans'
    
    input_str = f'{selection}\n'
    subprocess.run(cmd, shell=True, input=input_str,
                  capture_output=True, text=True)


def read_gro_coordinates(gro_file):
    """读取gro文件中的坐标"""
    coords = []
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        n_atoms = int(lines[1])
        
        for line in lines[2:2+n_atoms]:
            # gro格式: resnum, resname, atomname, atomnum, x, y, z
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            coords.extend([x, y, z])
    
    return np.array(coords)


def perform_pca(trajectory, topology, selection, n_components=10):
    """执行PCA分析"""
    
    print("提取轨迹坐标...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # 提取每一帧的坐标
        # 使用gmx covar计算协方差矩阵和特征向量
        
        eigenvec_file = os.path.join(tmpdir, 'eigenvec.trr')
        eigenval_file = os.path.join(tmpdir, 'eigenval.xvg')
        
        # 计算协方差矩阵
        cmd = f'gmx covar -s {topology} -f {trajectory} -o {eigenval_file} -v {eigenvec_file}'
        input_str = f'{selection}\n{selection}\n'
        subprocess.run(cmd, shell=True, input=input_str,
                      capture_output=True, text=True)
        
        # 读取特征值
        eigenvalues = np.loadtxt(eigenval_file, comments=['#', '@'])[:, 1]
        
        # 投影轨迹到主成分
        proj_file = os.path.join(tmpdir, 'proj.xvg')
        cmd = f'gmx anaeig -s {topology} -f {trajectory} -v {eigenvec_file} -proj {proj_file} -first 1 -last {n_components}'
        subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # 读取投影
        projections = np.loadtxt(proj_file, comments=['#', '@'])
        
    return eigenvalues, projections


def plot_pca(eigenvalues, projections, output='pca.png'):
    """绘制PCA结果"""
    
    fig = plt.figure(figsize=(16, 12))
    
    # 1. 特征值分布 (解释方差)
    ax1 = plt.subplot(2, 3, 1)
    total_var = np.sum(eigenvalues)
    explained_var = eigenvalues / total_var * 100
    
    ax1.bar(range(1, len(explained_var)+1), explained_var, alpha=0.7)
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Explained Variance (%)')
    ax1.set_title('PCA Explained Variance')
    ax1.grid(True, alpha=0.3)
    
    # 累积方差
    ax1b = ax1.twinx()
    cumulative_var = np.cumsum(explained_var)
    ax1b.plot(range(1, len(cumulative_var)+1), cumulative_var, 'r-o', markersize=4)
    ax1b.set_ylabel('Cumulative Variance (%)', color='r')
    ax1b.tick_params(axis='y', labelcolor='r')
    
    # 2. PC1-PC2投影 (颜色=时间)
    ax2 = plt.subplot(2, 3, 2)
    time_colors = np.arange(len(projections))
    scatter = ax2.scatter(projections[:, 0], projections[:, 1], 
                         c=time_colors, cmap='viridis', s=5, alpha=0.6)
    ax2.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)')
    ax2.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)')
    ax2.set_title('Trajectory Projection (PC1-PC2)')
    plt.colorbar(scatter, ax=ax2, label='Frame')
    ax2.grid(True, alpha=0.3)
    
    # 3. PC1-PC3投影
    ax3 = plt.subplot(2, 3, 3)
    scatter3 = ax3.scatter(projections[:, 0], projections[:, 2],
                          c=time_colors, cmap='plasma', s=5, alpha=0.6)
    ax3.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)')
    ax3.set_ylabel(f'PC3 ({explained_var[2]:.1f}%)')
    ax3.set_title('Trajectory Projection (PC1-PC3)')
    plt.colorbar(scatter3, ax=ax3, label='Frame')
    ax3.grid(True, alpha=0.3)
    
    # 4. 自由能形貌图 (PC1-PC2)
    ax4 = plt.subplot(2, 3, 4)
    
    # 2D直方图
    H, xedges, yedges = np.histogram2d(projections[:, 0], projections[:, 1], bins=50)
    
    # 计算自由能 (F = -kT * ln(P))
    H = H.T  # 转置以匹配imshow
    H = H / np.sum(H)  # 归一化
    H[H == 0] = np.min(H[H > 0])  # 避免log(0)
    F = -np.log(H)
    F = F - np.min(F)  # 归一化到0
    
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    im = ax4.imshow(F, extent=extent, origin='lower', cmap='jet', aspect='auto')
    ax4.set_xlabel(f'PC1 ({explained_var[0]:.1f}%)')
    ax4.set_ylabel(f'PC2 ({explained_var[1]:.1f}%)')
    ax4.set_title('Free Energy Landscape (PC1-PC2)')
    plt.colorbar(im, ax=ax4, label='Free Energy (kT)')
    
    # 5. PC1随时间变化
    ax5 = plt.subplot(2, 3, 5)
    time_ps = np.arange(len(projections)) * 10  # 假设每10ps一帧
    ax5.plot(time_ps, projections[:, 0], 'b-', linewidth=0.8, alpha=0.7)
    ax5.set_xlabel('Time (ps)')
    ax5.set_ylabel('PC1 Projection')
    ax5.set_title('PC1 vs Time')
    ax5.grid(True, alpha=0.3)
    
    # 6. PC2随时间变化
    ax6 = plt.subplot(2, 3, 6)
    ax6.plot(time_ps, projections[:, 1], 'r-', linewidth=0.8, alpha=0.7)
    ax6.set_xlabel('Time (ps)')
    ax6.set_ylabel('PC2 Projection')
    ax6.set_title('PC2 vs Time')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"PCA分析图已保存: {output}")
    
    # 打印统计信息
    print("\nPCA统计信息:")
    print(f"总方差解释:")
    for i in range(min(5, len(explained_var))):
        print(f"  PC{i+1}: {explained_var[i]:.2f}%")
    print(f"  前5PC累积: {np.sum(explained_var[:5]):.2f}%")
    print(f"  前10PC累积: {np.sum(explained_var[:10]):.2f}%")


def main():
    args = parse_args()
    
    print("PCA分析")
    print("="*60)
    print(f"轨迹: {args.trajectory}")
    print(f"拓扑: {args.topology}")
    print(f"选择: {args.selection}")
    print(f"主成分数: {args.n_components}")
    print("="*60)
    
    # 执行PCA
    eigenvalues, projections = perform_pca(
        args.trajectory, args.topology,
        args.selection, args.n_components
    )
    
    # 绘图
    print("\n生成图表...")
    plot_pca(eigenvalues, projections, args.output)
    
    print("\nPCA分析完成!")


if __name__ == '__main__':
    main()
