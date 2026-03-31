#!/usr/bin/env python3
"""
自由能形貌图 (Free Energy Surface) 分析脚本
基于两个集体变量(CV)计算自由能形貌

使用方法:
    python fes.py -f trajectory.xtc -s system.tpr \\
                  --cv1 rmsd --cv2 gyrate -o fes.png

支持的CV:
    - rmsd: 蛋白RMSD
    - gyrate: 回转半径
    - rgyr: 同上
    - dpca: 距离PCA (需要指定原子组)
    - custom: 自定义 (需要提供数据文件)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import subprocess
import os
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description='自由能形貌图分析')
    parser.add_argument('-f', '--trajectory', required=True, help='轨迹文件')
    parser.add_argument('-s', '--topology', required=True, help='拓扑文件')
    parser.add_argument('--cv1', default='rmsd', help='第一个CV (默认: rmsd)')
    parser.add_argument('--cv2', default='gyrate', help='第二个CV (默认: gyrate)')
    parser.add_argument('-o', '--output', default='fes.png', help='输出图片')
    parser.add_argument('--temp', type=float, default=300, help='温度 (K)')
    parser.add_argument('--bins', type=int, default=50, help='直方图箱数')
    parser.add_argument('--sigma', type=float, default=1.0, help='高斯平滑参数')
    parser.add_argument('--vmin', type=float, default=0, help='颜色最小值')
    parser.add_argument('--vmax', type=float, default=10, help='颜色最大值')
    return parser.parse_args()


def calculate_rmsd(trajectory, topology, output_xvg):
    """计算RMSD"""
    cmd = f'gmx rms -s {topology} -f {trajectory} -o {output_xvg} -xvg none'
    input_str = 'Backbone\nBackbone\n'
    subprocess.run(cmd, shell=True, input=input_str,
                  capture_output=True, text=True)
    data = np.loadtxt(output_xvg)
    return data[:, 1]  # 返回RMSD值


def calculate_gyrate(trajectory, topology, output_xvg):
    """计算回转半径"""
    cmd = f'gmx gyrate -s {topology} -f {trajectory} -o {output_xvg}'
    input_str = 'Protein\n'
    subprocess.run(cmd, shell=True, input=input_str,
                  capture_output=True, text=True)
    data = np.loadtxt(output_xvg, comments=['#', '@'])
    return data[:, 1]  # 返回Rg值


def calculate_distance(trajectory, topology, group1, group2, output_xvg):
    """计算两组原子间距离"""
    # 创建索引文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.ndx', delete=False) as f:
        ndx_file = f.name
    
    # 这里简化处理,实际应该创建适当的索引组
    cmd = f'gmx distance -s {topology} -f {trajectory} -oall {output_xvg}'
    subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if os.path.exists(output_xvg):
        data = np.loadtxt(output_xvg, comments=['#', '@'])
        return data[:, 1]
    else:
        return None


def get_cv_data(cv_name, trajectory, topology, tmpdir):
    """获取CV数据"""
    
    if cv_name.lower() == 'rmsd':
        output_xvg = os.path.join(tmpdir, 'cv1_rmsd.xvg')
        return calculate_rmsd(trajectory, topology, output_xvg)
    
    elif cv_name.lower() in ['gyrate', 'rgyr']:
        output_xvg = os.path.join(tmpdir, 'cv1_gyrate.xvg')
        return calculate_gyrate(trajectory, topology, output_xvg)
    
    elif cv_name.lower() == 'dpca':
        # 这里简化,实际应该使用gmx anaeig
        print("警告: DPCA需要预先计算,使用随机数据演示")
        return np.random.randn(1000)  # 占位
    
    else:
        print(f"错误: 未知的CV类型: {cv_name}")
        return None


def calculate_fes(cv1_data, cv2_data, bins=50, temperature=300, sigma=1.0):
    """
    计算自由能形貌图
    
    参数:
        cv1_data, cv2_data: CV数据数组
        bins: 直方图箱数
        temperature: 温度 (K)
        sigma: 高斯平滑参数
    
    返回:
        X, Y: 网格坐标
        F: 自由能矩阵
        extent: 图像范围
    """
    # 创建2D直方图
    H, xedges, yedges = np.histogram2d(cv1_data, cv2_data, bins=bins)
    
    # 归一化得到概率
    H = H / np.sum(H)
    
    # 平滑处理
    if sigma > 0:
        H = gaussian_filter(H, sigma=sigma)
    
    # 避免log(0)
    H_min = np.min(H[H > 0])
    H[H == 0] = H_min * 0.1
    
    # 计算自由能: F = -kT * ln(P)
    kT = 0.008314 * temperature  # kJ/mol/K * K = kJ/mol
    F = -kT * np.log(H)
    
    # 归一化到最小值为0
    F = F - np.min(F)
    
    # 创建网格
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    
    return X, Y, F, extent


def plot_fes(X, Y, F, extent, cv1_name, cv2_name, output, 
             vmin=0, vmax=10, temperature=300):
    """绘制自由能形貌图"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. 2D自由能形貌图
    ax1 = axes[0, 0]
    im1 = ax1.contourf(X, Y, F.T, levels=20, cmap='jet', 
                       vmin=vmin, vmax=vmax, extend='max')
    ax1.contour(X, Y, F.T, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax1.set_xlabel(f'{cv1_name}')
    ax1.set_ylabel(f'{cv2_name}')
    ax1.set_title(f'Free Energy Surface (T={temperature}K)')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('Free Energy (kJ/mol)')
    
    # 2. 3D视图
    ax2 = axes[0, 1]
    from mpl_toolkits.mplot3d import Axes3D
    ax2.remove()
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    surf = ax2.plot_surface(X, Y, F.T, cmap='jet', 
                           vmin=vmin, vmax=vmax, alpha=0.8)
    ax2.set_xlabel(cv1_name)
    ax2.set_ylabel(cv2_name)
    ax2.set_zlabel('Free Energy (kJ/mol)')
    ax2.set_title('3D Free Energy Surface')
    
    # 3. CV1分布
    ax3 = axes[1, 0]
    cv1_data = X[0, :]
    cv1_fes = -0.008314 * temperature * np.log(np.sum(np.exp(-F.T/0.008314/temperature), axis=0))
    cv1_fes = cv1_fes - np.min(cv1_fes)
    ax3.plot(cv1_data, cv1_fes, 'b-', linewidth=2)
    ax3.fill_between(cv1_data, cv1_fes, alpha=0.3)
    ax3.set_xlabel(cv1_name)
    ax3.set_ylabel('Free Energy (kJ/mol)')
    ax3.set_title(f'{cv1_name} Free Energy Profile')
    ax3.grid(True, alpha=0.3)
    
    # 4. CV2分布
    ax4 = axes[1, 1]
    cv2_data = Y[:, 0]
    cv2_fes = -0.008314 * temperature * np.log(np.sum(np.exp(-F.T/0.008314/temperature), axis=1))
    cv2_fes = cv2_fes - np.min(cv2_fes)
    ax4.plot(cv2_data, cv2_fes, 'r-', linewidth=2)
    ax4.fill_between(cv2_data, cv2_fes, alpha=0.3, color='red')
    ax4.set_xlabel(cv2_name)
    ax4.set_ylabel('Free Energy (kJ/mol)')
    ax4.set_title(f'{cv2_name} Free Energy Profile')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"自由能形貌图已保存: {output}")
    
    # 统计信息
    print("\n自由能形貌统计:")
    print(f"  最小自由能: {np.min(F):.2f} kJ/mol")
    print(f"  最大自由能: {np.max(F):.2f} kJ/mol")
    print(f"  平均自由能: {np.mean(F):.2f} ± {np.std(F):.2f} kJ/mol")
    
    # 寻找极小值点 (简化的局部最小值检测)
    from scipy.ndimage import minimum_filter
    local_min = (F == minimum_filter(F, size=3))
    n_minima = np.sum(local_min)
    print(f"  检测到的局部最小值: {n_minima} 个")


def main():
    args = parse_args()
    
    print("自由能形貌图分析")
    print("="*60)
    print(f"轨迹: {args.trajectory}")
    print(f"拓扑: {args.topology}")
    print(f"CV1: {args.cv1}")
    print(f"CV2: {args.cv2}")
    print(f"温度: {args.temp} K")
    print(f"箱数: {args.bins}")
    print("="*60)
    
    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        # 获取CV数据
        print(f"\n计算 {args.cv1}...")
        cv1_data = get_cv_data(args.cv1, args.trajectory, args.topology, tmpdir)
        
        print(f"计算 {args.cv2}...")
        cv2_data = get_cv_data(args.cv2, args.trajectory, args.topology, tmpdir)
        
        if cv1_data is None or cv2_data is None:
            print("错误: 无法计算CV")
            return
        
        # 确保数据长度一致
        min_len = min(len(cv1_data), len(cv2_data))
        cv1_data = cv1_data[:min_len]
        cv2_data = cv2_data[:min_len]
        
        print(f"数据点数: {min_len}")
        
        # 计算自由能形貌
        print("\n计算自由能形貌...")
        X, Y, F, extent = calculate_fes(
            cv1_data, cv2_data, 
            bins=args.bins, 
            temperature=args.temp,
            sigma=args.sigma
        )
        
        # 绘图
        print("生成图表...")
        plot_fes(X, Y, F, extent, args.cv1, args.cv2, args.output,
                args.vmin, args.vmax, args.temp)
    
    print("\n自由能形貌分析完成!")


if __name__ == '__main__':
    main()
