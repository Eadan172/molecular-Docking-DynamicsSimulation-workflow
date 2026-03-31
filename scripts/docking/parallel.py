#!/usr/bin/env python3
"""
并行分子对接脚本 (AutoDock Vina + multiprocessing)
适用于中等规模虚拟筛选

使用方法:
    python parallel.py --receptor protein.pdb --ligands ligands.sdf \\
                      --center 15.0 25.0 5.0 --box_size 60 60 60 \\
                      --exhaustiveness 32 --n_cpus 8 --output results/

性能加速:
    - 2核: ~1.8x
    - 4核: ~3.5x
    - 8核: ~6.5x
    - 16核: ~12x

计算时间估计 (exhaustiveness=32):
    - 10分子 + 8核: ~1-2分钟
    - 100分子 + 8核: ~10-15分钟
    - 1000分子 + 16核: ~1-2小时
"""

import argparse
import os
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

# 导入单核版本的函数
sys.path.insert(0, str(Path(__file__).parent))
from single_core import prepare_receptor, prepare_ligand_sdf, run_docking


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='并行分子对接 (AutoDock Vina)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 使用所有CPU核心
    python parallel.py -r protein.pdb -l ligands.sdf \\
                      -c 15.0 25.0 5.0 -s 60 60 60
    
    # 指定8核
    python parallel.py -r protein.pdb -l ligands.sdf \\
                      -c 15.0 25.0 5.0 -s 60 60 60 -n 8
    
    # 快速筛选
    python parallel.py -r protein.pdb -l ligands.sdf \\
                      -c 15.0 25.0 5.0 -s 60 60 60 \\
                      -e 16 -n 16 -o results_quick/
        """
    )
    
    parser.add_argument('-r', '--receptor', required=True,
                       help='受体蛋白PDB文件路径')
    parser.add_argument('-l', '--ligands', required=True,
                       help='配体分子SDF文件路径')
    parser.add_argument('-c', '--center', nargs=3, type=float, required=True,
                       metavar=('X', 'Y', 'Z'),
                       help='网格盒中心坐标')
    parser.add_argument('-s', '--box_size', nargs=3, type=float, 
                       default=[60, 60, 60], metavar=('SX', 'SY', 'SZ'),
                       help='网格盒大小 (默认: 60 60 60 Å)')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                       help='搜索详尽度 (默认: 32)')
    parser.add_argument('-n', '--n_cpus', type=int, default=None,
                       help='CPU核心数 (默认: 使用所有核心)')
    parser.add_argument('--num_modes', type=int, default=20,
                       help='输出构象数 (默认: 20)')
    parser.add_argument('--energy_range', type=float, default=3.0,
                       help='能量范围 (默认: 3.0 kcal/mol)')
    parser.add_argument('-o', '--output', default='results/docking_parallel',
                       help='输出目录')
    parser.add_argument('--keep_pdbqt', action='store_true',
                       help='保留中间PDBQT文件')
    parser.add_argument('--vina_exe', default='vina',
                       help='Vina可执行文件路径')
    
    return parser.parse_args()


def dock_single_worker(args):
    """
    单个对接任务的工作函数
    用于多进程池
    
    参数:
        args: (name, ligand_pdbqt, receptor_pdbqt, center, box_size, 
               exhaustiveness, num_modes, energy_range, output_dir, vina_exe)
    
    返回:
        result_dict: 包含对接结果的字典
    """
    (name, ligand_pdbqt, receptor_pdbqt, center, box_size,
     exhaustiveness, num_modes, energy_range, output_dir, vina_exe) = args
    
    output_file = os.path.join(output_dir, f'{name}_docked.pdbqt')
    log_file = os.path.join(output_dir, f'{name}.log')
    
    energy, time_taken = run_docking(
        receptor_pdbqt, ligand_pdbqt,
        center, box_size,
        exhaustiveness, num_modes, energy_range,
        output_file, log_file, vina_exe
    )
    
    return {
        'name': name,
        'energy': energy,
        'time': time_taken,
        'output': output_file,
        'log': log_file
    }


def main():
    """主函数"""
    args = parse_args()
    
    # 确定CPU核心数
    if args.n_cpus is None:
        args.n_cpus = cpu_count()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    print("="*70)
    print("并行分子对接 (AutoDock Vina)")
    print("="*70)
    print(f"受体: {args.receptor}")
    print(f"配体库: {args.ligands}")
    print(f"网格盒中心: {args.center}")
    print(f"网格盒大小: {args.box_size}")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"CPU核心数: {args.n_cpus}")
    print(f"输出目录: {args.output}")
    print("="*70)
    
    # 准备受体
    print("\n[1/3] 准备受体...")
    receptor_pdbqt = prepare_receptor(args.receptor, args.output)
    
    # 准备配体
    print("\n[2/3] 准备配体...")
    ligand_files = prepare_ligand_sdf(args.ligands, args.output)
    
    if not ligand_files:
        print("错误: 没有成功准备配体")
        sys.exit(1)
    
    # 准备任务参数
    print(f"\n[3/3] 运行并行对接 ({len(ligand_files)} 个分子, {args.n_cpus} 核)...")
    print("-"*70)
    
    task_args = [
        (name, ligand_pdbqt, receptor_pdbqt, args.center, args.box_size,
         args.exhaustiveness, args.num_modes, args.energy_range,
         args.output, args.vina_exe)
        for name, ligand_pdbqt in ligand_files
    ]
    
    # 运行并行对接
    total_start = time.time()
    
    with Pool(processes=args.n_cpus) as pool:
        results = pool.map(dock_single_worker, task_args)
    
    total_time = time.time() - total_start
    
    # 过滤失败的对接
    results = [r for r in results if r['energy'] is not None]
    
    # 排序
    results.sort(key=lambda x: x['energy'])
    
    # 输出结果
    print("\n" + "="*70)
    print("对接结果汇总")
    print("="*70)
    print(f"{'排名':<6}{'分子名称':<20}{'结合能(kcal/mol)':<20}{'耗时(秒)':<10}")
    print("-"*70)
    
    for i, r in enumerate(results, 1):
        print(f"{i:<6}{r['name']:<20}{r['energy']:<20.2f}{r['time']:<10.1f}")
    
    print("="*70)
    print(f"成功: {len(results)}/{len(ligand_files)} 个分子")
    print(f"总耗时: {total_time/60:.1f} 分钟")
    print(f"平均耗时: {total_time/len(ligand_files):.1f} 秒/分子")
    print(f"并行效率: {(sum(r['time'] for r in results)/total_time):.1f}x")
    print("="*70)
    
    # 保存CSV结果
    csv_file = os.path.join(args.output, 'docking_results.csv')
    with open(csv_file, 'w') as f:
        f.write('rank,name,binding_energy_kcal_mol,time_sec,output_file\n')
        for i, r in enumerate(results, 1):
            f.write(f"{i},{r['name']},{r['energy']:.2f},{r['time']:.1f},{r['output']}\n")
    
    print(f"\n结果已保存: {csv_file}")
    
    # 保存前3名分子用于MD模拟
    if len(results) >= 3:
        print("\n前3名分子:")
        for i, r in enumerate(results[:3], 1):
            print(f"  {i}. {r['name']}: {r['energy']:.2f} kcal/mol")
            # 复制到单独目录
            md_dir = os.path.join(args.output, 'top3_for_md')
            os.makedirs(md_dir, exist_ok=True)
            import shutil
            shutil.copy(r['output'], os.path.join(md_dir, f'top{i}_{r["name"]}.pdbqt'))
        print(f"\n已保存到: {md_dir}")
    
    # 清理临时文件
    if not args.keep_pdbqt:
        import shutil
        ligand_dir = os.path.join(args.output, 'ligands_pdbqt')
        if os.path.exists(ligand_dir):
            shutil.rmtree(ligand_dir)
            print(f"\n清理临时文件: {ligand_dir}")
    
    print("\n并行对接完成!")


if __name__ == '__main__':
    main()
