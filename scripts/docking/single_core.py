#!/usr/bin/env python3
"""
单核分子对接脚本 (AutoDock Vina)
适用于测试和小批量对接

使用方法:
    python single_core.py --receptor protein.pdb --ligands ligands.sdf \
                         --center 15.0 25.0 5.0 --box_size 60 60 60 \
                         --exhaustiveness 32 --output results/

计算时间估计:
    - 10分子: ~5-10分钟 (exhaustiveness=32)
    - 100分子: ~1-2小时
    - 1000分子: ~10-20小时

调参建议:
    - 快速测试: exhaustiveness=8-16
    - 标准筛选: exhaustiveness=32 (推荐)
    - 高精度: exhaustiveness=64-128
"""

import argparse
import os
import sys
import subprocess
import time
from pathlib import Path

# 添加父目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装,将使用外部工具转换格式")


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='单核分子对接 (AutoDock Vina)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 基本使用
    python single_core.py -r protein.pdb -l ligands.sdf \\
                         -c 15.0 25.0 5.0 -s 60 60 60 \\
                         -e 32 -o results/
    
    # 快速测试
    python single_core.py -r protein.pdb -l ligands.sdf \\
                         -c 15.0 25.0 5.0 -s 60 60 60 \\
                         -e 8 -o results_quick/
        """
    )
    
    parser.add_argument('-r', '--receptor', required=True,
                       help='受体蛋白PDB文件路径')
    parser.add_argument('-l', '--ligands', required=True,
                       help='配体分子SDF文件路径')
    parser.add_argument('-c', '--center', nargs=3, type=float, required=True,
                       metavar=('X', 'Y', 'Z'),
                       help='网格盒中心坐标 (例如: 15.0 25.0 5.0)')
    parser.add_argument('-s', '--box_size', nargs=3, type=float, 
                       default=[60, 60, 60], metavar=('SX', 'SY', 'SZ'),
                       help='网格盒大小 (默认: 60 60 60 Å)')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32,
                       help='搜索详尽度 (默认: 32, 推荐: 8-128)')
    parser.add_argument('-n', '--num_modes', type=int, default=20,
                       help='输出构象数 (默认: 20)')
    parser.add_argument('--energy_range', type=float, default=3.0,
                       help='能量范围 (默认: 3.0 kcal/mol)')
    parser.add_argument('-o', '--output', default='results/docking',
                       help='输出目录 (默认: results/docking)')
    parser.add_argument('--keep_pdbqt', action='store_true',
                       help='保留中间PDBQT文件')
    parser.add_argument('--vina_exe', default='vina',
                       help='Vina可执行文件路径 (默认: vina)')
    
    return parser.parse_args()


def prepare_receptor(pdb_file, output_dir):
    """
    准备受体PDBQT文件
    
    参数:
        pdb_file: 输入PDB文件
        output_dir: 输出目录
    
    返回:
        pdbqt_file: 输出PDBQT文件路径
    
    计算时间: ~5-10秒
    """
    pdbqt_file = os.path.join(output_dir, 'receptor.pdbqt')
    
    # 使用OpenBabel或MGLTools转换
    # 这里使用简单的命令行工具
    cmd = f'obabel {pdb_file} -O {pdbqt_file} -p 7.4'
    
    try:
        subprocess.run(cmd, shell=True, check=True, 
                      capture_output=True, timeout=30)
        print(f"受体准备完成: {pdbqt_file}")
        return pdbqt_file
    except subprocess.TimeoutExpired:
        print("错误: 受体准备超时")
        sys.exit(1)
    except subprocess.CalledProcessError:
        # 如果obabel失败,尝试使用prepare_receptor4.py
        cmd = f'prepare_receptor4.py -r {pdb_file} -o {pdbqt_file}'
        try:
            subprocess.run(cmd, shell=True, check=True, timeout=60)
            print(f"受体准备完成 (MGLTools): {pdbqt_file}")
            return pdbqt_file
        except:
            print("错误: 无法准备受体文件")
            print("请确保安装OpenBabel或MGLTools")
            sys.exit(1)


def prepare_ligand_sdf(sdf_file, output_dir):
    """
    将SDF文件中的多个分子转换为单独的PDBQT文件
    
    参数:
        sdf_file: 输入SDF文件
        output_dir: 输出目录
    
    返回:
        ligand_files: PDBQT文件列表
    
    计算时间: ~1-2秒/分子
    """
    ligand_dir = os.path.join(output_dir, 'ligands_pdbqt')
    os.makedirs(ligand_dir, exist_ok=True)
    
    ligand_files = []
    
    if RDKIT_AVAILABLE:
        # 使用RDKit读取SDF
        supplier = Chem.SDMolSupplier(sdf_file)
        
        for i, mol in enumerate(supplier):
            if mol is None:
                print(f"警告: 无法读取分子 {i+1}")
                continue
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 生成3D坐标
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            
            # 保存为PDBQT
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'ligand_{i+1}'
            pdbqt_file = os.path.join(ligand_dir, f'{name}.pdbqt')
            
            # 使用OpenBabel转换
            pdb_file = os.path.join(ligand_dir, f'{name}.pdb')
            Chem.MolToPDBFile(mol, pdb_file)
            
            cmd = f'obabel {pdb_file} -O {pdbqt_file}'
            subprocess.run(cmd, shell=True, capture_output=True)
            
            if os.path.exists(pdbqt_file):
                ligand_files.append((name, pdbqt_file))
                os.remove(pdb_file)  # 清理临时文件
    else:
        # 使用OpenBabel直接转换
        cmd = f'obabel {sdf_file} -O {ligand_dir}/ligand_.pdbqt -m'
        subprocess.run(cmd, shell=True, check=True)
        
        # 收集生成的文件
        for f in os.listdir(ligand_dir):
            if f.endswith('.pdbqt'):
                name = f.replace('.pdbqt', '')
                ligand_files.append((name, os.path.join(ligand_dir, f)))
    
    print(f"配体准备完成: {len(ligand_files)} 个分子")
    return ligand_files


def run_docking(receptor_pdbqt, ligand_pdbqt, center, box_size, 
                exhaustiveness, num_modes, energy_range, 
                output_file, log_file, vina_exe='vina'):
    """
    运行单个分子对接
    
    参数:
        receptor_pdbqt: 受体PDBQT文件
        ligand_pdbqt: 配体PDBQT文件
        center: 网格盒中心 [x, y, z]
        box_size: 网格盒大小 [sx, sy, sz]
        exhaustiveness: 搜索详尽度
        num_modes: 输出构象数
        energy_range: 能量范围
        output_file: 输出文件
        log_file: 日志文件
        vina_exe: Vina可执行文件
    
    返回:
        best_energy: 最佳结合能
        time_taken: 计算时间
    
    计算时间: 
        - exhaustiveness=8: ~10s
        - exhaustiveness=32: ~40s
        - exhaustiveness=128: ~160s
    """
    start_time = time.time()
    
    cmd = [
        vina_exe,
        '--receptor', receptor_pdbqt,
        '--ligand', ligand_pdbqt,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(box_size[0]),
        '--size_y', str(box_size[1]),
        '--size_z', str(box_size[2]),
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
        '--energy_range', str(energy_range),
        '--out', output_file,
        '--log', log_file
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            print(f"对接失败: {result.stderr}")
            return None, None
        
        # 解析日志文件获取最佳能量
        best_energy = None
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if line.startswith('REMARK VINA RESULT:'):
                        parts = line.split()
                        if len(parts) >= 4:
                            best_energy = float(parts[3])
                        break
        
        time_taken = time.time() - start_time
        return best_energy, time_taken
        
    except subprocess.TimeoutExpired:
        print("对接超时")
        return None, None
    except Exception as e:
        print(f"对接错误: {e}")
        return None, None


def main():
    """主函数"""
    args = parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    print("="*60)
    print("单核分子对接 (AutoDock Vina)")
    print("="*60)
    print(f"受体: {args.receptor}")
    print(f"配体库: {args.ligands}")
    print(f"网格盒中心: {args.center}")
    print(f"网格盒大小: {args.box_size}")
    print(f"Exhaustiveness: {args.exhaustiveness}")
    print(f"输出目录: {args.output}")
    print("="*60)
    
    # 准备受体
    print("\n[1/3] 准备受体...")
    receptor_pdbqt = prepare_receptor(args.receptor, args.output)
    
    # 准备配体
    print("\n[2/3] 准备配体...")
    ligand_files = prepare_ligand_sdf(args.ligands, args.output)
    
    if not ligand_files:
        print("错误: 没有成功准备配体")
        sys.exit(1)
    
    # 运行对接
    print(f"\n[3/3] 运行对接 ({len(ligand_files)} 个分子)...")
    print("-"*60)
    
    results = []
    total_start = time.time()
    
    for i, (name, ligand_pdbqt) in enumerate(ligand_files, 1):
        print(f"\n[{i}/{len(ligand_files)}] 对接: {name}")
        
        output_file = os.path.join(args.output, f'{name}_docked.pdbqt')
        log_file = os.path.join(args.output, f'{name}.log')
        
        energy, time_taken = run_docking(
            receptor_pdbqt, ligand_pdbqt,
            args.center, args.box_size,
            args.exhaustiveness, args.num_modes, args.energy_range,
            output_file, log_file, args.vina_exe
        )
        
        if energy is not None:
            print(f"  结合能: {energy:.2f} kcal/mol")
            print(f"  耗时: {time_taken:.1f} 秒")
            results.append({
                'name': name,
                'energy': energy,
                'time': time_taken,
                'output': output_file,
                'log': log_file
            })
        else:
            print(f"  失败")
    
    total_time = time.time() - total_start
    
    # 排序并输出结果
    results.sort(key=lambda x: x['energy'])
    
    print("\n" + "="*60)
    print("对接结果汇总")
    print("="*60)
    print(f"{'排名':<6}{'分子名称':<20}{'结合能(kcal/mol)':<20}{'耗时(秒)':<10}")
    print("-"*60)
    
    for i, r in enumerate(results, 1):
        print(f"{i:<6}{r['name']:<20}{r['energy']:<20.2f}{r['time']:<10.1f}")
    
    print("="*60)
    print(f"总计: {len(results)} 个分子")
    print(f"总耗时: {total_time/60:.1f} 分钟")
    print(f"平均耗时: {total_time/len(results):.1f} 秒/分子")
    print("="*60)
    
    # 保存CSV结果
    csv_file = os.path.join(args.output, 'docking_results.csv')
    with open(csv_file, 'w') as f:
        f.write('rank,name,binding_energy_kcal_mol,time_sec,output_file\n')
        for i, r in enumerate(results, 1):
            f.write(f"{i},{r['name']},{r['energy']:.2f},{r['time']:.1f},{r['output']}\n")
    
    print(f"\n结果已保存: {csv_file}")
    
    # 清理临时文件
    if not args.keep_pdbqt:
        import shutil
        ligand_dir = os.path.join(args.output, 'ligands_pdbqt')
        if os.path.exists(ligand_dir):
            shutil.rmtree(ligand_dir)
            print(f"清理临时文件: {ligand_dir}")
    
    print("\n对接完成!")


if __name__ == '__main__':
    main()
