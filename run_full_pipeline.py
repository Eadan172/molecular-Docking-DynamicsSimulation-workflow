#!/usr/bin/env python3
"""
完整分子模拟流程脚本
1. AutoDock Vina分子对接
2. Meeko+RDKit相互作用分析和可视化
3. Gromacs分子动力学模拟准备
4. 动力学结果分析（RMSD、PCA、FES）

使用方法:
    python run_full_pipeline.py --protein data/protein.pdb --ligand data/ligands.sdf \\
                               --center 24 22 17 --box_size 20 20 20
"""

import argparse
import os
import sys
import subprocess
import time
import shutil
from pathlib import Path
from typing import Dict, List, Optional


def check_environment() -> Dict[str, bool]:
    """检查环境"""
    print("="*70)
    print("环境检查")
    print("="*70)
    
    tools = {
        'vina': False,
        'gromacs': False,
        'obabel': False,
        'python_packages': False
    }
    
    # 检查AutoDock Vina
    try:
        result = subprocess.run(['vina', '--version'], capture_output=True, text=True, timeout=10)
        if result.returncode == 0 or 'AutoDock Vina' in result.stdout + result.stderr:
            tools['vina'] = True
            print("✓ AutoDock Vina 可用")
    except:
        print("✗ AutoDock Vina 未找到")
    
    # 检查GROMACS
    try:
        result = subprocess.run(['gmx', '--version'], capture_output=True, text=True, timeout=10)
        if result.returncode == 0 or 'GROMACS' in result.stdout + result.stderr:
            tools['gromacs'] = True
            print("✓ GROMACS 可用")
    except:
        print("✗ GROMACS 未找到")
    
    # 检查OpenBabel
    try:
        result = subprocess.run(['obabel', '--version'], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            tools['obabel'] = True
            print("✓ OpenBabel 可用")
    except:
        print("✗ OpenBabel 未找到")
    
    # 检查Python包
    required_packages = ['meeko', 'rdkit', 'numpy', 'matplotlib', 'seaborn', 'scikit-learn']
    missing_packages = []
    for pkg in required_packages:
        try:
            __import__(pkg)
        except ImportError:
            missing_packages.append(pkg)
    
    if not missing_packages:
        tools['python_packages'] = True
        print("✓ 所有必需Python包已安装")
    else:
        print(f"✗ 缺少Python包: {', '.join(missing_packages)}")
    
    print("="*70)
    return tools


def run_docking(protein_file: str, ligand_file: str, center: List[float], 
               box_size: List[float], output_dir: str, exhaustiveness: int = 8) -> bool:
    """运行分子对接"""
    print("\n" + "="*70)
    print("步骤 1/4: 分子对接 (AutoDock Vina)")
    print("="*70)
    
    script_path = Path(__file__).parent / 'scripts' / 'docking' / 'single_core.py'
    
    cmd = [
        sys.executable, str(script_path),
        '-r', protein_file,
        '-l', ligand_file,
        '-c', str(center[0]), str(center[1]), str(center[2]),
        '-s', str(box_size[0]), str(box_size[1]), str(box_size[2]),
        '-e', str(exhaustiveness),
        '-o', os.path.join(output_dir, 'docking')
    ]
    
    print(f"运行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("✓ 分子对接完成")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ 分子对接失败: {e}")
        return False


def run_interaction_analysis(output_dir: str) -> bool:
    """运行相互作用分析和可视化"""
    print("\n" + "="*70)
    print("步骤 2/4: 相互作用分析 (Meeko+RDKit)")
    print("="*70)
    
    script_path = Path(__file__).parent / 'scripts' / 'visualization' / 'visualize_interactions.py'
    
    cmd = [sys.executable, str(script_path)]
    
    print(f"运行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, cwd=str(Path(__file__).parent))
        print("✓ 相互作用分析完成")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ 相互作用分析失败: {e}")
        return False


def prepare_md_system(protein_file: str, ligand_file: str, output_dir: str) -> bool:
    """准备Gromacs MD系统"""
    print("\n" + "="*70)
    print("步骤 3/4: 准备Gromacs MD系统")
    print("="*70)
    
    script_path = Path(__file__).parent / 'scripts' / 'md' / 'setup.sh'
    
    if not script_path.exists():
        print("警告: MD设置脚本不存在，跳过此步骤")
        print("提示: 在Linux/Mac环境下可以运行完整MD流程")
        return False
    
    cmd = [
        'bash', str(script_path),
        '-p', protein_file,
        '-l', ligand_file,
        '-o', os.path.join(output_dir, 'md')
    ]
    
    print(f"运行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("✓ MD系统准备完成")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ MD系统准备失败: {e}")
        print("提示: 确保在Linux/Mac环境或WSL中运行")
        return False


def create_batch_scripts(output_dir: str) -> None:
    """创建批处理脚本"""
    print("\n" + "="*70)
    print("创建批处理脚本")
    print("="*70)
    
    # 创建Windows批处理脚本
    bat_script = Path(output_dir) / 'run_docking.bat'
    with open(bat_script, 'w') as f:
        f.write('@echo off\n')
        f.write('echo 运行分子对接...\n')
        f.write(f'cd /d "{Path(__file__).parent}"\n')
        f.write('python run_full_pipeline.py --protein data/protein.pdb --ligand data/ligands.sdf --center 24 22 17 --box_size 20 20 20\n')
        f.write('pause\n')
    
    # 创建Linux Shell脚本
    sh_script = Path(output_dir) / 'run_docking.sh'
    with open(sh_script, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('echo "运行分子对接..."\n')
        f.write(f'cd "{Path(__file__).parent}"\n')
        f.write('python run_full_pipeline.py --protein data/protein.pdb --ligand data/ligands.sdf --center 24 22 17 --box_size 20 20 20\n')
    
    os.chmod(sh_script, 0o755)
    
    # 创建集群提交脚本
    cluster_script = Path(output_dir) / 'submit_md.sh'
    with open(cluster_script, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --job-name=md_simulation\n')
        f.write('#SBATCH --partition=gpu\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --ntasks-per-node=8\n')
        f.write('#SBATCH --gres=gpu:1\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write('#SBATCH --output=md_%j.out\n')
        f.write('#SBATCH --error=md_%j.err\n\n')
        f.write('echo "Job started at $(date)"\n')
        f.write('cd results/md\n')
        f.write('gmx grompp -f ../../params/nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 2\n')
        f.write('gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp 8\n')
        f.write('gmx grompp -f ../../params/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 2\n')
        f.write('gmx mdrun -v -deffnm npt -ntmpi 1 -ntomp 8\n')
        f.write('gmx grompp -f ../../params/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 2\n')
        f.write('gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 8\n')
        f.write('echo "Job finished at $(date)"\n')
    
    os.chmod(cluster_script, 0o755)
    
    print(f"✓ Windows批处理脚本: {bat_script}")
    print(f"✓ Linux Shell脚本: {sh_script}")
    print(f"✓ 集群提交脚本: {cluster_script}")


def create_analysis_scripts(output_dir: str) -> None:
    """创建分析脚本"""
    print("\n" + "="*70)
    print("创建分析脚本")
    print("="*70)
    
    # RMSD分析脚本
    rmsd_script = Path(output_dir) / 'analyze_rmsd.py'
    with open(rmsd_script, 'w') as f:
        f.write('#!/usr/bin/env python3\n')
        f.write('"""RMSD分析脚本"""\n')
        f.write('import sys\n')
        f.write('from pathlib import Path\n')
        f.write('sys.path.insert(0, str(Path(__file__).parent.parent))\n')
        f.write('from scripts.analysis.rmsd import main\n')
        f.write('if __name__ == "__main__":\n')
        f.write('    main()\n')
    
    # PCA分析脚本
    pca_script = Path(output_dir) / 'analyze_pca.py'
    with open(pca_script, 'w') as f:
        f.write('#!/usr/bin/env python3\n')
        f.write('"""PCA分析脚本"""\n')
        f.write('import sys\n')
        f.write('from pathlib import Path\n')
        f.write('sys.path.insert(0, str(Path(__file__).parent.parent))\n')
        f.write('from scripts.analysis.pca import main\n')
        f.write('if __name__ == "__main__":\n')
        f.write('    main()\n')
    
    # FES分析脚本
    fes_script = Path(output_dir) / 'analyze_fes.py'
    with open(fes_script, 'w') as f:
        f.write('#!/usr/bin/env python3\n')
        f.write('"""自由能形貌图分析脚本"""\n')
        f.write('import sys\n')
        f.write('from pathlib import Path\n')
        f.write('sys.path.insert(0, str(Path(__file__).parent.parent))\n')
        f.write('from scripts.analysis.fes import main\n')
        f.write('if __name__ == "__main__":\n')
        f.write('    main()\n')
    
    os.chmod(rmsd_script, 0o755)
    os.chmod(pca_script, 0o755)
    os.chmod(fes_script, 0o755)
    
    print(f"✓ RMSD分析脚本: {rmsd_script}")
    print(f"✓ PCA分析脚本: {pca_script}")
    print(f"✓ FES分析脚本: {fes_script}")


def main():
    parser = argparse.ArgumentParser(
        description='完整分子模拟流程',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 基本使用（使用示例数据）
    python run_full_pipeline.py
    
    # 自定义参数
    python run_full_pipeline.py --protein data/protein.pdb \\
                               --ligand data/ligands.sdf \\
                               --center 24 22 17 \\
                               --box_size 20 20 20
    
    # 快速测试
    python run_full_pipeline.py --exhaustiveness 4
        """
    )
    
    parser.add_argument('--protein', default='data/protein.pdb',
                       help='受体蛋白PDB文件 (默认: data/protein.pdb)')
    parser.add_argument('--ligand', default='data/ligands.sdf',
                       help='配体SDF文件 (默认: data/ligands.sdf)')
    parser.add_argument('--center', nargs=3, type=float, default=[24.0, 22.0, 17.0],
                       metavar=('X', 'Y', 'Z'),
                       help='网格盒中心坐标 (默认: 24.0 22.0 17.0)')
    parser.add_argument('--box_size', nargs=3, type=float, default=[20.0, 20.0, 20.0],
                       metavar=('SX', 'SY', 'SZ'),
                       help='网格盒大小 (默认: 20.0 20.0 20.0 Å)')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                       help='搜索详尽度 (默认: 8, 推荐: 8-32)')
    parser.add_argument('--output_dir', default='results',
                       help='输出目录 (默认: results)')
    parser.add_argument('--skip_docking', action='store_true',
                       help='跳过对接步骤')
    parser.add_argument('--skip_interaction', action='store_true',
                       help='跳过相互作用分析')
    parser.add_argument('--skip_md', action='store_true',
                       help='跳过MD准备')
    
    args = parser.parse_args()
    
    print("="*70)
    print("完整分子模拟流程")
    print("="*70)
    print(f"蛋白文件: {args.protein}")
    print(f"配体文件: {args.ligand}")
    print(f"网格中心: {args.center}")
    print(f"网格大小: {args.box_size}")
    print(f"搜索详尽度: {args.exhaustiveness}")
    print(f"输出目录: {args.output_dir}")
    print("="*70)
    
    # 检查环境
    env = check_environment()
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # 步骤1: 分子对接
    docking_success = True
    if not args.skip_docking and env['vina']:
        docking_success = run_docking(
            args.protein, args.ligand,
            args.center, args.box_size,
            args.output_dir, args.exhaustiveness
        )
    
    # 步骤2: 相互作用分析
    interaction_success = True
    if not args.skip_interaction and env['python_packages']:
        interaction_success = run_interaction_analysis(args.output_dir)
    
    # 步骤3: MD系统准备
    md_success = True
    if not args.skip_md and env['gromacs'] and env['obabel']:
        docked_ligand = Path(args.output_dir) / 'docking' / 'docking_result.pdbqt'
        if docked_ligand.exists():
            md_success = prepare_md_system(args.protein, str(docked_ligand), args.output_dir)
        else:
            print("警告: 未找到对接结果，跳过MD准备")
    
    # 创建批处理和分析脚本
    create_batch_scripts(args.output_dir)
    create_analysis_scripts(args.output_dir)
    
    # 总结
    print("\n" + "="*70)
    print("流程完成总结")
    print("="*70)
    print(f"✓ 分子对接: {'完成' if docking_success else '跳过/失败'}")
    print(f"✓ 相互作用分析: {'完成' if interaction_success else '跳过/失败'}")
    print(f"✓ MD系统准备: {'完成' if md_success else '跳过/失败'}")
    print(f"✓ 批处理脚本: 已创建")
    print(f"✓ 分析脚本: 已创建")
    print(f"\n输出目录: {output_dir.absolute()}")
    print("="*70)
    print("\n下一步:")
    print("  1. 查看对接结果和相互作用图")
    print("  2. 在Linux/Mac环境下运行MD模拟")
    print("  3. 使用分析脚本处理MD轨迹")
    print("\n提示: 查看 README.md 获取详细使用说明")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
