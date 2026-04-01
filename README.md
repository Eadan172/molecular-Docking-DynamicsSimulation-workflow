# 分子对接与动力学模拟 Demo

一个完整的分子对接和分子动力学模拟工具包，适用于面试演示和实际科研工作。

## 📋 项目功能

### 1. 分子对接 (AutoDock Vina)
- 从 SDF 配体和 PDB 蛋白直接运行对接
- 支持自定义网格中心、大小和搜索详尽度
- 单核心和并行对接模式

### 2. 相互作用分析 (Meeko + RDKit)
- 可视化蛋白质-配体相互作用
- 显示相互作用力类型（氢键、疏水、π-π等）
- 显示能量大小和参与结合残基类别
- 生成综合分析图表

### 3. 分子动力学模拟 (Gromacs)
- 自动生成全部 MD 输入文件
- 支持力场选择（AMBER99SB-ILDN、CHARMM36等）
- 完整的溶剂化和离子化流程
- 提供集群提交脚本

### 4. 动力学结果分析
- **RMSD 分析**: 蛋白和配体稳定性评估
- **PCA 分析**: 主要运动模式可视化
- **自由能形貌图 (FES)**: 能量景观分析

## 📁 项目结构

```
molecular_docking_demo/
├── data/                          # 输入数据
│   ├── protein.pdb               # 受体蛋白
│   └── ligands.sdf               # 配体库
├── scripts/                       # 脚本目录
│   ├── docking/                  # 对接脚本
│   │   ├── single_core.py        # 单核心对接
│   │   └── parallel.py           # 并行对接
│   ├── analysis/                 # 分析脚本
│   │   ├── rmsd.py               # RMSD分析
│   │   ├── pca.py                # PCA分析
│   │   └── fes.py                # 自由能形貌图
│   ├── md/                       # MD脚本
│   │   └── setup.sh              # MD系统构建
│   └── visualization/            # 可视化脚本
│       └── visualize_interactions.py  # 相互作用可视化
├── params/                       # 参数文件
│   ├── nvt.mdp                   # NVT平衡参数
│   ├── npt.mdp                   # NPT平衡参数
│   └── md.mdp                    # MD生产模拟参数
├── results/                      # 输出结果
│   ├── docking/                  # 对接结果
│   ├── visualization/            # 可视化结果
│   └── md/                       # MD结果
├── examples/                     # 示例数据
├── docs/                         # 文档
├── run_full_pipeline.py          # 完整流程脚本
└── README.md
```

## 🛠️ 安装依赖

### 必需软件
- **AutoDock Vina**: 分子对接
- **GROMACS**: 分子动力学模拟
- **OpenBabel**: 格式转换

### Python 包
```bash
pip install meeko rdkit numpy matplotlib seaborn scikit-learn
```

## 🚀 快速开始

### 方法1: 使用完整流程脚本
```bash
# 运行完整流程（使用示例数据）
python run_full_pipeline.py

# 自定义参数
python run_full_pipeline.py \
    --protein data/protein.pdb \
    --ligand data/ligands.sdf \
    --center 24 22 17 \
    --box_size 20 20 20 \
    --exhaustiveness 16
```

### 方法2: 单独运行各模块

#### 1. 分子对接
```bash
python scripts/docking/single_core.py \
    -r data/protein.pdb \
    -l data/ligands.sdf \
    -c 24 22 17 \
    -s 20 20 20 \
    -o results/docking
```

#### 2. 相互作用分析
```bash
python scripts/visualization/visualize_interactions.py
```

#### 3. MD 系统准备 (Linux/Mac)
```bash
cd scripts/md
bash setup.sh -p ../../data/protein.pdb -l ../../results/docking/docking_result.pdbqt -o ../../results/md
```

#### 4. 轨迹分析
```bash
# RMSD 分析
python scripts/analysis/rmsd.py

# PCA 分析
python scripts/analysis/pca.py

# 自由能形貌图
python scripts/analysis/fes.py
```

## 📊 输出结果

### 分子对接结果
- `results/docking/docking_result.pdbqt`: 最佳对接构象
- `results/docking/docking_log.txt`: 对接日志

### 相互作用可视化
- `results/visualization/interaction_analysis_overview.png`: 综合分析图
- `results/visualization/interaction_types.png`: 相互作用类型分布
- `results/visualization/docking_energies.png`: 对接能量分布
- `results/visualization/interaction_network.png`: 相互作用网络

### MD 结果
- `results/md/topol.top`: 拓扑文件
- `results/md/em.gro`: 能量最小化结构
- `results/md/nvt.gro`: NVT平衡结构
- `results/md/npt.gro`: NPT平衡结构

## 🎯 集群使用

### 提交 MD 作业
```bash
# 使用 SLURM
sbatch results/submit_md.sh

# 使用 PBS
qsub results/submit_md.sh
```

## 📝 注意事项

1. **Windows 用户**: MD 模拟需要在 Linux/Mac 或 WSL 中运行
2. **网格参数**: 根据你的蛋白和配体位置调整网格中心和大小
3. **计算资源**: 完整的 MD 模拟可能需要数小时到数天
4. **力场选择**: 根据研究体系选择合适的力场

## 📄 许可证

MIT License
