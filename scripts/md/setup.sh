#!/bin/bash
# MD系统准备脚本
# 包含: 拓扑生成、溶剂化、离子化、能量最小化

set -e  # 遇到错误立即退出

# 默认参数
PROTEIN=""
LIGAND=""
OUTPUT="results/md"
FORCEFIELD="amber99sb-ildn"
WATER="tip3p"
BOX_SIZE=1.0  # nm
ION_CONC=0.15  # M

# 显示帮助
show_help() {
    cat << EOF
MD系统准备脚本 (Gromacs)

使用方法:
    bash setup.sh -p protein.pdb -l ligand.pdbqt [选项]

必需参数:
    -p, --protein       受体蛋白PDB文件
    -l, --ligand        配体PDBQT文件 (对接后的最佳构象)

可选参数:
    -o, --output        输出目录 (默认: results/md)
    -f, --forcefield    力场 (默认: amber99sb-ildn)
    -w, --water         水模型 (默认: tip3p)
    -b, --box_size      盒子边界 (默认: 1.0 nm)
    -i, --ion_conc      离子浓度 (默认: 0.15 M)
    -h, --help          显示帮助

示例:
    # 基本使用
    bash setup.sh -p protein.pdb -l top1_docked.pdbqt
    
    # 自定义参数
    bash setup.sh -p protein.pdb -l ligand.pdbqt \\
                  -o md_results -f charmm36-jul2022 -b 1.2

计算时间估计:
    - 拓扑生成: ~1-2分钟
    - 溶剂化: ~30秒
    - 离子化: ~30秒
    - 能量最小化: ~10-30分钟
    - 总计: ~15-35分钟
EOF
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--protein)
            PROTEIN="$2"
            shift 2
            ;;
        -l|--ligand)
            LIGAND="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -f|--forcefield)
            FORCEFIELD="$2"
            shift 2
            ;;
        -w|--water)
            WATER="$2"
            shift 2
            ;;
        -b|--box_size)
            BOX_SIZE="$2"
            shift 2
            ;;
        -i|--ion_conc)
            ION_CONC="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "未知参数: $1"
            show_help
            exit 1
            ;;
    esac
done

# 检查必需参数
if [[ -z "$PROTEIN" ]] || [[ -z "$LIGAND" ]]; then
    echo "错误: 必须提供蛋白和配体文件"
    show_help
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT"

echo "========================================"
echo "MD系统准备 (Gromacs)"
echo "========================================"
echo "蛋白: $PROTEIN"
echo "配体: $LIGAND"
echo "力场: $FORCEFIELD"
echo "水模型: $WATER"
echo "盒子边界: $BOX_SIZE nm"
echo "离子浓度: $ION_CONC M"
echo "输出目录: $OUTPUT"
echo "========================================"

# 步骤1: 生成蛋白拓扑
echo -e "\n[1/6] 生成蛋白拓扑..."
gmx pdb2gmx -f "$PROTEIN" -o "$OUTPUT/protein.gro" -p "$OUTPUT/topol.top" \
    -i "$OUTPUT/posre.itp" -ff "$FORCEFIELD" -water "$WATER" -ignh << EOF
1
EOF

# 步骤2: 转换配体格式
echo -e "\n[2/6] 准备配体..."
# 将PDBQT转换为PDB
obabel "$LIGAND" -O "$OUTPUT/ligand.pdb"

# 使用ACPYPE生成配体拓扑
echo "生成配体拓扑 (GAFF)..."
mkdir -p "$OUTPUT/ligand"
cd "$OUTPUT/ligand"

# 检查acpype是否可用
if command -v acpype &> /dev/null; then
    acpype -i ../ligand.pdb -c bcc -a gaff2
    cd ../..
    
    # 整合拓扑文件
    cp "$OUTPUT/ligand/ligand.acpype/ligand_GMX.itp" "$OUTPUT/ligand.itp"
    cp "$OUTPUT/ligand/ligand.acpype/ligand_GMX.gro" "$OUTPUT/ligand.gro"
else
    echo "警告: acpype未安装,使用简化方法"
    echo "建议安装: pip install acpype"
    cd ../..
    # 这里可以添加替代方案
fi

# 步骤3: 合并蛋白和配体
echo -e "\n[3/6] 合并蛋白-配体复合物..."
# 合并gro文件
head -n -1 "$OUTPUT/protein.gro" > "$OUTPUT/complex.gro"
cat "$OUTPUT/ligand.gro" >> "$OUTPUT/complex.gro"

# 更新拓扑文件
cat >> "$OUTPUT/topol.top" << EOF

; Ligand topology
#include "ligand.itp"

[ molecules ]
; Compound        nmols
Protein_chain_A     1
LIG                 1
EOF

# 步骤4: 溶剂化
echo -e "\n[4/6] 溶剂化..."
gmx editconf -f "$OUTPUT/complex.gro" -o "$OUTPUT/complex_box.gro" \
    -c -d "$BOX_SIZE" -bt cubic

gmx solvate -cp "$OUTPUT/complex_box.gro" -cs spc216.gro \
    -o "$OUTPUT/complex_solv.gro" -p "$OUTPUT/topol.top"

# 步骤5: 添加离子
echo -e "\n[5/6] 添加离子..."
# 创建ions.mdp
 cat > "$OUTPUT/ions.mdp" << EOF
; 离子化参数
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

; 邻域列表
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12
rcoulomb        = 1.0
rvdw            = 1.0
EOF

gmx grompp -f "$OUTPUT/ions.mdp" -c "$OUTPUT/complex_solv.gro" \
    -p "$OUTPUT/topol.top" -o "$OUTPUT/ions.tpr" -maxwarn 2

gmx genion -s "$OUTPUT/ions.tpr" -o "$OUTPUT/complex_ions.gro" \
    -p "$OUTPUT/topol.top" -pname NA -nname CL -neutral \
    -conc "$ION_CONC" << EOF
13
EOF

# 步骤6: 能量最小化
echo -e "\n[6/6] 能量最小化..."
cat > "$OUTPUT/em.mdp" << EOF
; 能量最小化参数
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

; 邻域列表
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.12
rcoulomb        = 1.0
rvdw            = 1.0

; 约束
constraints     = h-bonds
constraint_algorithm = LINCS
EOF

gmx grompp -f "$OUTPUT/em.mdp" -c "$OUTPUT/complex_ions.gro" \
    -p "$OUTPUT/topol.top" -o "$OUTPUT/em.tpr" -maxwarn 2

echo "运行能量最小化..."
gmx mdrun -v -deffnm "$OUTPUT/em" -c "$OUTPUT/em.gro"

# 检查能量最小化是否收敛
if [ -f "$OUTPUT/em.gro" ]; then
    echo -e "\n能量最小化完成!"
    
    # 分析能量
    echo "能量分析:"
    gmx energy -f "$OUTPUT/em.edr" -o "$OUTPUT/potential.xvg" << EOF
10
0
EOF
    
    echo "势能图: $OUTPUT/potential.xvg"
else
    echo "错误: 能量最小化失败"
    exit 1
fi

echo -e "\n========================================"
echo "系统准备完成!"
echo "========================================"
echo "输出文件:"
echo "  拓扑: $OUTPUT/topol.top"
echo "  结构: $OUTPUT/em.gro"
echo "  配体拓扑: $OUTPUT/ligand.itp"
echo "========================================"
echo -e "\n下一步: 运行NVT平衡"
echo "  bash run_md.sh --input $OUTPUT/em.gro --topology $OUTPUT/topol.top"
