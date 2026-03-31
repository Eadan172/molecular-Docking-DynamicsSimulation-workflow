@echo off
chcp 65001 >nul
echo ========================================
echo 分子对接与动力学模拟 Demo - 快速开始
echo ========================================
echo.

REM 检查Python
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [错误] 未找到Python，请先安装Python 3.8+
    pause
    exit /b 1
)
echo [√] Python已安装

REM 检查依赖
python -c "import meeko" >nul 2>&1
if %errorlevel% neq 0 (
    echo [警告] 缺少依赖包，正在安装...
    pip install -r requirements.txt
)

echo.
echo 选择操作:
echo 1. 查看项目说明 (README.md)
echo 2. 运行相互作用可视化演示
echo 3. 运行完整流程 (需要AutoDock Vina)
echo 4. 退出
echo.
set /p choice=请输入选项 (1-4): 

if "%choice%"=="1" (
    notepad README.md
    goto end
)

if "%choice%"=="2" (
    echo.
    echo 正在运行相互作用可视化...
    python scripts\visualization\visualize_interactions.py
    goto end
)

if "%choice%"=="3" (
    echo.
    echo 正在运行完整流程...
    python run_full_pipeline.py
    goto end
)

if "%choice%"=="4" (
    goto end
)

echo [错误] 无效选项
:end
echo.
echo 按任意键退出...
pause >nul
