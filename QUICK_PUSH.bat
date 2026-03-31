@echo off
chcp 65001 >nul
echo ========================================
echo 快速推送到 GitHub
echo ========================================
echo.
echo 当前目录: %CD%
echo.
echo 正在尝试推送到 GitHub...
echo.
git push -u origin main
echo.
if %errorlevel% neq 0 (
    echo.
    echo ========================================
    echo 推送失败！可能原因：
    echo 1. 网络连接问题（需要VPN）
    echo 2. GitHub认证失败
    echo 3. 仓库不存在
    echo.
    echo 请参考 UPLOAD_GITHUB.md 获取帮助
    echo ========================================
    echo.
    echo 替代方案：
    echo 1. 使用 GitHub Desktop（推荐）
    echo 2. 使用 VS Code 的源代码管理
    echo 3. 先解决网络问题再重试
    echo.
) else (
    echo.
    echo ========================================
    echo ✓ 推送成功！
    echo ========================================
    echo.
    echo 项目已上传到：
    echo https://github.com/Eadan172/molecular-Docking-DynamicsSimulation-workflow
    echo.
)
pause
