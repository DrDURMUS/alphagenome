#!/bin/bash
# ==============================================================================
#  AlphaGenome Studio - Başlatıcı Betiği
# ==============================================================================

# Script'in bulunduğu dizin
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$SCRIPT_DIR" || exit 1

clear
echo "================================================================"
echo "  🧬 ALPHAGENOME STUDIO"
echo "  Atlas & Model Entegre Varyant Analiz Platformu"
echo "================================================================"
echo "  🚀 Başlatılıyor..."
echo "  🌐 Web Arayüzü: http://127.0.0.1:8000"
echo "  🛑 Durdurmak için: Klavyeden Ctrl + C tuşlarına basın"
echo "================================================================"
echo ""

# Python yolunu dinamik tespit et
PYTHON_CMD=""

# 1. alphagenome conda ortamı yolları kontrolü
CANDIDATE_PATHS=(
    "$HOME/miniconda3/envs/alphagenome/bin/python"
    "$HOME/anaconda3/envs/alphagenome/bin/python"
    "$HOME/miniforge3/envs/alphagenome/bin/python"
    "/opt/homebrew/Caskroom/miniconda/base/envs/alphagenome/bin/python"
)

for p in "${CANDIDATE_PATHS[@]}"; do
    if [ -x "$p" ]; then
        PYTHON_CMD="$p"
        break
    fi
done

# 2. Eğer bulunamadıysa aktif conda ortamı üzerinden kontrol et
if [ -z "$PYTHON_CMD" ] && command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook 2>/dev/null)"
    conda activate alphagenome 2>/dev/null
    if [ "$CONDA_DEFAULT_ENV" = "alphagenome" ]; then
        PYTHON_CMD="python"
    fi
fi

# 3. Bulunamazsa sistemdeki python / python3
if [ -z "$PYTHON_CMD" ]; then
    if command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        PYTHON_CMD="python3"
    fi
fi

exec "$PYTHON_CMD" web_app.py --port 8000
