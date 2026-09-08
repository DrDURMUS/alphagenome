#!/bin/bash
# ==============================================================================
#  AlphaGenome Studio - Masaüstü Başlatıcı Oluşturucu (macOS)
# ==============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DESKTOP_DIR="$HOME/Desktop"
APP_NAME="AlphaGenome Studio.app"
TARGET_APP="$DESKTOP_DIR/$APP_NAME"

echo "AlphaGenome Studio macOS masaüstü uygulaması oluşturuluyor..."

if command -v osacompile &> /dev/null; then
    rm -rf "$TARGET_APP"
    osacompile -o "$TARGET_APP" -e "tell application \"Terminal\"
        activate
        do script \"bash \\\"$SCRIPT_DIR/start_studio.sh\\\"\"
    end tell"
    
    if [ -f "$SCRIPT_DIR/assets/app.icns" ]; then
        cp "$SCRIPT_DIR/assets/app.icns" "$TARGET_APP/Contents/Resources/applet.icns"
        touch "$TARGET_APP"
    fi
    echo "✅ Masaüstünüze uygulama simgesi başarıyla eklendi: $TARGET_APP"
else
    cp "$SCRIPT_DIR/AlphaGenome Studio.command" "$DESKTOP_DIR/"
    chmod +x "$DESKTOP_DIR/AlphaGenome Studio.command"
    echo "✅ Masaüstünüze komut kısayolu eklendi: $DESKTOP_DIR/AlphaGenome Studio.command"
fi
