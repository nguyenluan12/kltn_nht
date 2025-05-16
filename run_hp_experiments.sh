#!/bin/bash

# Thư mục cố định
BASE_DIR="/Users/NguyenLuan/pbwt-comparison"
PY_SCRIPT="$BASE_DIR/run_hp_pbwt.py"

# Kiểm tra dotnet
if ! command -v dotnet &> /dev/null; then
    echo "❌ dotnet không được cài đặt trong PATH"
    exit 1
fi

# Thông báo bắt đầu
echo "🚀 Bắt đầu chạy thực nghiệm HP-PBWT..."

# Chạy file Python
python3 "$PY_SCRIPT"

# Kiểm tra kết quả
if [[ $? -eq 0 ]]; then
    echo "✅ Đã hoàn tất. Kết quả được lưu tại: $BASE_DIR/pbwt_vs_hp_pbwt_results.csv"
else
    echo "❌ Gặp lỗi trong quá trình thực nghiệm."
fi
