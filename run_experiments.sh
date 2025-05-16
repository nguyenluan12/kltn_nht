#!/bin/bash

# Đường dẫn cố định
BASE_DIR="/Users/NguyenLuan/pbwt-comparison"
DATA_DIR="$BASE_DIR/random_data"
RESULT_CSV="$BASE_DIR/pbwt_vs_p_pbwt_results.csv"
PBWT_EXE="$BASE_DIR/pbwt-src/pbwt"
P_PBWT_EXE="$BASE_DIR/pbwt-src/p_pbwt"
PY_SCRIPT="$BASE_DIR/experiment_pipeline.py"

# Kiểm tra tồn tại binary
if [[ ! -f "$PBWT_EXE" ]]; then
    echo "❌ Không tìm thấy file: $PBWT_EXE"
    exit 1
fi

if [[ ! -f "$P_PBWT_EXE" ]]; then
    echo "❌ Không tìm thấy file: $P_PBWT_EXE"
    exit 1
fi

# Xoá file CSV cũ nếu tồn tại
if [[ -f "$RESULT_CSV" ]]; then
    echo "🗑️ Đang xoá file kết quả cũ: $RESULT_CSV"
    rm "$RESULT_CSV"
fi

# In thông báo bắt đầu
echo "🚀 Chạy thực nghiệm PBWT vs P-PBWT..."

# Chạy script thực nghiệm Python
python3 "$PY_SCRIPT" \
    --pbwt "$PBWT_EXE" \
    --p_pbwt "$P_PBWT_EXE" \
    --output_csv "$RESULT_CSV"

# Thông báo kết thúc
if [[ $? -eq 0 ]]; then
    echo "✅ Đã hoàn tất thực nghiệm. Kết quả lưu tại: $RESULT_CSV"
else
    echo "❌ Có lỗi xảy ra trong quá trình thực nghiệm."
fi