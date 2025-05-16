#!/usr/bin/env bash
set -e

# Xác định thư mục gốc
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
VCF="$ROOT_DIR/random_data.vcf"
THREADS=(1 2 4 8)

PBWT="$ROOT_DIR/pbwt"
PPBWT="$ROOT_DIR/p_pbwt"

# 🔍 Tự động tìm file DLL dù là net7.0 hay net8.0
HPPBWT_DLL=$(find "$ROOT_DIR/hp-pbwt/source_IO_Excluded/bin/Release" -name "HP-PBWT.dll" | head -n 1)

if [ ! -f "$HPPBWT_DLL" ]; then
  echo "❌ Không tìm thấy HP-PBWT.dll trong thư mục Release!"
  exit 1
fi

RESULTS="$ROOT_DIR/results/results.csv"
mkdir -p "$ROOT_DIR/results_random_data"
echo "project,threads,time" > "$RESULTS"

for T in "${THREADS[@]}"; do
  echo "🧪 Threads = $T"

  echo -n "pbwt,$T," >> "$RESULTS"
  gtime -f "%e" "$PBWT" > /dev/null 2>> "$RESULTS"

  echo -n "p_pbwt,$T," >> "$RESULTS"
  gtime -f "%e" "$PPBWT" --vcf "$VCF" --threads "$T" > /dev/null 2>> "$RESULTS"

  # echo -n "hp-pbwt,$T," >> "$RESULTS"
  # gtime -f "%e" dotnet "$HPPBWT_DLL" "$VCF" "$T" 4 1 0 > /dev/null 2>> "$RESULTS"
done

echo "✅ Benchmark hoàn tất. Kết quả lưu tại: $RESULTS"
