import os
import subprocess
import csv
from pathlib import Path
import time

# Cấu hình
BASE_DIR = Path("/Users/NguyenLuan/pbwt-comparison")
DATA_DIR = BASE_DIR / "random_data"
RESULT_CSV = BASE_DIR / "pbwt_vs_hp_pbwt_results.csv"
LOG_FILE = Path("hp_pbwt_result.log")  # File log trong BASE_DIR

# Tham số chạy
LLM_LEN = 4
N_RUN = 1
COOLDOWN_MS = 0

# Di chuyển vào BASE_DIR để chệ chương trình ghi log đúng nơi
os.chdir(BASE_DIR)

# Xóa file CSV cũ (nếu có)
if RESULT_CSV.exists():
    RESULT_CSV.unlink()

# Header CSV
with open(RESULT_CSV, mode='w', newline='\n') as f:
    writer = csv.writer(f)
    writer.writerow(["algorithm", "threads", "data_size", "runtime_seconds"])

# Lặp qua tất cả file VCF trong folder
for vcf_path in sorted(DATA_DIR.glob("*.vcf")):
    data_size = vcf_path.stem.replace("random_", "")

    # Lặp qua số luồng từ 1, 2, 4, 6, 8
    for threads in [1, 2, 4, 6, 8]:
        print(f"\n[+] Running on {vcf_path.name} with {threads} thread(s)...")

        # Xóa log cũ (nếu có)
        if LOG_FILE.exists():
            LOG_FILE.unlink()

        cmd = [
            "dotnet", "run", "--project", str(BASE_DIR / "hp-pbwt/source_IO_Excluded/HP-PBWT.csproj"), "--",
            str(vcf_path), str(threads), str(LLM_LEN), str(N_RUN), str(COOLDOWN_MS)
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if LOG_FILE.exists():
                with open(LOG_FILE, 'r') as log:
                    last_line = log.readlines()[-1].strip()
                    parts = last_line.split('\t')
                    ms = int(parts[-1])
                    sec = ms / 1000
                    algorithm = "pbwt" if threads == 1 else "hp_pbwt"

                    with open(RESULT_CSV, mode='a', newline='\n') as f:
                        writer = csv.writer(f)
                        writer.writerow([algorithm, threads, data_size, sec])

                    print(f"[\u2713] Done: {sec:.3f} seconds")
            else:
                print("[!] Không tìm thấy file log")

        except subprocess.TimeoutExpired:
            print("[!] Timeout")
        except Exception as e:
            print(f"[!] Error: {e}")
