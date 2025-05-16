import os
import subprocess
import re
import csv
import argparse

# Thư mục chứa dữ liệu VCF và đường dẫn kết quả
DATA_DIR = "/Users/NguyenLuan/pbwt-comparison/random_data"
RESULT_CSV = "/Users/NguyenLuan/pbwt-comparison/pbwt_vs_p_pbwt_results.csv"

# Lấy thời gian từ dòng stdout (us → giây)
def extract_total_time(stdout_text):
    match = re.search(r"total time: (\d+) us", stdout_text, re.IGNORECASE)
    if match:
        return int(match.group(1)) / 1e6  # convert to seconds
    return None

# Chạy chương trình và lấy stdout
def run_binary_get_output(binary_path, input_file):
    result = subprocess.run([binary_path, input_file], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return result.stdout

# Ghi kết quả vào CSV
def append_result(csv_path, algorithm, threads, data_size, runtime):
    with open(csv_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([algorithm, threads, data_size, f"{runtime:.6f}"])

# Hàm chính
def main():
    parser = argparse.ArgumentParser(description="Run PBWT and P-PBWT experiments and collect results")
    parser.add_argument("--pbwt", required=True, help="Path to pbwt executable")
    parser.add_argument("--p_pbwt", required=True, help="Path to p_pbwt executable")
    parser.add_argument("--output_csv", default=RESULT_CSV, help="Path to save CSV results")
    args = parser.parse_args()

    if not os.path.exists(args.output_csv):
        with open(args.output_csv, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["algorithm", "threads", "data_size", "runtime_seconds"])

    for filename in sorted(os.listdir(DATA_DIR)):
        if not filename.endswith(".vcf"):
            continue
        filepath = os.path.join(DATA_DIR, filename)
        data_size = filename.replace("random_", "").replace(".vcf", "")

        # ----- Chạy PBWT -----
        output = run_binary_get_output(args.pbwt, filepath)
        runtime = extract_total_time(output)
        if runtime:
            append_result(args.output_csv, "pbwt", 1, data_size, runtime)

        # ----- Chạy P-PBWT -----
        output = run_binary_get_output(args.p_pbwt, filepath)
        for match in re.finditer(r"Num thread:(\d+).*?total time: (\d+) us", output, re.DOTALL):
            thread = int(match.group(1))
            runtime = int(match.group(2)) / 1e6
            append_result(args.output_csv, "p_pbwt", thread, data_size, runtime)

    print(f"[✓] Kết quả đã lưu vào {args.output_csv}")

if __name__ == "__main__":
    main()
