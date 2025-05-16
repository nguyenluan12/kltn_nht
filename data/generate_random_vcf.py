import random
import os
import glob

def random_genotype():
    return f"{random.randint(0,1)}|{random.randint(0,1)}"

def generate_vcf(num_variants, num_samples):
    output_dir = "/Users/NguyenLuan/pbwt-comparison/random_data"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"random_{num_samples}x{num_variants}.vcf")

    with open(output_file, "w") as f:
        # Write VCF headers
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for i in range(num_samples):
            f.write(f"\tsample{i+1}")
        f.write("\n")

        # Write variant data
        for i in range(num_variants):
            f.write(f"chr21\t{10000 + i}\t.\tA\tG\t.\tPASS\t.\tGT")
            for _ in range(num_samples):
                f.write(f"\t{random_genotype()}")
            f.write("\n")

    print(f"[✓] File '{output_file}' generated with {num_variants} variants and {num_samples} samples.")

if __name__ == "__main__":
    output_dir = "/Users/NguyenLuan/pbwt-comparison/random_data"

    # 🧹 Xoá toàn bộ file .vcf cũ
    old_files = glob.glob(os.path.join(output_dir, "*.vcf"))
    for file in old_files:
        os.remove(file)
    print(f"🧹 Đã xoá {len(old_files)} file .vcf cũ trong thư mục {output_dir}")

    # Tạo các tập dữ liệu mới
    sizes = [
        (7, 5),
        
    ]
    for num_variants, num_samples in sizes:
        generate_vcf(num_variants, num_samples)
