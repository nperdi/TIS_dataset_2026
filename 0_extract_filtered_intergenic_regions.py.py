"""
extract_filtered_intergenic_regions.py

Extract intergenic genomic regions from an annotation BED file by:
1. keeping standard chromosomes,
2. sorting intervals,
3. adding flanking slack,
4. merging overlapping/nearby annotations,
5. computing the complement,
6. filtering resulting regions by distance.
"""

import os
import sys
import subprocess

################## Set working dir same as script dir ##################

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

################## Helper for shell commands ###########################

def run_shell_command(cmd):
    print("[CMD]", cmd)
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise RuntimeError("Command failed:\n" + cmd)

################## Default vars ########################################

annotation_file = "data/annotations/refseq/Human_GRCh38_RefSeq_Curated.bed"
min_distance = 20000
max_distance = 50000
use_minmax = True

print("len(sys.argv) =", len(sys.argv))

if len(sys.argv) == 3:
    annotation_file = sys.argv[1]
    min_distance = int(sys.argv[2])
    use_minmax = False

elif len(sys.argv) == 4:
    annotation_file = sys.argv[1]
    min_distance = int(sys.argv[2])
    max_distance = int(sys.argv[3])
    use_minmax = True

elif len(sys.argv) != 1:
    print("Usage:")
    print("  python extract_filtered_intergenic_regions.py")
    print("  python extract_filtered_intergenic_regions.py <annotationFile> <minDistance>")
    print("  python extract_filtered_intergenic_regions.py <annotationFile> <minDistance> <maxDistance>")
    sys.exit(1)

if use_minmax and max_distance <= min_distance:
    raise ValueError("maxDistance must be greater than minDistance")

################## Paths / files #######################################

print("------", os.getcwd())
print(
    "\n#############     #############          #############\n "
    + os.path.abspath(__file__)
    + "\n #############     #############       ###############\n"
)

genome_file = "data/GENOME/hg38.fa"
chrom_sizes_file = "data/GENOME/hg38.chrom.size.clean"

allowed_chromosomes = set(
    ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
)

if not os.path.exists(annotation_file):
    raise FileNotFoundError("Annotation file not found: " + annotation_file)

if not os.path.exists(chrom_sizes_file):
    raise FileNotFoundError("Chrom sizes file not found: " + chrom_sizes_file)

base_name = os.path.splitext(os.path.basename(annotation_file))[0]

# General output root
output_root = "output"

# Run-specific output directory
if use_minmax:
    run_output_dir = os.path.join(
        output_root,
        f"dist_GT_{min_distance}_LT_{max_distance}"
    )
    final_output_name = (
        f"6_{base_name}_distGT_{min_distance}_LT_{max_distance}"
        f"_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
    )
else:
    run_output_dir = os.path.join(
        output_root,
        f"dist_GT_{min_distance}"
    )
    final_output_name = (
        f"6_{base_name}_distGT_{min_distance}"
        f"_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
    )

# Step-specific subdirectories
step1_dir = os.path.join(run_output_dir, "step1_chr_filter")
step2_dir = os.path.join(run_output_dir, "step2_slop")
step3_dir = os.path.join(run_output_dir, "step3_merge")
step4_dir = os.path.join(run_output_dir, "step4_complement")
step5_dir = os.path.join(run_output_dir, "step5_chr_filter")
final_dir = os.path.join(run_output_dir, "final")

for directory in [
    output_root,
    run_output_dir,
    step1_dir,
    step2_dir,
    step3_dir,
    step4_dir,
    step5_dir,
    final_dir,
]:
    os.makedirs(directory, exist_ok=True)

chr_annotation_bed_file = os.path.join(
    step1_dir, f"1_{base_name}_chrFilter.bed"
)

sorted_annotation_bed_file = os.path.join(
    step1_dir, f"1_{base_name}_sorted_chrFilter.bed"
)

slop_sorted_annotation_bed_file = os.path.join(
    step2_dir, f"2_{base_name}_slop_sorted_chrFilter.bed"
)

merge_slop_sorted_annotation_bed_file = os.path.join(
    step3_dir, f"3_{base_name}_merge_slop_sorted_chrFilter.bed"
)

complement_merge_slop_sorted_annotation_bed_file = os.path.join(
    step4_dir, f"4_{base_name}_complement_merge_slop_sorted_chrFilter.bed"
)

chr_filter_complement_merge_slop_sorted_annotation_bed_file = os.path.join(
    step5_dir, f"5_{base_name}_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
)

distance_filtered_output_file = os.path.join(
    final_dir, final_output_name
)

################## Initialize ##########################################

initialize = True

if initialize:
    with open(chr_annotation_bed_file, "w") as out_f:
        with open(annotation_file, "r") as in_f:
            for line in in_f:
                line = line.rstrip()
                if line == "":
                    continue

                cols = line.split("\t")
                if len(cols) < 3:
                    continue

                if cols[0] not in allowed_chromosomes:
                    continue

                out_f.write(line + "\n")

    full_command = (
        "bedtools sort -faidx "
        + chrom_sizes_file
        + " -i "
        + chr_annotation_bed_file
        + " > "
        + sorted_annotation_bed_file
    )
    run_shell_command(full_command)

    full_command = (
        "bedtools slop -b 2000 -g "
        + chrom_sizes_file
        + " -i "
        + sorted_annotation_bed_file
        + " > "
        + slop_sorted_annotation_bed_file
    )
    run_shell_command(full_command)

    full_command = (
        "bedtools merge -i "
        + slop_sorted_annotation_bed_file
        + " > "
        + merge_slop_sorted_annotation_bed_file
    )
    run_shell_command(full_command)

    full_command = (
        "bedtools complement -i "
        + merge_slop_sorted_annotation_bed_file
        + " -g "
        + chrom_sizes_file
        + " > "
        + complement_merge_slop_sorted_annotation_bed_file
    )
    run_shell_command(full_command)

    with open(chr_filter_complement_merge_slop_sorted_annotation_bed_file, "w") as out_f:
        with open(complement_merge_slop_sorted_annotation_bed_file, "r") as in_f:
            for line in in_f:
                line = line.rstrip()
                if line == "":
                    continue

                cols = line.split("\t")
                if len(cols) < 3:
                    continue

                if cols[0] not in allowed_chromosomes:
                    continue

                out_f.write(line + "\n")

################## Distance filter #####################################

count = 0

with open(distance_filtered_output_file, "w") as out_f:
    with open(chr_filter_complement_merge_slop_sorted_annotation_bed_file, "r") as in_f:
        for line in in_f:
            line = line.rstrip()
            if line == "":
                continue

            cols = line.split("\t")
            if len(cols) < 3:
                continue

            start = int(cols[1])
            end = int(cols[2])
            dist = end - start

            if use_minmax:
                if dist > min_distance and dist < max_distance:
                    count += 1
                    new_cols = cols + [
                        f"region_GT_{min_distance}_LT_{max_distance}_{count}",
                        "0",
                    ]
                    out_f.write("\t".join(new_cols) + "\n")
            else:
                if dist > min_distance:
                    count += 1
                    new_cols = cols + [
                        f"region_GT_{min_distance}_{count}",
                        "0",
                    ]
                    out_f.write("\t".join(new_cols) + "\n")

print(count)
print("Final output:", distance_filtered_output_file)