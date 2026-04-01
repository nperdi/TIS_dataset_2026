#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

########################################################################
# Set working directory = script directory
########################################################################

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

########################################################################
# Parameters
########################################################################

gtf_file = "../data/ensembl/Homo_sapiens.GRCh38.116.gtf"
chrom_sizes_file = "../data/ensembl/hg38.ensembl.chrom.size"

feature_type = "gene"
slop_bp = 2000

min_distance = 10000  # min length of extracted regions
max_distance = 150000 # max length of extracted regions
use_minmax = True     # if use max lengrh restriction

allowed_chromosomes = set([str(i) for i in range(1, 23)] + ["X", "Y"])

########################################################################
# Small helper for bedtools commands
########################################################################

def run_command(cmd_list, output_file=None):
    print("[CMD]", " ".join(cmd_list))

    if output_file is not None:
        with open(output_file, "w") as fout:
            result = subprocess.run(cmd_list, stdout=fout, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd_list, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(cmd_list)
            + "\n\nstderr:\n"
            + result.stderr
        )

########################################################################
# Checks
########################################################################

if not os.path.exists(gtf_file):
    raise FileNotFoundError("GTF file not found: " + gtf_file)

if not os.path.exists(chrom_sizes_file):
    raise FileNotFoundError("Chrom sizes file not found: " + chrom_sizes_file)

if use_minmax and max_distance <= min_distance:
    raise ValueError("max_distance must be greater than min_distance")

########################################################################
# Output paths
########################################################################

base_name = os.path.splitext(os.path.basename(gtf_file))[0]
output_root = "output"

if use_minmax:
    run_output_dir = os.path.join(output_root, f"dist_GT_{min_distance}_LT_{max_distance}")
    final_output_name = (
        f"6_{base_name}_intergenic_regions_distGT_{min_distance}_LT_{max_distance}.bed"
    )
else:
    run_output_dir = os.path.join(output_root, f"dist_GT_{min_distance}")
    final_output_name = (
        f"6_{base_name}_intergenic_regions__distGT_{min_distance}.bed"
    )

step0_dir = os.path.join(run_output_dir, "step0_gtf_to_bed")
step1_dir = os.path.join(run_output_dir, "step1_sort")
step2_dir = os.path.join(run_output_dir, "step2_slop")
step3_dir = os.path.join(run_output_dir, "step3_merge")
step4_dir = os.path.join(run_output_dir, "step4_complement")
step5_dir = os.path.join(run_output_dir, "step5_chr_filter")
final_dir = os.path.join(run_output_dir, "final")

for d in [output_root, run_output_dir, step0_dir, step1_dir, step2_dir, step3_dir, step4_dir, step5_dir, final_dir]:
    os.makedirs(d, exist_ok=True)

gtf_annotation_bed_file = os.path.join(step0_dir, f"0_{base_name}_{feature_type}.bed")
sorted_annotation_bed_file = os.path.join(step1_dir, f"1_{base_name}_{feature_type}_sorted.bed")
slop_sorted_annotation_bed_file = os.path.join(step2_dir, f"2_{base_name}_{feature_type}_slop_sorted.bed")
merge_slop_sorted_annotation_bed_file = os.path.join(step3_dir, f"3_{base_name}_{feature_type}_merge_slop_sorted.bed")
complement_merge_slop_sorted_annotation_bed_file = os.path.join(step4_dir, f"4_{base_name}_{feature_type}_complement_merge_slop_sorted.bed")
chr_filter_complement_merge_slop_sorted_annotation_bed_file = os.path.join(step5_dir, f"5_{base_name}_{feature_type}_chrFilter_complement_merge_slop_sorted.bed")
distance_filtered_output_file = os.path.join(final_dir, final_output_name)

########################################################################
# Step 0: GTF -> BED
########################################################################

print("Converting GTF to BED...")

count = 0

with open(gtf_file, "r") as in_f, open(gtf_annotation_bed_file, "w") as out_f:
    for line in in_f:
        if line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        chrom = cols[0]
        feature = cols[2]

        if feature != feature_type:  ##feature_type==gene
            continue

        if chrom not in allowed_chromosomes: ## 1,2,3,4....
            continue

        start_1based = int(cols[3])
        end_1based = int(cols[4])
        score = cols[5]
        strand = cols[6]
        attributes = cols[8]

        start_0based = start_1based - 1
        end_bed = end_1based

        if start_0based < 0:
            start_0based = 0

        gene_name = ""
        gene_id = ""

        attr_fields = attributes.strip().split(";")
        for field in attr_fields:
            field = field.strip()
            if field.startswith("gene_name "):
                gene_name = field.replace("gene_name ", "").strip().strip('"')
            elif field.startswith("gene_id "):
                gene_id = field.replace("gene_id ", "").strip().strip('"')

        name = gene_name if gene_name != "" else gene_id
        if name == "":
            name = f"gene_{count+1}"

        if score == ".":
            score = "0"

        out_f.write("\t".join([chrom, str(start_0based), str(end_bed), name, score, strand ]) + "\n")
        count += 1

print("Extracted intervals:", count)

########################################################################
# Step 1: sort
########################################################################

run_command(["bedtools", "sort", "-faidx", chrom_sizes_file, "-i", gtf_annotation_bed_file], output_file=sorted_annotation_bed_file )

########################################################################
# Step 2: slop
########################################################################

run_command(["bedtools", "slop", "-b", str(slop_bp), "-g", chrom_sizes_file, "-i", sorted_annotation_bed_file],output_file=slop_sorted_annotation_bed_file)

########################################################################
# Step 3: merge
########################################################################

run_command(["bedtools", "merge", "-i", slop_sorted_annotation_bed_file],output_file=merge_slop_sorted_annotation_bed_file)

########################################################################
# Step 4: complement
########################################################################

run_command(["bedtools", "complement", "-i", merge_slop_sorted_annotation_bed_file, "-g", chrom_sizes_file], output_file=complement_merge_slop_sorted_annotation_bed_file)

########################################################################
# Step 5: chromosome filter
########################################################################

with open(complement_merge_slop_sorted_annotation_bed_file, "r") as in_f, \
     open(chr_filter_complement_merge_slop_sorted_annotation_bed_file, "w") as out_f:

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

########################################################################
# Step 6: distance filter
########################################################################

count = 0

with open(chr_filter_complement_merge_slop_sorted_annotation_bed_file, "r") as in_f, \
     open(distance_filtered_output_file, "w") as out_f:

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
            if min_distance < dist < max_distance:
                count += 1
                new_cols = cols + [f"region_GT_{min_distance}_LT_{max_distance}_{count}", "0"]
                out_f.write("\t".join(new_cols) + "\n")
        else:
            if dist > min_distance:
                count += 1
                new_cols = cols + [f"region_GT_{min_distance}_{count}", "0"]
                out_f.write("\t".join(new_cols) + "\n")

print("Final region count:", count)
print("Final output:", distance_filtered_output_file)