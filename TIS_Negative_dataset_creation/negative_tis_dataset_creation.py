#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

###########################################################
# Parameters
###########################################################

genome_file = "../data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
annotation_file_name = "../intergenic_region_parsing/output/dist_GT_10000_LT_150000/final/6_Homo_sapiens.GRCh38.116_intergenic_regions_distGT_10000_LT_150000.bed"
working_dir = "output/intergenic_TIS/"

num_of_samples = 15000 # maximum number of samples
relative_start = 0     # in each intergenic region only a subregion is processed.This is the start processing position 
relative_end = 80000   # in each intergenic region only a subregion is processed.This is the end processing position 

upstream_context_length_list = [100, 300, 500] # the upstream length. run for a list of lengths. Each one produce diferect file
downstream_context_length = 500                # the downstream length

forbid_in_frame_stop_in_downstream_window = True # flag for check in frame stop codon
stop_codons = {"TAA", "TAG", "TGA"}              #stop codon list  

###########################################################
# Prepare output directory
###########################################################

os.makedirs(working_dir, exist_ok=True)

intergenic_subregion_bed_file_name = working_dir + "intergenic_subregions.bed"
intergenic_subregion_fa_file_name = working_dir  + "intergenic_subregions.fa"

###########################################################
# Step 1: Build BED file with extracted intergenic regions
# for each subregion from relativ start to relative end
###########################################################

count_regions = 0

with open(annotation_file_name, "r") as annotation_file, open(intergenic_subregion_bed_file_name, "w") as intergenic_region_bed_file:
    for line in annotation_file:
        if not line.strip():
            continue

        cols = line.rstrip().split("\t")
        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])

        if len(cols) > 3:
            name = cols[3]
        else:
            name = f"region_{count_regions + 1}"

        current_relative_end = min(relative_end, end - start)
        extract_region_start = start + relative_start
        extract_region_end = start + current_relative_end

        region_name = f"{name}:{chrom}:{extract_region_start}:{extract_region_end}:0:+"
        intergenic_region_bed_file.write(f"{chrom}\t{extract_region_start}\t{extract_region_end}\t{region_name}\n")

        count_regions += 1
        if count_regions >= num_of_samples:
            break

print("Intergenic regions written:", count_regions)
print("BED:", intergenic_subregion_bed_file_name)

###########################################################
# Step 2: Extract FASTA from BED
###########################################################

command = "bedtools getfasta -name -s -fi " + genome_file + " -bed " + intergenic_subregion_bed_file_name + " -fo " + intergenic_subregion_fa_file_name
print("[CMD]", command)
result = subprocess.run(command, shell=True)
if result.returncode != 0:
    raise RuntimeError("Command failed:\n" + command)

###########################################################
# Step 3: Scan for ATG and create candidate negative TIS
###########################################################

for upstream_context_length in upstream_context_length_list:
    count_tis = 0

    intergenic_tis_fa_file_name = working_dir +  "intergenic_negative_TIS_upstream-" + str(upstream_context_length) + "_downstream-" + str(downstream_context_length) + ".fa"
    intergenic_tis_bed_file_name = working_dir + "intergenic_negative_TIS_upstream-" + str(upstream_context_length) + "_downstream-" + str(downstream_context_length) + ".bed"

    with open(intergenic_subregion_fa_file_name, "r") as intergenic_region_fa_file, \
         open(intergenic_tis_fa_file_name, "w") as intergenic_tis_fa_file, \
         open(intergenic_tis_bed_file_name, "w") as intergenic_tis_bed_file:

        header_name = ""
        header_chrom = ""
        header_start = 0
        header_end = 0
        header_score = "0"
        header_strand = "+"

        stop_processing = False

        for line in intergenic_region_fa_file:
            line = line.rstrip().upper()

            if not line:
                continue

            if line.startswith(">"):
                header = line[1:]
                header_cols = header.split(":")
                header_name = header_cols[0]
                header_chrom = header_cols[1]
                header_start = int(header_cols[2])
                header_end = int(header_cols[3])
                header_score = header_cols[4]
                header_strand = header_cols[5]
                continue

            seq = line
            search_pos = 0 #start from the begining

            while True:
                atg_pos = seq.find("ATG", search_pos)
                if atg_pos == -1:
                    break

                window_start = atg_pos - upstream_context_length
                window_end = atg_pos + 3 + downstream_context_length

                if window_start < 0 or window_end > len(seq):
                    search_pos = atg_pos + 3
                    continue

                tis_window_seq = seq[window_start:window_end]

                if "N" in tis_window_seq:
                    search_pos = atg_pos + 3
                    continue

                if forbid_in_frame_stop_in_downstream_window:
                    has_in_frame_stop = False
                    scan_end = min(atg_pos + 3 + downstream_context_length, len(seq))

                    for i in range(atg_pos + 3, scan_end - 2, 3):
                        codon = seq[i:i+3]
                        if codon in stop_codons:
                            has_in_frame_stop = True
                            break

                    if has_in_frame_stop:
                        search_pos = atg_pos + 3
                        continue

                if header_strand == "+":
                    tis_start = header_start + window_start
                    tis_end = header_start + window_end
                elif header_strand == "-":
                    tis_start = header_end - window_end
                    tis_end = header_end - window_start
                else:
                    search_pos = atg_pos + 3
                    continue

                tis_name = header_name + ":" + header_chrom + ":" + str(tis_start) + ":" + str(tis_end) + ":ATGpos=" + str(atg_pos)

                intergenic_tis_fa_file.write(">" + tis_name + ":" + header_score + ":" + header_strand + "\n")
                intergenic_tis_fa_file.write(tis_window_seq + "\n")

                intergenic_tis_bed_file.write(header_chrom + "\t" + str(tis_start) + "\t" + str(tis_end) + "\t" + tis_name + "\t" + header_score + "\t" + header_strand + "\n")

                count_tis += 1
                if count_tis >= num_of_samples:
                    stop_processing = True
                    break

                search_pos = window_end

            if stop_processing:
                break

    print("Finished upstream =", upstream_context_length)
    print("  FASTA:", intergenic_tis_fa_file_name)
    print("  BED  :", intergenic_tis_bed_file_name)
    print("  Count:", count_tis)