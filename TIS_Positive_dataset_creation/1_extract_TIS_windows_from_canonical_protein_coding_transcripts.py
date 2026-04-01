#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import os
import re
import subprocess
from collections import defaultdict

###############################################################################
# ΠΑΡΑΜΕΤΡΟΙ
###############################################################################

GTF_FILE = "../data/ensembl/Homo_sapiens.GRCh38.116.gtf"
GENOME_FASTA = "../data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

UPSTREAM = 500
DOWNSTREAM = 500

OUT_PREFIX = "output/canonical_trancripts_tis_window"

###############################################################################
# ΒΟΗΘΗΤΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ
###############################################################################

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_gtf_attributes(attr_text):
    attrs = {}
    tags = []

    for m in re.finditer(r'(\S+)\s+"([^"]*)";', attr_text):
        key, value = m.group(1), m.group(2)
        if key == "tag":
            tags.append(value)
        else:
            attrs[key] = value

    attrs["tag"] = tags
    return attrs


def reverse_complement(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def wrap_fasta(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def load_fai_lengths(fai_file):
    chrom_lengths = {}
    with open(fai_file, "r", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            length = int(fields[1])
            chrom_lengths[chrom] = length
    return chrom_lengths


def fetch_sequence_samtools(fasta_file, chrom, start_1based, end_1based):
    region = f"{chrom}:{start_1based}-{end_1based}"

    result = subprocess.run(
        ["samtools", "faidx", fasta_file, region],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"samtools faidx failed for region {region}\n"
            f"stderr:\n{result.stderr}"
        )

    seq_lines = result.stdout.splitlines()
    seq = "".join(line.strip() for line in seq_lines if not line.startswith(">"))
    return seq.upper()


###############################################################################
# ΚΥΡΙΑ ΛΟΓΙΚΗ
###############################################################################

def main():
    genome_fai = GENOME_FASTA + ".fai"

    output_dir = os.path.dirname(OUT_PREFIX)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    output_tsv = OUT_PREFIX + ".tsv"
    output_bed = OUT_PREFIX + ".bed"
    output_fasta = OUT_PREFIX + ".fa"

    output_no_overlap_tsv = OUT_PREFIX + ".no_overlap.tsv"
    output_no_overlap_bed = OUT_PREFIX + ".no_overlap.bed"
    output_no_overlap_fasta = OUT_PREFIX + ".no_overlap.fa"

    chrom_lengths = load_fai_lengths(genome_fai)
    fasta_keys = set(chrom_lengths.keys())

    print(f"[INFO] FASTA sequences in index: {len(fasta_keys)}")
    print(f"[INFO] Example FASTA names: {list(chrom_lengths.keys())[:10]}")

    canonical_transcripts = {}
    skipped_no_chrom_match = 0

    with open_maybe_gzip(GTF_FILE) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr_text = fields

            if feature != "transcript":
                continue

            attrs = parse_gtf_attributes(attr_text)

            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            transcript_id = attrs.get("transcript_id", "")
            transcript_name = attrs.get("transcript_name", "")
            gene_biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            transcript_biotype = attrs.get("transcript_biotype", attrs.get("transcript_type", ""))
            tags = attrs.get("tag", [])

            if gene_biotype != "protein_coding":
                continue
            if transcript_biotype != "protein_coding":
                continue
            if "Ensembl_canonical" not in tags:
                continue

            if chrom is None:
                skipped_no_chrom_match += 1
                continue

            canonical_transcripts[transcript_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "transcript_id": transcript_id,
                "transcript_name": transcript_name,
                "chrom_gtf": chrom,
                "chrom_fasta": chrom,
                "strand": strand,
            }

    print(f"[INFO] Canonical protein-coding transcripts kept: {len(canonical_transcripts)}")
    print(f"[INFO] Skipped transcripts due to chromosome mismatch: {skipped_no_chrom_match}")

    start_codons = defaultdict(list)

    with open_maybe_gzip(GTF_FILE) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attr_text = fields

            if feature != "start_codon":
                continue

            attrs = parse_gtf_attributes(attr_text)
            transcript_id = attrs.get("transcript_id", "")

            if transcript_id not in canonical_transcripts:
                continue

            start_codons[transcript_id].append((int(start), int(end)))

    records = []
    skipped_no_start = 0
    skipped_bad_seq = 0

    for transcript_id, meta in canonical_transcripts.items():
        coords = start_codons.get(transcript_id, [])
        if not coords:
            skipped_no_start += 1
            continue

        chrom_gtf = meta["chrom_gtf"]
        chrom_fasta = meta["chrom_fasta"]
        strand = meta["strand"]

        starts = [x[0] for x in coords]
        ends = [x[1] for x in coords]

        if strand == "+":
            tis = min(starts)
            window_start = tis - UPSTREAM
            window_end = tis + DOWNSTREAM - 1
        elif strand == "-":
            tis = max(ends)
            window_start = tis - DOWNSTREAM + 1
            window_end = tis + UPSTREAM
        else:
            skipped_bad_seq += 1
            continue

        chrom_len = chrom_lengths[chrom_fasta]

        window_start = max(1, window_start)
        window_end = min(chrom_len, window_end)

        if window_start > window_end:
            skipped_bad_seq += 1
            continue

        try:
            seq = fetch_sequence_samtools(GENOME_FASTA, chrom_fasta, window_start, window_end)
        except Exception as e:
            print(f"[WARN] Could not fetch sequence for {transcript_id}: {e}")
            skipped_bad_seq += 1
            continue

        if strand == "-":
            seq = reverse_complement(seq)

        bed_start = window_start - 1
        bed_end = window_end
        name = f"{meta['transcript_id']}|{meta['gene_name']}"

        records.append({
            "gene_id": meta["gene_id"],
            "gene_name": meta["gene_name"],
            "transcript_id": meta["transcript_id"],
            "transcript_name": meta["transcript_name"],
            "chrom_gtf": chrom_gtf,
            "chrom_fasta": chrom_fasta,
            "strand": strand,
            "tis": tis,
            "window_start": window_start,
            "window_end": window_end,
            "bed_start": bed_start,
            "bed_end": bed_end,
            "upstream": UPSTREAM,
            "downstream": DOWNSTREAM,
            "seq": seq,
            "name": name
        })

    with open(output_tsv, "w", encoding="utf-8") as tsv, open(output_bed, "w", encoding="utf-8") as bed, open(output_fasta, "w", encoding="utf-8") as fasta:
        tsv.write("\t".join([
            "gene_id",
            "gene_name",
            "transcript_id",
            "transcript_name",
            "chrom_gtf",
            "chrom_fasta",
            "strand",
            "tis_genomic_pos",
            "window_start_1based",
            "window_end_1based",
            "bed_start_0based",
            "bed_end",
            "upstream",
            "downstream",
            "window_length"
        ]) + "\n")

        for rec in records:
            tsv.write("\t".join(map(str, [
                rec["gene_id"],
                rec["gene_name"],
                rec["transcript_id"],
                rec["transcript_name"],
                rec["chrom_gtf"],
                rec["chrom_fasta"],
                rec["strand"],
                rec["tis"],
                rec["window_start"],
                rec["window_end"],
                rec["bed_start"],
                rec["bed_end"],
                rec["upstream"],
                rec["downstream"],
                len(rec["seq"])
            ])) + "\n")

            bed.write("\t".join(map(str, [
                rec["chrom_fasta"],
                rec["bed_start"],
                rec["bed_end"],
                rec["name"],
                0,
                rec["strand"]
            ])) + "\n")

            fasta_header = f">{rec['gene_id']}|{rec['gene_name']}|{rec['transcript_id']}|{rec['transcript_name']}|{rec['chrom_fasta']}:{rec['window_start']}-{rec['window_end']}({rec['strand']})|TIS={rec['tis']}"
            fasta.write(fasta_header + "\n")
            fasta.write(wrap_fasta(rec["seq"]) + "\n")

    records_sorted = sorted(records, key=lambda x: (x["chrom_fasta"], x["bed_start"], x["bed_end"]))

    filtered_records = []
    last_chrom = None
    last_end = -1

    for rec in records_sorted:
        chrom = rec["chrom_fasta"]
        start = rec["bed_start"]
        end = rec["bed_end"]

        if chrom != last_chrom:
            filtered_records.append(rec)
            last_chrom = chrom
            last_end = end
        else:
            if start >= last_end:
                filtered_records.append(rec)
                last_end = end

    with open(output_no_overlap_tsv, "w", encoding="utf-8") as tsv, open(output_no_overlap_bed, "w", encoding="utf-8") as bed, open(output_no_overlap_fasta, "w", encoding="utf-8") as fasta:
        tsv.write("\t".join([
            "gene_id",
            "gene_name",
            "transcript_id",
            "transcript_name",
            "chrom_gtf",
            "chrom_fasta",
            "strand",
            "tis_genomic_pos",
            "window_start_1based",
            "window_end_1based",
            "bed_start_0based",
            "bed_end",
            "upstream",
            "downstream",
            "window_length"
        ]) + "\n")

        for rec in filtered_records:
            tsv.write("\t".join(map(str, [
                rec["gene_id"],
                rec["gene_name"],
                rec["transcript_id"],
                rec["transcript_name"],
                rec["chrom_gtf"],
                rec["chrom_fasta"],
                rec["strand"],
                rec["tis"],
                rec["window_start"],
                rec["window_end"],
                rec["bed_start"],
                rec["bed_end"],
                rec["upstream"],
                rec["downstream"],
                len(rec["seq"])
            ])) + "\n")

            bed.write("\t".join(map(str, [
                rec["chrom_fasta"],
                rec["bed_start"],
                rec["bed_end"],
                rec["name"],
                0,
                rec["strand"]
            ])) + "\n")

            fasta_header = f">{rec['gene_id']}|{rec['gene_name']}|{rec['transcript_id']}|{rec['transcript_name']}|{rec['chrom_fasta']}:{rec['window_start']}-{rec['window_end']}({rec['strand']})|TIS={rec['tis']}"
            fasta.write(fasta_header + "\n")
            fasta.write(wrap_fasta(rec["seq"]) + "\n")

    print(f"[INFO] Written entries (all): {len(records)}")
    print(f"[INFO] Written entries (no overlap): {len(filtered_records)}")
    print(f"[INFO] Skipped (no start_codon): {skipped_no_start}")
    print(f"[INFO] Skipped (bad/failed sequence): {skipped_bad_seq}")
    print(f"[INFO] TSV:             {output_tsv}")
    print(f"[INFO] BED:             {output_bed}")
    print(f"[INFO] FASTA:           {output_fasta}")
    print(f"[INFO] TSV no overlap:  {output_no_overlap_tsv}")
    print(f"[INFO] BED no overlap:  {output_no_overlap_bed}")
    print(f"[INFO] FASTA no overlap:{output_no_overlap_fasta}")


if __name__ == "__main__":
    main()