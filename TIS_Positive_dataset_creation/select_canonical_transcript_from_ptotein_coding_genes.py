#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import re
import csv
from pathlib import Path

###############################################################################
# ΡΥΘΜΙΣΕΙΣ
###############################################################################

# Βάλε εδώ το GTF της Ensembl που θέλεις.
# Παράδειγμα για άνθρωπο:
GTF_FILE = "../data/ensembl/Homo_sapiens.GRCh38.116.gtf.gz"

# Αρχείο εξόδου TSV
OUTPUT_TSV = "ensembl_protein_coding_canonical_transcripts.tsv"


###############################################################################
# ΒΟΗΘΗΤΙΚΕΣ ΣΥΝΑΡΤΗΣΕΙΣ
###############################################################################

def parse_gtf_attributes(attr_text: str) -> dict:
    """
    Κάνει parse το 9ο πεδίο του GTF και επιστρέφει dictionary.
    Υποστηρίζει και πολλαπλά tags.
    """
    attrs = {}
    tags = []

    # GTF μορφή: key "value";
    for match in re.finditer(r'(\S+)\s+"([^"]*)";', attr_text):
        key, value = match.group(1), match.group(2)
        if key == "tag":
            tags.append(value)
        else:
            attrs[key] = value

    attrs["tag"] = tags
    return attrs


def open_maybe_gzip(filepath: str):
    """
    Ανοίγει κανονικό ή gzipped αρχείο.
    """
    path = Path(filepath)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


###############################################################################
# ΚΥΡΙΑ ΛΟΓΙΚΗ
###############################################################################

def extract_protein_coding_canonical_transcripts(gtf_file: str, output_tsv: str):
    """
    Διαβάζει Ensembl GTF και κρατά μόνο canonical protein-coding transcripts.
    """

    kept = 0
    total_transcripts = 0

    with open_maybe_gzip(gtf_file) as fin, open(output_tsv, "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout, delimiter="\t")

        writer.writerow([
            "gene_id",
            "gene_name",
            "gene_biotype",
            "transcript_id",
            "transcript_name",
            "transcript_biotype",
            "chromosome",
            "start",
            "end",
            "strand",
            "is_ensembl_canonical",
            "tags"
        ])

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            # Θέλουμε μόνο transcript γραμμές
            if feature != "transcript":
                continue

            total_transcripts += 1

            attrs = parse_gtf_attributes(attributes)

            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            gene_biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            transcript_id = attrs.get("transcript_id", "")
            transcript_name = attrs.get("transcript_name", "")
            transcript_biotype = attrs.get("transcript_biotype", attrs.get("transcript_type", ""))
            tags = attrs.get("tag", [])

            # Φιλτράρουμε μόνο protein-coding genes / transcripts
            if gene_biotype != "protein_coding":
                continue
            if transcript_biotype != "protein_coding":
                continue

            # Κρατάμε μόνο canonical transcripts
            is_canonical = "Ensembl_canonical" in tags
            if not is_canonical:
                continue

            writer.writerow([
                gene_id,
                gene_name,
                gene_biotype,
                transcript_id,
                transcript_name,
                transcript_biotype,
                chrom,
                start,
                end,
                strand,
                1,
                ",".join(tags)
            ])

            kept += 1

    print(f"Συνολικά transcript records: {total_transcripts}")
    print(f"Canonical protein-coding transcripts που κρατήθηκαν: {kept}")
    print(f"Αποθηκεύτηκαν στο: {output_tsv}")


if __name__ == "__main__":
    extract_protein_coding_canonical_transcripts(GTF_FILE, OUTPUT_TSV)