#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import os
import re
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================

genome_file = "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

OUT_DIR = Path("output_orfs_parallel")
TMP_DIR = OUT_DIR / "tmp_by_chrom"

OUT_TSV = OUT_DIR / "human_GRCh38_all_ORFs.tsv.gz"
OUT_HIST = OUT_DIR / "human_GRCh38_ORF_lengths_histogram.png"

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

MIN_ORF_LENGTH = 0
MAX_ORF_LENGTH = None

BIN_SIZE = 30

RESET_ON_AMBIGUOUS = True

# Βάλε π.χ. 4, 8, 12.
# Αν το αφήσεις None, χρησιμοποιεί όλα τα cores - 1.
N_WORKERS = None

PRINT_PROGRESS = True
PROGRESS_EVERY_ORFS = 100000


# ============================================================
# HELPERS
# ============================================================

def open_text_maybe_gzip(path, mode="wt"):
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")

    return open(path, mode, encoding="utf-8", newline="")


def read_fasta_chromosomes(fasta_file):
    """
    Streaming FASTA reader.
    Yields one chromosome/contig at a time:
        chrom, sequence
    """
    chrom = None
    seq_parts = []

    with open(fasta_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(seq_parts).upper()

                chrom = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if chrom is not None:
        yield chrom, "".join(seq_parts).upper()


def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def is_valid_codon(codon):
    return len(codon) == 3 and all(base in "ACGT" for base in codon)


def update_histogram(hist_counts, length):
    bin_start = (length // BIN_SIZE) * BIN_SIZE
    hist_counts[bin_start] += 1


def passes_length_filter(length):
    if length < MIN_ORF_LENGTH:
        return False

    if MAX_ORF_LENGTH is not None and length > MAX_ORF_LENGTH:
        return False

    return True


def safe_chrom_name(chrom):
    """
    Κάνει ασφαλές filename από chromosome name.
    """
    return re.sub(r"[^A-Za-z0-9_.-]", "_", chrom)


def scan_sequence_for_orfs(
    seq,
    chrom,
    strand,
    chrom_len,
    writer,
    hist_counts,
    stats
):
    """
    Finds all candidate ORFs:
        ATG ... in-frame STOP

    Για + strand:
        coordinates είναι direct.

    Για - strand:
        seq είναι reverse complement και τα coordinates
        μετατρέπονται πίσω στο original chromosome.
    """

    n = len(seq)

    for frame in range(3):
        open_starts = []
        i = frame

        while i + 3 <= n:
            codon = seq[i:i + 3]

            if RESET_ON_AMBIGUOUS and not is_valid_codon(codon):
                open_starts = []
                i += 3
                continue

            if codon == START_CODON:
                open_starts.append(i)

            elif codon in STOP_CODONS:
                if open_starts:
                    orf_end = i + 3

                    for orf_start in open_starts:
                        length = orf_end - orf_start

                        if not passes_length_filter(length):
                            continue

                        orf_seq = seq[orf_start:orf_end]

                        if strand == "+":
                            genomic_start = orf_start
                            genomic_stop = orf_end
                        else:
                            genomic_start = chrom_len - orf_end
                            genomic_stop = chrom_len - orf_start

                        writer.writerow({
                            "chrom": chrom,
                            "start": genomic_start,
                            "stop": genomic_stop,
                            "strand": strand,
                            "length": length,
                            "sequence": orf_seq,
                        })

                        update_histogram(hist_counts, length)

                        stats["total_orfs"] += 1
                        stats["orfs_by_strand"][strand] += 1

                        if (
                            PRINT_PROGRESS
                            and stats["total_orfs"] % PROGRESS_EVERY_ORFS == 0
                        ):
                            print(
                                f"[{chrom} {strand}] "
                                f"{stats['total_orfs']:,} ORFs"
                            )

                    # Μετά από stop codon κλείνουν όλα τα ORFs αυτού του frame
                    open_starts = []

            i += 3


def process_one_chromosome(args):
    """
    Worker function.
    Κάθε process επεξεργάζεται ένα chromosome και γράφει δικό του temp TSV.gz.
    """
    chrom, seq = args

    chrom_len = len(seq)

    tmp_file = TMP_DIR / f"{safe_chrom_name(chrom)}.orfs.tsv.gz"

    columns = [
        "chrom",
        "start",
        "stop",
        "strand",
        "length",
        "sequence",
    ]

    hist_counts = defaultdict(int)

    stats = {
        "chrom": chrom,
        "chrom_len": chrom_len,
        "total_orfs": 0,
        "orfs_by_strand": {
            "+": 0,
            "-": 0,
        },
        "tmp_file": str(tmp_file),
    }

    if PRINT_PROGRESS:
        print(f"Processing {chrom} length={chrom_len:,}")

    with open_text_maybe_gzip(tmp_file, "wt") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=columns,
            delimiter="\t"
        )
        writer.writeheader()

        # + strand
        scan_sequence_for_orfs(
            seq=seq,
            chrom=chrom,
            strand="+",
            chrom_len=chrom_len,
            writer=writer,
            hist_counts=hist_counts,
            stats=stats
        )

        # - strand
        rc_seq = reverse_complement(seq)

        scan_sequence_for_orfs(
            seq=rc_seq,
            chrom=chrom,
            strand="-",
            chrom_len=chrom_len,
            writer=writer,
            hist_counts=hist_counts,
            stats=stats
        )

    stats["hist_counts"] = dict(hist_counts)

    if PRINT_PROGRESS:
        print(
            f"Finished {chrom}: "
            f"{stats['total_orfs']:,} ORFs "
            f"(+ {stats['orfs_by_strand']['+']:,}, "
            f"- {stats['orfs_by_strand']['-']:,})"
        )

    return stats


def merge_temp_files(results, out_tsv):
    """
    Ενώνει όλα τα chromosome-level temp files σε ένα τελικό TSV.gz.
    """
    columns = [
        "chrom",
        "start",
        "stop",
        "strand",
        "length",
        "sequence",
    ]

    with open_text_maybe_gzip(out_tsv, "wt") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(columns)

        for r in results:
            tmp_file = r["tmp_file"]

            with open_text_maybe_gzip(tmp_file, "rt") as inp:
                reader = csv.reader(inp, delimiter="\t")

                # skip header
                next(reader, None)

                for row in reader:
                    writer.writerow(row)


def plot_histogram(hist_counts, out_png):
    if not hist_counts:
        print("No ORFs found. Histogram was not created.")
        return

    x = sorted(hist_counts.keys())
    y = [hist_counts[k] for k in x]

    plt.figure(figsize=(12, 6))
    plt.bar(x, y, width=BIN_SIZE, align="edge")

    plt.xlabel("ORF length (nt)")
    plt.ylabel("Number of ORFs")
    plt.title("Distribution of ORF lengths")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def remove_temp_files(results):
    for r in results:
        tmp_file = Path(r["tmp_file"])
        if tmp_file.exists():
            tmp_file.unlink()


# ============================================================
# MAIN
# ============================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    if N_WORKERS is None:
        workers = max(1, cpu_count() - 1)
    else:
        workers = N_WORKERS

    print("Scanning genome for ORFs in parallel")
    print(f"Genome: {genome_file}")
    print(f"Workers: {workers}")
    print(f"Output TSV: {OUT_TSV}")
    print(f"Histogram: {OUT_HIST}")
    print(f"MIN_ORF_LENGTH: {MIN_ORF_LENGTH}")
    print(f"MAX_ORF_LENGTH: {MAX_ORF_LENGTH}")
    print()

    results = []

    with Pool(processes=workers) as pool:
        for result in pool.imap_unordered(
            process_one_chromosome,
            read_fasta_chromosomes(genome_file)
        ):
            results.append(result)

    # Για πιο σταθερό output order: 1,2,3,...,X,Y αν γίνεται
    results.sort(key=lambda r: r["chrom"])

    print()
    print("Merging chromosome files...")
    merge_temp_files(results, OUT_TSV)

    print("Combining histograms...")
    final_hist_counts = Counter()

    total_orfs = 0
    total_plus = 0
    total_minus = 0

    for r in results:
        final_hist_counts.update(r["hist_counts"])
        total_orfs += r["total_orfs"]
        total_plus += r["orfs_by_strand"]["+"]
        total_minus += r["orfs_by_strand"]["-"]

    plot_histogram(final_hist_counts, OUT_HIST)

    print("Removing temporary files...")
    remove_temp_files(results)

    print()
    print("Done.")
    print(f"Chromosomes/contigs processed: {len(results):,}")
    print(f"Total ORFs: {total_orfs:,}")
    print(f"+ strand ORFs: {total_plus:,}")
    print(f"- strand ORFs: {total_minus:,}")
    print(f"TSV file: {OUT_TSV}")
    print(f"Histogram: {OUT_HIST}")


if __name__ == "__main__":
    main()