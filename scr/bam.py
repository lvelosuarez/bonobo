#!/usr/bin/env python3
"""
bam.py

Read BAM files from a folder, filter on mapping quality, and compare each
sample to a negative-control BAM **per reference** using z-scores and p-values,
then annotate references with GTDB taxonomy.

For each BAM:
- Only mapped reads with mapping_quality >= min_mapq are counted.
- Reads are counted separately for each reference (contig / chromosome),
  where the reference name is an accession (e.g. GCA_000027065.2).

If a negative-control BAM is provided:
    A single BAM is treated as the negative control.
    For each other BAM and each reference, we compute:

        mu = neg_read_count + pseudocount
        z  = (read_count - mu) / sqrt(mu)

    and a one-sided p-value:

        p = P(Z >= z) = 1 - Φ(z)

If NO negative-control BAM is provided:
    - neg_read_count is taken as 0 for all references
    - z and p are computed vs 0 + pseudocount.

Taxonomy:
- We load a GTDB taxonomy table (e.g. gtdb_taxonomy.rds) with at least:
    - accessions
    - gtdb_taxonomy
  and join it on reference = accessions.

Additional behaviour:
- Rows with read_count == 0 in the sample are **dropped** (even if the negative
  has reads for that accession), for all BAMs including the negative control.

Usage example with negative:
    python compare_bams_per_reference_to_negative_with_taxonomy.py \
        --bam-dir /path/to/bams \
        --neg-name NEGATIVE_SAMPLE.bam \
        --taxonomy /path/to/gtdb_taxonomy.rds \
        --output per_reference_sparki_compare.tsv \
        --min-mapq 25

Usage example without negative:
    python compare_bams_per_reference_to_negative_with_taxonomy.py \
        --bam-dir /path/to/bams \
        --taxonomy /path/to/gtdb_taxonomy.rds \
        --output per_reference_counts.tsv \
        --min-mapq 25
"""

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

import polars as pl
import pysam

# ---------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------


def normal_cdf(x: float) -> float:
    """Standard normal CDF Φ(x) using error function."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def z_and_p(count_sample: int, count_neg: int, pseudocount: float = 0.5):
    """
    Compute z-score and one-sided p-value (enrichment vs negative),
    allowing neg_read_count == 0 by adding a small pseudocount.

    Effective background mean:
        mu = count_neg + pseudocount

    z = (count_sample - mu) / sqrt(mu)
    p = 1 - Φ(z)  (one-sided, P(Z >= observed))

    Special case (both zero):
        z = 0, p = 1
    """
    if count_sample == 0 and count_neg == 0:
        return 0.0, 1.0

    mu = count_neg + pseudocount
    z = (count_sample - mu) / math.sqrt(mu)
    p = 1.0 - normal_cdf(z)
    return z, p


# ---------------------------------------------------------------------
# BAM processing
# ---------------------------------------------------------------------


def count_reads_per_reference(bam_path: Path, min_mapq: int):
    """
    Count paired-end reads per reference with mapping quality >= min_mapq.

    A read pair is counted ONLY if:
      - the read is paired
      - the pair is properly paired
      - both mates are mapped
      - both mates map to the SAME reference
      - this is read1 (to count the pair once)

    Returns:
        dict: {reference_name: paired_read_count}
    """
    counts = defaultdict(int)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            # Must be paired
            if not read.is_paired:
                continue

            # Count each pair only once
            if not read.is_read1:
                continue

            # Must be properly paired
            if not read.is_proper_pair:
                continue

            # Both mates must be mapped
            if read.is_unmapped or read.mate_is_unmapped:
                continue

            # MAPQ filter (apply to read1; mate usually similar)
            if read.mapping_quality < min_mapq:
                continue

            # Both mates must map to the same reference
            if read.reference_id != read.next_reference_id:
                continue

            ref_name = bam.get_reference_name(read.reference_id)
            counts[ref_name] += 1

    return dict(counts)


# ---------------------------------------------------------------------
# Taxonomy loading (polars)
# ---------------------------------------------------------------------


def load_taxonomy(taxonomy_path: Path) -> pl.DataFrame:
    """
    Load taxonomy table from .rds, .tsv, .txt or .csv into a Polars DataFrame.

    Expected columns:
      - accessions
      - gtdb_taxonomy
      - optionally: ncbi_taxid, gtdb_genome_representative, gtdb_representative
    """
    suffix = taxonomy_path.suffix.lower()

    if suffix == ".rds":
        try:
            import pyreadr  # type: ignore
        except ImportError as err:
            raise SystemExit("ERROR: pyreadr is required to read .rds files.\n" "Install it with: pip install pyreadr") from err

        result = pyreadr.read_r(str(taxonomy_path))
        if not result:
            raise SystemExit(f"ERROR: No object found in RDS file: {taxonomy_path}")
        pandas_df = next(iter(result.values()))
        tax_df = pl.from_pandas(pandas_df)

    elif suffix in {".tsv", ".txt"}:
        tax_df = pl.read_csv(taxonomy_path, separator="\t")

    elif suffix == ".csv":
        tax_df = pl.read_csv(taxonomy_path)

    else:
        raise SystemExit(f"ERROR: Unsupported taxonomy file extension: {suffix}\n" f"Use .rds, .tsv, .txt or .csv")

    if "accessions" not in tax_df.columns:
        raise SystemExit("ERROR: taxonomy table must contain a column named 'accessions'.")

    keep_cols = [
        c
        for c in [
            "accessions",
            "gtdb_taxonomy",
            "ncbi_taxid",
            "gtdb_genome_representative",
            "gtdb_representative",
        ]
        if c in tax_df.columns
    ]

    return tax_df.select(keep_cols).unique(subset=["accessions"])


# ---------------------------------------------------------------------
# File / sample helpers
# ---------------------------------------------------------------------


def sample_num_from_bam(path: Path) -> int:
    """
    Extract trailing number from BAM filename for numeric sorting.

    Example:
        MG-Run35-Ech-12.bam -> 12
        MG-Run35-Ech-2.bam  -> 2
    """
    m = re.search(r"(\d+)$", path.stem)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare BAM files to an optional negative sample using MAPQ-filtered "
            "paired read counts per reference, and annotate with GTDB taxonomy."
        )
    )
    p.add_argument(
        "--bam-dir",
        required=True,
        help="Directory containing BAM files.",
    )
    p.add_argument(
        "--neg-name",
        required=False,
        default=None,
        help=(
            "Optional file name of the negative-control BAM within bam-dir "
            "(e.g. 'MG-Run35-Ech-11.bam'). If omitted, background is 0 "
            "for all references."
        ),
    )
    p.add_argument(
        "--taxonomy",
        required=True,
        help=("GTDB taxonomy file (.rds, .tsv, .txt or .csv) with at least " "an 'accessions' column and ideally 'gtdb_taxonomy'."),
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=25,
        help="Minimum mapping quality to keep a read pair (default: 25).",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output TSV file with per-reference counts, z-scores, p-values and taxonomy.",
    )
    p.add_argument(
        "--ignore-samples",
        help=("Optional text file with one sample name per line to ignore " "(either 'MG-Run35-Ech-12' or 'MG-Run35-Ech-12.bam')."),
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    args = parse_args()

    bam_dir = Path(args.bam_dir)
    if not bam_dir.is_dir():
        raise SystemExit(f"bam-dir is not a directory: {bam_dir}")

    taxonomy_path = Path(args.taxonomy)
    if not taxonomy_path.exists():
        raise SystemExit(f"Taxonomy file not found: {taxonomy_path}")

    # Load ignore list (normalise to stems, without .bam)
    ignore_samples: set[str] = set()
    if args.ignore_samples:
        ignore_path = Path(args.ignore_samples)
        if not ignore_path.exists():
            raise SystemExit(f"Ignore-samples file not found: {ignore_path}")
        with ignore_path.open() as f:
            for line in f:
                name = line.strip()
                if not name:
                    continue
                # strip optional .bam suffix
                if name.endswith(".bam"):
                    name = name[:-4]
                ignore_samples.add(name)
        print(f"[INFO] Ignoring {len(ignore_samples)} samples: " f"{', '.join(sorted(ignore_samples))}")

    # Collect BAMs and sort numerically by sample suffix
    bam_files = sorted(
        bam_dir.glob("*.bam"),
        key=sample_num_from_bam,
    )
    if not bam_files:
        raise SystemExit(f"No BAM files found in directory: {bam_dir}")

    # Negative-control handling (optional)
    neg_path: Path | None = None
    neg_counts: dict[str, int] = {}

    if args.neg_name:
        neg_path = bam_dir / args.neg_name
        if not neg_path.exists():
            raise SystemExit(f"Negative-control BAM not found: {neg_path}\n" f"Make sure --neg-name matches a BAM in {bam_dir}")
        print(f"[INFO] Using negative control BAM: {neg_path.name}")

        # Per-reference counts for negative control
        print(f"[INFO] Counting paired reads per reference for negative control " f"(MAPQ >= {args.min_mapq})...")
        neg_counts = count_reads_per_reference(neg_path, args.min_mapq)
        print(f"[INFO] Negative control references: {len(neg_counts)}")
    else:
        print("[INFO] No negative-control BAM provided. " "Using 0 background counts for all references.")
        neg_counts = {}

    # Load taxonomy
    print(f"[INFO] Loading taxonomy from: {taxonomy_path}")
    tax_df = load_taxonomy(taxonomy_path)
    print(f"[INFO] Taxonomy table loaded with {tax_df.height} accessions; " f"columns: {', '.join(tax_df.columns)}")

    # Process all BAMs (skip negative in output if defined)
    rows: list[dict] = []

    for bam_path in bam_files:
        sample = bam_path.stem

        # Skip the negative control from per-sample output if present
        if neg_path is not None and bam_path == neg_path:
            print(f"[INFO] Skipping negative control sample in output: {sample}")
            continue

        # Skip ignored samples
        if sample in ignore_samples:
            print(f"[INFO] Skipping ignored sample: {sample}")
            continue

        print(f"[INFO] Processing sample: {sample} " f"(MAPQ >= {args.min_mapq})")
        sample_counts = count_reads_per_reference(bam_path, args.min_mapq)

        # Union of references from this sample and the negative control (if any)
        references = set(sample_counts.keys()) | set(neg_counts.keys())

        for ref in references:
            count_sample = sample_counts.get(ref, 0)
            count_neg = neg_counts.get(ref, 0)

            # Eliminate references with 0 reads in the sample
            if count_sample == 0:
                continue

            z, p = z_and_p(count_sample, count_neg)

            rows.append(
                {
                    "sample": sample,
                    # "bam_path": str(bam_path),
                    "reference": ref,  # accession ID (e.g. GCA_...)
                    "mapq_min": args.min_mapq,
                    "read_count": count_sample,
                    "neg_read_count": count_neg,
                    "z_score": z,
                    "p_value": p,
                }
            )

    if not rows:
        raise SystemExit("No rows produced (all samples empty or ignored).")

    df = pl.DataFrame(rows)

    # Merge taxonomy: reference (accession) -> GTDB taxonomy
    df = df.join(
        tax_df,
        how="left",
        left_on="reference",
        right_on="accessions",
    )

    # Add numeric sample index for sorting: MG-Run35-Ech-12 -> 12
    df = df.with_columns(pl.col("sample").str.extract(r"(\d+)$").cast(pl.Int32).alias("sample_num"))

    # Sort: by numeric sample order, then read_count descending
    df = (
        df.sort(
            by=["sample_num", "read_count"],
            descending=[False, True],
        )
        .drop("sample_num")
        .drop("gtdb_genome_representative")
        .drop("gtdb_representative")
    )

    # Save
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(out_path, separator="\t")
    print(f"[INFO] Written results to {out_path}")


if __name__ == "__main__":
    main()
