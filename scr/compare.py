#!/usr/bin/env python3
"""
compare.py

Merge:
- CUB output (Kraken2/minimizer stats) with columns including:
    sample, taxid, name, z_score, p_value, p_value_neg, ...
- bam.py output with columns including:
    sample, reference, read_count, neg_read_count, z_score, p_value,
    gtdb_taxonomy, ncbi_taxid

Join key: (sample, NCBI taxid)  ->  cub.taxid  <->  bam.ncbi_taxid

Output: one TSV with
- sample
- taxid
- name (from CUB, when available)
- gtdb_taxonomy (from BAM, when available)
- all CUB metrics prefixed with `cub_`
- all BAM metrics prefixed with `bam_`

Rows are sorted:
- by sample "prefix" (non-numeric part) then numeric suffix, i.e.
  MG-Run35-Ech-1, MG-Run35-Ech-2, ..., MG-Run35-Ech-10
- within each sample, by bam_read_count (descending)
"""

import argparse
from pathlib import Path

import polars as pl

# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(description="Compare CUB TSV with bam.py TSV using sample + NCBI taxid (polars).")
    p.add_argument("--cub", required=True, help="CUB TSV file")
    p.add_argument("--bam", required=True, help="bam.py TSV file")
    p.add_argument("--output", required=True, help="Output merged TSV")

    # Column names – adjust if your files differ
    p.add_argument(
        "--cub-sample-col",
        default="sample",
        help="Sample column in CUB TSV (default: 'sample')",
    )
    p.add_argument(
        "--cub-taxid-col",
        default="taxid",
        help="NCBI taxid column in CUB TSV (default: 'taxid')",
    )
    p.add_argument(
        "--bam-sample-col",
        default="sample",
        help="Sample column in bam TSV (default: 'sample')",
    )
    p.add_argument(
        "--bam-taxid-col",
        default="ncbi_taxid",
        help="NCBI taxid column in bam TSV (default: 'ncbi_taxid')",
    )

    return p.parse_args()


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def clean_taxid(df: pl.DataFrame, col: str, out_col: str = "taxid") -> pl.DataFrame:
    """
    Normalize taxid / ncbi_taxid so that:
    - '43769', '43769.0', 43769, 43769.0 -> 43769 (Int64)
    We store as Int64 for join.
    """
    if col not in df.columns:
        raise SystemExit(f"Column '{col}' not found in dataframe when cleaning taxid.")

    return df.with_columns(
        pl.col(col)
        .cast(pl.Float64, strict=False)  # strings like '43769.0' -> float
        .round(0)  # 43769.0 -> 43769.0
        .cast(pl.Int64, strict=False)  # -> 43769 (Int64)
        .alias(out_col)
    )


def load_cub(path: str, sample_col: str, taxid_col: str) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", infer_schema_length=10000)

    if sample_col not in df.columns:
        raise SystemExit(f"CUB: sample column '{sample_col}' not found in {path}")
    if taxid_col not in df.columns:
        raise SystemExit(f"CUB: taxid column '{taxid_col}' not found in {path}")

    # Normalize taxid
    df = clean_taxid(df, taxid_col, out_col="taxid")

    # Rename sample column to 'sample' if needed
    if sample_col != "sample":
        df = df.rename({sample_col: "sample"})

    # Ensure sample is string
    df = df.with_columns(pl.col("sample").cast(pl.Utf8))

    # Key columns
    key_cols = ["sample", "taxid"]
    if "name" in df.columns:
        key_cols.append("name")

    # Metric columns = everything else
    metric_cols = [c for c in df.columns if c not in key_cols]

    # Order columns: keys first, then metrics
    df = df.select(key_cols + metric_cols)

    # Prefix metric columns with cub_
    rename_map = {c: f"cub_{c}" for c in metric_cols}
    df = df.rename(rename_map)

    return df


def load_bam(path: str, sample_col: str, taxid_col: str) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", infer_schema_length=10000)

    if sample_col not in df.columns:
        raise SystemExit(f"BAM: sample column '{sample_col}' not found in {path}")
    if taxid_col not in df.columns:
        raise SystemExit(f"BAM: taxid column '{taxid_col}' not found in {path}")

    # Normalize taxid
    df = clean_taxid(df, taxid_col, out_col="taxid")

    # Rename sample column to 'sample' if needed
    if sample_col != "sample":
        df = df.rename({sample_col: "sample"})

    # Ensure sample is string
    df = df.with_columns(pl.col("sample").cast(pl.Utf8))

    # Key columns
    key_cols = ["sample", "taxid"]
    if "gtdb_taxonomy" in df.columns:
        key_cols.append("gtdb_taxonomy")

    # Metric columns = everything else
    metric_cols = [c for c in df.columns if c not in key_cols]

    # Order columns: keys first, then metrics
    df = df.select(key_cols + metric_cols)

    # Prefix metric columns with bam_
    rename_map = {c: f"bam_{c}" for c in metric_cols}
    df = df.rename(rename_map)

    return df


def natural_sort(merged: pl.DataFrame) -> pl.DataFrame:
    """
    Sort the merged dataframe by:
    - sample prefix (non-numeric end)
    - numeric suffix of sample (if present)
    - bam_read_count (descending)

    Example:
      MG-Run35-Ech-1
      MG-Run35-Ech-2
      MG-Run35-Ech-10
    """

    # Add columns to help sorting
    # Extract numeric suffix; if no match -> null
    merged = merged.with_columns(
        [
            pl.col("sample").cast(pl.Utf8).str.extract(r"(\d+)$", 1).cast(pl.Int64, strict=False).alias("sample_num"),
            pl.col("sample").cast(pl.Utf8).str.replace(r"\d+$", "", literal=False).alias("sample_prefix"),
        ]
    )

    # Sorting by bam_read_count (if present), else use 0
    if "bam_read_count" in merged.columns:
        merged = merged.with_columns(pl.col("bam_read_count").fill_null(0).alias("_sort_reads"))
    else:
        merged = merged.with_columns(pl.lit(0).alias("_sort_reads"))

    merged = merged.sort(
        by=["sample_prefix", "sample_num", "_sort_reads"],
        descending=[False, False, True],
    )

    # Drop helper columns
    merged = merged.drop(["sample_prefix", "sample_num", "_sort_reads"])

    return merged


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main():
    args = parse_args()

    cub_df = load_cub(
        args.cub,
        sample_col=args.cub_sample_col,
        taxid_col=args.cub_taxid_col,
    )

    bam_df = load_bam(
        args.bam,
        sample_col=args.bam_sample_col,
        taxid_col=args.bam_taxid_col,
    )

    # Outer join on (sample, taxid)
    merged = cub_df.join(bam_df, on=["sample", "taxid"], how="full")

    # Natural-ish sort
    merged = natural_sort(merged)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write as TSV
    merged.write_csv(out_path, separator="\t")


if __name__ == "__main__":
    main()
