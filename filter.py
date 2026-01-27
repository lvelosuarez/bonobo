#!/usr/bin/env python3
"""
summarize_cub.py

Create a clinician-friendly summary from a cub TSV.

Features:
- Uses negative-control corrected stats if present (fdr_neg, p_value_neg).
- Filters by p_value and genome coverage.
- Always keeps taxa from a pathogen panel (species-level).
- Drops host taxa.
- Outputs a compact, readable TSV.

Input:  *.tsv
Output: summary.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize cub results into a table."
    )
    p.add_argument("--input", required=True, help="cub TSV")
    p.add_argument("--output", required=True, help="Output TSV")
    p.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Top N taxa per sample (default: keep all)",
    )
    p.add_argument("--max-p", type=float, default=0.1, help="Max p_value to keep taxa")
    p.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Min genome_coverage_estimate (0–1)",
    )
    p.add_argument(
        "--pathogen-panel",
        help="TSV/CSV with columns: taxid and/or organism_name",
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def choose_stats_columns(df):
    """Prefer negative-control stats if present."""
    if {"fdr_neg", "p_value_neg", "z_score_neg"}.issubset(df.columns):
        return "fdr_neg", "p_value_neg", "z_score_neg"
    return "fdr", "p_value", "z_score"


def load_pathogen_panel(path):
    """
    Load pathogen panel.
    Accepts columns:
      - taxid
      - organism_name (or organism / name)
    """
    if path is None:
        return set(), set()

    panel = pd.read_csv(path, sep=None, engine="python")

    taxids = set()
    names = set()

    if "taxid" in panel.columns:
        taxids = {str(x) for x in panel["taxid"].dropna().astype(str)}

    for col in ["organism_name", "organism", "name"]:
        if col in panel.columns:
            names |= {
                str(x).strip().lower()
                for x in panel[col].dropna().astype(str)
            }

    if not taxids and not names:
        raise RuntimeError(
            "Pathogen panel must contain 'taxid' or 'organism_name'"
        )

    return taxids, names


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # 1) Drop host taxa
    if "is_host" in df.columns:
        df = df[~df["is_host"]].copy()

    # 2) Choose stats columns (negative-control vs global)
    fdr_col, p_col, z_col = choose_stats_columns(df)

    # 3) Ensure coverage column exists
    if "genome_coverage_estimate" not in df.columns:
        df["genome_coverage_estimate"] = np.nan

    # 4) Load pathogen panel
    panel_taxids, panel_names = load_pathogen_panel(args.pathogen_panel)

    df["on_pathogen_panel"] = False
    if panel_taxids:
        df.loc[df["taxid"].astype(str).isin(panel_taxids), "on_pathogen_panel"] = True
    if panel_names:
        df.loc[
            df["name"].str.lower().isin(panel_names),
            "on_pathogen_panel",
        ] = True

    # 5) Automatic filtering (stats-based, using p-value)
    df_auto = df.copy()
    df_auto = df_auto[df_auto[p_col].notna()]
    df_auto = df_auto[df_auto[p_col] <= args.max_p]
    df_auto = df_auto[
        df_auto["genome_coverage_estimate"].fillna(0.0) >= args.min_coverage
    ]

    # 6) Force-keep pathogen panel taxa
    df_panel = df[df["on_pathogen_panel"]].copy()

    df = (
        pd.concat([df_auto, df_panel], ignore_index=True)
        .drop_duplicates(subset=["sample", "taxid"])
    )

    if df.empty:
        print("No taxa to report.")
        return

    # 7) Add clinician-friendly columns (no confidence categories)
    df["relative_abundance"] = df["minimizer_proportion_sample"]
    df["genome_coverage_percent"] = df["genome_coverage_estimate"] * 100.0

    df["fdr_used"] = df[fdr_col]
    df["p_value_used"] = df[p_col]
    df["z_score_used"] = df[z_col]

    # 8) Ordering: samples by numeric index, taxa per sample by number of reads
    #    - samples like 'sample1', 'sample2', ..., 'sample10' will be ordered 1,2,3,...,10
    s = df["sample"].astype(str)
    num = s.str.extract(r"(\d+)$", expand=False)
    df["_sample_num"] = pd.to_numeric(num, errors="coerce")

    # Sort by numeric sample index, then sample name, then reads (descending)
    df = df.sort_values(
        ["_sample_num", "sample", "z_score"],
        ascending=[True, True, False],
    )

    # 9) Optional per-sample top N (default: keep all)
    if args.top_n is not None:
        df = df.groupby("sample", group_keys=False).head(args.top_n)

    # 10) Final clinician table
    out = df[
        [
            "sample",
            "name",
            "rank_code",
            "taxid",
            "fragments_clade",
            "relative_abundance",
            "genome_coverage_percent",
            "z_score_used",
            "p_value_used",
            "fdr_used",
            "on_pathogen_panel",
            "_sample_num",
        ]
    ].copy()

    out = out.rename(
        columns={
            "name": "organism",
            "rank_code": "rank",
            "fragments_clade": "reads_in_clade",
            "relative_abundance": "relative_abundance_minimizers",
            "genome_coverage_percent": "genome_coverage_percent_minimizers",
            "z_score_used": "z_score",
            "p_value_used": "p_value",
            "fdr_used": "fdr",
        }
    )

    # Drop internal sort key from final output
    out = out.drop(columns=["_sample_num"])

    # 11) Rounding for readability
    out["relative_abundance_minimizers"] = out[
        "relative_abundance_minimizers"
    ].round(10)
    out["genome_coverage_percent_minimizers"] = out[
        "genome_coverage_percent_minimizers"
    ].round(5)
    out["z_score"] = out["z_score"].round(2)
    out["p_value"] = out["p_value"].apply(lambda x: float(f"{x:.3g}"))
    out["fdr"] = out["fdr"].apply(lambda x: float(f"{x:.3g}"))

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote summary to {args.output}")


if __name__ == "__main__":
    main()
