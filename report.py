#!/usr/bin/env python3
"""
report.py

Generate an HTML report from SPARKI-like TSV output.

Uses NEGATIVE-CONTROL P-VALUES (p_value_neg) when available.

Filtering, ranking and confidence are based on p-values (NOT FDR).

Dependencies:
    - pandas
    - numpy
"""

import argparse
from pathlib import Path
import html

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate an HTML report from SPARKI-like results."
    )
    p.add_argument("--input", required=True, help="SPARKI-like TSV file")
    p.add_argument("--output", required=True, help="Output HTML report")
    p.add_argument("--top-n", type=int, default=500)
    p.add_argument(
        "--max-pvalue",
        type=float,
        default=0.05,
        help="Maximum p-value to keep taxa (default: 0.05)",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Minimum genome_coverage_estimate (0–1)",
    )
    p.add_argument("--pathogen-panel", help="Pathogen panel TSV/CSV")
    return p.parse_args()


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def choose_stats_columns(df):
    """
    Prefer negative-control p-values if present.
    """
    if "p_value_neg" in df.columns:
        return {
            "p": "p_value_neg",
            "z": "z_score_neg" if "z_score_neg" in df.columns else None,
            "fdr": "fdr_neg" if "fdr_neg" in df.columns else None,
        }
    return {
        "p": "p_value",
        "z": "z_score",
        "fdr": "fdr",
    }


def load_pathogen_panel(path):
    if path is None:
        return set(), set()

    panel = pd.read_csv(path, sep=None, engine="python")

    taxids = set(panel["taxid"].astype(str)) if "taxid" in panel.columns else set()
    names = set()

    for col in ["organism_name", "organism", "name"]:
        if col in panel.columns:
            names |= {
                str(x).strip().lower()
                for x in panel[col].dropna().astype(str)
            }

    return taxids, names


def add_confidence(df, p_col):
    """
    Confidence based on p-value + coverage.
    """
    p = df[p_col]
    cov = df["genome_coverage_estimate"].fillna(0.0)

    conditions = [
        (p <= 1e-6) & (cov >= 0.20),
        (p <= 1e-4) & (cov >= 0.05),
        (p <= 0.01) & (cov >= 0.01),
    ]
    choices = ["High", "Medium", "Low"]

    df["confidence"] = np.select(
        conditions,
        choices,
        default="Very low / background",
    )
    return df


def confidence_class(conf):
    return {
        "High": "conf-high",
        "Medium": "conf-med",
        "Low": "conf-low",
    }.get(conf, "conf-verylow")


# ---------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------

def make_html_report(df, output_path):
    css = """
    <style>
    body { font-family: sans-serif; margin: 1.5rem; }
    table { border-collapse: collapse; width: 100%; }
    th, td { border: 1px solid #ddd; padding: 0.3rem; font-size: 0.85rem; }
    th { background-color: #eee; }
    .pathogen { background-color: #fff4e0; font-weight: bold; }
    .conf-high { background-color: #d9fdd3; }
    .conf-med { background-color: #fff5cc; }
    .conf-low { background-color: #ffe4d6; }
    .conf-verylow { color: #777; }
    </style>
    """

    html_out = [
        "<!DOCTYPE html>",
        "<html><head><meta charset='utf-8'>",
        "<title>Metagenomic Report</title>",
        css,
        "</head><body>",
        "<h1>Metagenomic Report</h1>",
    ]

    for sample, sub in df.groupby("sample"):
        html_out.append(f"<h2>Sample: {html.escape(str(sample))}</h2>")
        html_out.append("<table>")
        html_out.append(
            "<tr>"
            "<th>Organism</th><th>Rank</th><th>Reads</th>"
            "<th>Coverage %</th><th>p-value</th><th>Confidence</th><th>Panel</th>"
            "</tr>"
        )

        for _, r in sub.iterrows():
            cls = confidence_class(r["confidence"])
            if r["on_pathogen_panel"]:
                cls += " pathogen"

            html_out.append(
                f"<tr class='{cls}'>"
                f"<td>{html.escape(str(r['organism']))}</td>"
                f"<td>{r['rank']}</td>"
                f"<td>{int(r['reads_in_clade'])}</td>"
                f"<td>{r['genome_coverage_percent_minimizers']:.1f}</td>"
                f"<td>{r['p_value']:.3g}</td>"
                f"<td>{r['confidence']}</td>"
                f"<td>{'Yes' if r['on_pathogen_panel'] else ''}</td>"
                "</tr>"
            )

        html_out.append("</table>")

    html_out.append("</body></html>")
    output_path.write_text("\n".join(html_out), encoding="utf-8")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Drop host
    if "is_host" in df.columns:
        df = df[~df["is_host"]]

    stats = choose_stats_columns(df)
    p_col = stats["p"]

    if "genome_coverage_estimate" not in df.columns:
        df["genome_coverage_estimate"] = 0.0

    # Pathogen panel
    taxids, names = load_pathogen_panel(args.pathogen_panel)
    df["on_pathogen_panel"] = False
    if taxids:
        df.loc[df["taxid"].astype(str).isin(taxids), "on_pathogen_panel"] = True
    if names:
        df.loc[df["name"].str.lower().isin(names), "on_pathogen_panel"] = True

    # Filter by p-value + coverage
    df_auto = df[
        (df[p_col].notna()) &
        (df[p_col] <= args.max_pvalue) &
        (df["genome_coverage_estimate"] >= args.min_coverage)
    ]

    # Force keep pathogen taxa
    df_panel = df[df["on_pathogen_panel"]]
    df = pd.concat([df_auto, df_panel]).drop_duplicates(["sample", "taxid"])

    if df.empty:
        print("No taxa to report.")
        return

    # Final columns
    df["p_value"] = df[p_col]
    df["genome_coverage_percent_minimizers"] = df["genome_coverage_estimate"] * 100
    df = add_confidence(df, "p_value")

    df = df.sort_values(
        ["sample", "on_pathogen_panel", "p_value"],
        ascending=[True, False, True],
    ).groupby("sample", group_keys=False).head(args.top_n)

    df = df.rename(
        columns={
            "name": "organism",
            "rank_code": "rank",
            "fragments_clade": "reads_in_clade",
        }
    )

    make_html_report(df, Path(args.output))
    print(f"HTML report written to {args.output}")


if __name__ == "__main__":
    main()
