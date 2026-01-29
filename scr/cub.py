#!/usr/bin/env python3
"""
cub.py

Kraken2 analysis:

- Reads Kraken2 standard reports produced with --report-minimizer-data.
- Reads Kraken2 DB inspect.txt (from --reference).
- Reads Kraken2 taxonomy (names.dmp, nodes.dmp) to:
    * resolve --organism (host species) to a taxid (supports name or taxID),
    * find all descendant taxids of that host,
    * mark those taxa as is_host = True,
    * exclude host taxa from:
        - minimizer proportion denominators,
        - statistical tests (z, p, FDR).

Outputs a single TSV with per-sample / per-taxon metrics.

Dependencies:
    - Python 3.8+
    - pandas
    - numpy
    - scipy

Install (if needed):
    pip install --user pandas numpy scipy
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(description=("Statistical analysis of Kraken2 reports using minimizer data, inspect.txt and host organism."))

    parser.add_argument(
        "--std-reports",
        required=True,
        help=("Directory containing Kraken2 standard reports produced with --report-minimizer-data."),
    )

    parser.add_argument(
        "--organism",
        help=(
            "Host organism. May be a scientific name (e.g. 'Homo sapiens') "
            "or a taxID (e.g. '9606'). If provided and taxonomy files are "
            "found, host taxa will be identified and excluded from "
            "statistical tests and denominators."
        ),
    )

    parser.add_argument(
        "--reference",
        required=True,
        help=("Path to Kraken2 DB inspect.txt (output of kraken2-inspect). Used to get per-taxon minimizer counts from the DB."),
    )

    parser.add_argument(
        "--domain",
        help=("Comma-separated list of domains of interest (e.g. 'Viruses,Bacteria'). Currently not used to filter."),
    )

    parser.add_argument(
        "--metadata",
        help="Optional metadata CSV with sample-level information.",
    )
    parser.add_argument(
        "--sample-col",
        help="Column in metadata file containing sample IDs.",
    )
    parser.add_argument(
        "--columns",
        help="Comma-separated list of metadata columns to keep.",
    )
    parser.add_argument(
        "--samples-to-remove",
        help="Text file with one sample ID per line to exclude.",
    )
    parser.add_argument(
        "--prefix",
        default="cub",
        help="Prefix for output files (default: 'cub').",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory (created if needed).",
    )
    parser.add_argument(
        "--min-prop",
        type=float,
        default=0.0,
        help=("Minimum non-host minimizer proportion to keep rows in final output (default: 0.0)."),
    )
    parser.add_argument(
        "--negative-control",
        help=(
            "Sample name corresponding to a negative control. "
            "If provided, taxa present in this sample will be used "
            "as a background contamination model."
        ),
    )

    args = parser.parse_args()
    return args


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------


def read_samples_to_remove(path):
    if path is None:
        return set()
    to_remove = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                to_remove.add(line)
    return to_remove


def guess_sample_id(path: Path) -> str:
    """
    Guess a sample ID from the filename: basename without extension.
    Adapt this if your naming scheme is more complex.
    """
    return path.stem


# ---------------------------------------------------------------------------
# Reading Kraken2 reports
# ---------------------------------------------------------------------------


def read_kraken_report(report_path: Path, sample_id: str) -> pd.DataFrame:
    """
    Read a Kraken2 standard report with --report-minimizer-data.

    Expected format (tab-separated):

        1.  percentage of fragments in clade
        2.  fragments in clade
        3.  fragments assigned directly to taxon
        4.  number of minimizers in read data associated with this taxon
        5.  estimate of distinct minimizers in read data associated with this taxon
        6.  rank code
        7.  taxid
        8.  indented scientific name
    """
    cols = [
        "pct_fragments",
        "fragments_clade",
        "fragments_direct",
        "minimizers",
        "distinct_minimizers",
        "rank_code",
        "taxid",
        "name",
    ]

    try:
        df = pd.read_csv(
            report_path,
            sep="\t",
            header=None,
            names=cols,
            dtype={
                "pct_fragments": float,
                "fragments_clade": float,
                "fragments_direct": float,
                "minimizers": float,
                "distinct_minimizers": float,
                "rank_code": str,
                "taxid": str,
                "name": str,
            },
            engine="python",
        )
    except Exception as e:
        raise RuntimeError(f"Error reading report {report_path}: {e}") from e

    if df.shape[1] != 8:
        raise RuntimeError(
            f"Report {report_path} does not look like a --report-minimizer-data report "
            f"(expected 8 columns, got {df.shape[1]}). "
            f"Did you run Kraken2 with --report-minimizer-data?"
        )

    df["sample"] = sample_id
    df["name"] = df["name"].astype(str).str.strip()

    return df


def load_all_reports(std_reports_dir: str, samples_to_remove: set) -> pd.DataFrame:
    """
    Walk through std-reports directory and load all Kraken2 reports.
    Returns a concatenated DataFrame with a 'sample' column.
    """
    std_dir = Path(std_reports_dir)
    if not std_dir.is_dir():
        raise RuntimeError(f"--std-reports directory not found: {std_reports_dir}")

    all_dfs = []
    n_files = 0

    for path in sorted(std_dir.iterdir()):
        if not path.is_file():
            continue
        # Accept .txt, .tsv, .report, or no extension
        if path.suffix.lower() not in {".txt", ".tsv", ".report", ""}:
            continue

        sample_id = guess_sample_id(path)
        if sample_id in samples_to_remove:
            continue

        n_files += 1
        df = read_kraken_report(path, sample_id)
        all_dfs.append(df)

    if n_files == 0:
        raise RuntimeError(f"No suitable report files found in directory: {std_reports_dir}")
    if not all_dfs:
        raise RuntimeError("All reports were excluded (possibly due to --samples-to-remove).")

    combined = pd.concat(all_dfs, ignore_index=True)
    return combined


# ---------------------------------------------------------------------------
# Reading inspect.txt
# ---------------------------------------------------------------------------


def read_inspect(reference_path: str) -> pd.DataFrame:
    """
    Read Kraken2 DB inspect.txt (output of `kraken2-inspect`).

    Typical format:

        pct  clade_minimizers  taxon_minimizers  rank_code  taxid  name...

    We parse:
        - pct               -> db_pct
        - clade_minimizers  -> db_clade_minimizers
        - taxon_minimizers  -> db_taxon_minimizers
        - rank_code         -> db_rank_code
        - taxid             -> taxid (string)
        - name              -> db_name
    """
    rows = []
    reference_path = Path(reference_path)
    if not reference_path.is_file():
        raise RuntimeError(f"--reference file not found: {reference_path}")

    with open(reference_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(maxsplit=5)
            if len(parts) < 5:
                continue

            if len(parts) == 5:
                pct_str, clade_str, taxon_str, rank_code, taxid = parts
                name = ""
            else:
                pct_str, clade_str, taxon_str, rank_code, taxid, name = parts

            try:
                db_pct = float(pct_str)
            except ValueError:
                # header / malformed lines
                continue

            try:
                clade_min = float(clade_str)
                taxon_min = float(taxon_str)
            except ValueError:
                continue

            rows.append(
                {
                    "taxid": str(taxid),
                    "db_pct": db_pct,
                    "db_clade_minimizers": clade_min,
                    "db_taxon_minimizers": taxon_min,
                    "db_rank_code": rank_code,
                    "db_name": name.strip(),
                }
            )

    if not rows:
        raise RuntimeError(f"Could not parse any rows from inspect.txt at {reference_path} – format may be different than expected.")

    inspect_df = pd.DataFrame(rows)
    return inspect_df


# ---------------------------------------------------------------------------
# Taxonomy: names.dmp & nodes.dmp for host handling
# ---------------------------------------------------------------------------


def get_taxonomy_paths(reference_path: Path):
    """
    From inspect.txt, infer the taxonomy directory containing
    names.dmp and nodes.dmp:

        db_root = reference_path.parent
        taxonomy_dir = db_root / "taxonomy"
    """
    db_root = reference_path.parent
    taxonomy_dir = db_root / "taxonomy"
    names_path = taxonomy_dir / "names.dmp"
    nodes_path = taxonomy_dir / "nodes.dmp"
    return names_path, nodes_path


def find_host_taxid(organism_name: str, names_path: Path) -> str:
    """
    Find the taxid corresponding to the given host organism scientific name
    using names.dmp.

    names.dmp lines are like:
        tax_id | name_txt | unique_name | name_class |

    We look for a line where:
        - name_class == 'scientific name'
        - name_txt (case-insensitive) == organism_name (case-insensitive)
    """
    if not names_path.is_file():
        raise RuntimeError(f"names.dmp not found at {names_path}")

    target = organism_name.strip().lower()
    candidate_taxid = None

    with open(names_path, "r") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            tax_id, name_txt, unique_name, name_class = parts[:4]
            if name_class != "scientific name":
                continue
            if name_txt.lower() == target:
                candidate_taxid = tax_id
                break

    if candidate_taxid is None:
        raise RuntimeError(f"Could not find scientific name '{organism_name}' in names.dmp at {names_path}")

    return str(candidate_taxid)


def resolve_organism_to_taxid(organism: str, names_path: Path) -> str:
    """
    Resolve --organism to a taxid.

    - If organism looks like an integer (all digits), treat it directly as taxID.
    - Otherwise, resolve it in names.dmp as a scientific name.
    """
    organism = organism.strip()

    # Case 1: taxID provided directly
    if organism.isdigit():
        return organism

    # Case 2: scientific name → taxID via names.dmp
    return find_host_taxid(organism, names_path)


def load_parent_children(nodes_path: Path):
    """
    Parse nodes.dmp and build a parent -> children mapping.

    nodes.dmp lines look like:
        tax_id | parent_tax_id | rank | ...
    """
    if not nodes_path.is_file():
        raise RuntimeError(f"nodes.dmp not found at {nodes_path}")

    parent_to_children = {}

    with open(nodes_path, "r") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            tax_id, parent_tax_id, rank = parts[:3]
            parent_tax_id = str(parent_tax_id)
            tax_id = str(tax_id)
            if parent_tax_id not in parent_to_children:
                parent_to_children[parent_tax_id] = []
            parent_to_children[parent_tax_id].append(tax_id)

    return parent_to_children


def get_descendants(root_taxid: str, parent_to_children: dict) -> set:
    """
    Given a root_taxid and a parent->children mapping, return the full set
    of descendant taxids, including the root itself.
    """
    root_taxid = str(root_taxid)
    descendants = set([root_taxid])
    stack = [root_taxid]

    while stack:
        current = stack.pop()
        children = parent_to_children.get(current, [])
        for child in children:
            if child not in descendants:
                descendants.add(child)
                stack.append(child)

    return descendants


def mark_host_taxa(df: pd.DataFrame, organism: str, reference_path: Path) -> pd.DataFrame:
    """
    Use names.dmp and nodes.dmp to:
        - resolve organism (name or taxID) -> host_taxid
        - find all descendant taxids (host clade)
        - add a boolean 'is_host' column to df

    If taxonomy files are missing or organism cannot be resolved,
    a warning is printed and no host marking is applied.
    """
    if organism is None:
        df["is_host"] = False
        return df

    names_path, nodes_path = get_taxonomy_paths(reference_path)

    try:
        host_taxid = resolve_organism_to_taxid(organism, names_path)
        parent_to_children = load_parent_children(nodes_path)
        host_taxids = get_descendants(host_taxid, parent_to_children)
        host_taxids = {str(t) for t in host_taxids}
    except Exception as e:
        print(
            f"[warning] Failed to resolve host organism '{organism}': {e}",
            file=sys.stderr,
        )
        df["is_host"] = False
        return df

    df["is_host"] = df["taxid"].astype(str).isin(host_taxids)
    return df


# ---------------------------------------------------------------------------
# Calculations
# ---------------------------------------------------------------------------


def compute_minimizer_proportions(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each sample, compute total distinct minimizers *excluding host taxa*
    (if 'is_host' exists and is True), then the proportion for each taxon:

        minimizer_proportion_sample =
            distinct_minimizers_taxon / sum(distinct_minimizers_over_non_host_taxa_in_sample)

    If all minimizers in a sample are host, the denominator is 0 and
    proportions are set to 0.
    """
    if "is_host" in df.columns:
        non_host = df[~df["is_host"]].copy()
    else:
        non_host = df.copy()
        non_host["is_host"] = False

    totals = non_host.groupby("sample")["distinct_minimizers"].sum().rename("total_distinct_minimizers_non_host")

    df = df.merge(totals, on="sample", how="left")
    df["total_distinct_minimizers_non_host"] = df["total_distinct_minimizers_non_host"].fillna(0.0)

    denom = df["total_distinct_minimizers_non_host"].replace(0, np.nan)
    df["minimizer_proportion_sample"] = df["distinct_minimizers"] / denom
    df["minimizer_proportion_sample"] = df["minimizer_proportion_sample"].fillna(0.0)

    return df


def apply_negative_control(df: pd.DataFrame, neg_sample: str) -> pd.DataFrame:
    """
    Use a negative control sample to estimate background contamination.

    Adds:
        - neg_control_prop
        - signal_over_neg = minimizer_proportion_sample - neg_control_prop
    """
    if neg_sample not in set(df["sample"]):
        raise RuntimeError(f"Negative control sample '{neg_sample}' not found in samples.")

    # Extract negative control proportions
    neg_df = df[df["sample"] == neg_sample][["taxid", "minimizer_proportion_sample"]].rename(
        columns={"minimizer_proportion_sample": "neg_control_prop"}
    )

    # Merge onto full table
    df = df.merge(neg_df, on="taxid", how="left")

    # Taxa absent from neg control → background = 0
    df["neg_control_prop"] = df["neg_control_prop"].fillna(0.0)

    # Background-corrected signal
    df["signal_over_neg"] = df["minimizer_proportion_sample"] - df["neg_control_prop"]

    # Do not allow negative signal (pure contaminants)
    df.loc[df["signal_over_neg"] < 0, "signal_over_neg"] = 0.0

    return df


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR correction.
    Returns array of adjusted p-values (same shape as input).
    """
    p = np.asarray(p_values, dtype=float)
    n = p.shape[0]
    order = np.argsort(p)
    ranks = np.arange(1, n + 1, dtype=float)

    p_sorted = p[order]
    p_adj_sorted = p_sorted * n / ranks
    # Enforce monotonicity
    p_adj_sorted = np.minimum.accumulate(p_adj_sorted[::-1])[::-1]
    p_adj_sorted = np.clip(p_adj_sorted, 0.0, 1.0)

    p_adj = np.empty_like(p_adj_sorted)
    p_adj[order] = p_adj_sorted
    return p_adj


def compute_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute, for each non-host sample–taxon:
      - z-score of minimizer_proportion_sample relative to that taxon's
        distribution across non-host samples
      - one-sided p-value (higher proportion -> smaller p)
      - Benjamini–Hochberg FDR across all non-host tests

    Host rows get z_score, p_value, fdr = NaN.
    """
    if "is_host" in df.columns:
        work_df = df[~df["is_host"]].copy()
    else:
        work_df = df.copy()
        work_df["is_host"] = False

    if work_df.empty:
        # No non-host taxa: just add NaNs
        df["taxon_mean_prop"] = np.nan
        df["taxon_sd_prop"] = np.nan
        df["n_samples_with_taxon"] = 0
        df["z_score"] = np.nan
        df["p_value"] = np.nan
        df["fdr"] = np.nan
        return df

    group_cols = ["taxid", "name"]

    stats = (
        work_df.groupby(group_cols)["minimizer_proportion_sample"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "taxon_mean_prop",
                "std": "taxon_sd_prop",
                "count": "n_samples_with_taxon",
            }
        )
    )

    work_df = work_df.merge(stats, on=group_cols, how="left")

    sd = work_df["taxon_sd_prop"].replace(0, np.nan)
    z = (work_df["minimizer_proportion_sample"] - work_df["taxon_mean_prop"]) / sd
    z = z.fillna(0.0)
    work_df["z_score"] = z

    work_df["p_value"] = 1.0 - norm.cdf(work_df["z_score"])
    work_df["fdr"] = benjamini_hochberg(work_df["p_value"].values)

    # Merge these statistics back onto the full df
    df = df.merge(
        work_df[
            group_cols
            + [
                "sample",
                "taxon_mean_prop",
                "taxon_sd_prop",
                "n_samples_with_taxon",
                "z_score",
                "p_value",
                "fdr",
            ]
        ],
        on=group_cols + ["sample"],
        how="left",
    )

    return df


def compute_statistics_neg_control(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute z-scores, p-values and FDR using background-corrected signal
    (signal_over_neg), excluding host taxa.
    """
    work_df = df[~df["is_host"]].copy()

    if work_df.empty:
        df["z_score_neg"] = np.nan
        df["p_value_neg"] = np.nan
        df["fdr_neg"] = np.nan
        return df

    group_cols = ["taxid", "name"]

    stats = (
        work_df.groupby(group_cols)["signal_over_neg"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "signal_mean_neg",
                "std": "signal_sd_neg",
                "count": "n_samples_with_taxon_neg",
            }
        )
    )

    work_df = work_df.merge(stats, on=group_cols, how="left")

    sd = work_df["signal_sd_neg"].replace(0, np.nan)
    z = (work_df["signal_over_neg"] - work_df["signal_mean_neg"]) / sd
    z = z.fillna(0.0)

    work_df["z_score_neg"] = z
    work_df["p_value_neg"] = 1.0 - norm.cdf(z)
    work_df["fdr_neg"] = benjamini_hochberg(work_df["p_value_neg"].values)

    df = df.merge(
        work_df[group_cols + ["sample", "z_score_neg", "p_value_neg", "fdr_neg"]],
        on=group_cols + ["sample"],
        how="left",
    )

    return df


def attach_metadata(df: pd.DataFrame, args) -> pd.DataFrame:
    """
    Optionally attach sample-level metadata.
    """
    if args.metadata is None:
        return df

    if args.sample_col is None:
        raise RuntimeError("If --metadata is provided, you must also provide --sample-col.")

    meta = pd.read_csv(args.metadata)

    cols_to_keep = None
    if args.columns:
        cols_to_keep = [col.strip() for col in args.columns.split(",") if col.strip()]
        if args.sample_col not in cols_to_keep:
            cols_to_keep.append(args.sample_col)

        missing = [c for c in cols_to_keep if c not in meta.columns]
        if missing:
            raise RuntimeError(f"Columns specified in --columns not found in metadata: {missing}")
        meta = meta[cols_to_keep]

    if args.sample_col not in meta.columns:
        raise RuntimeError(f"--sample-col='{args.sample_col}' not found in metadata columns.")

    meta = meta.rename(columns={args.sample_col: "sample"})
    df = df.merge(meta, on="sample", how="left")

    return df


def filter_to_best_kraken_taxa(df: pd.DataFrame, virus_direct_ratio: float = 0.8) -> pd.DataFrame:
    """
    Select 'best' Kraken taxa for downstream analysis, in a way that
    mimics what you care about clinically and is similar in spirit to
    Pavian's clade aggregation.

    - Non-viruses:
        * Keep only species-level rows (rank_code == 'S').
        * Use fragments_clade as the species-level abundance, which
          already includes S1/S2/descendants.

    - Viruses:
        * Consider species-like ranks (S, S1, S2).
        * Keep only nodes where direct reads dominate:
              fragments_direct / fragments_clade >= virus_direct_ratio
          This typically keeps the deepest node with real read support
          (e.g. 'Human adenovirus 5' at rank S2) and drops its ancestors.
    """
    if "rank_code" not in df.columns:
        return df

    # Mark viruses (simple heuristic; you can refine later with taxonomy)
    is_virus = df["name"].str.contains("virus", case=False, na=False)

    # ----------------- Non-viruses -----------------
    non_virus = df[~is_virus].copy()
    # keep species only
    non_virus = non_virus[non_virus["rank_code"] == "S"].copy()

    # ----------------- Viruses -----------------
    viruses = df[is_virus].copy()
    # keep only species-like ranks
    viruses = viruses[viruses["rank_code"].isin(["S", "S1", "S2"])].copy()

    if {"fragments_clade", "fragments_direct"}.issubset(viruses.columns):
        valid = viruses["fragments_clade"] > 0
        ratio = viruses.loc[valid, "fragments_direct"] / viruses.loc[valid, "fragments_clade"]

        keep_mask = pd.Series(False, index=viruses.index)
        keep_mask.loc[valid] = ratio >= virus_direct_ratio

        viruses = viruses[keep_mask].copy()

    # Combine back
    out = pd.concat([non_virus, viruses], ignore_index=True)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    samples_to_remove = read_samples_to_remove(args.samples_to_remove)

    # 1) Load all Kraken reports
    print("Loading Kraken2 reports...", file=sys.stderr)
    df = load_all_reports(args.std_reports, samples_to_remove)

    # 2) Read inspect.txt from the Kraken2 DB
    reference_path = Path(args.reference)
    print(f"Reading inspect.txt from {reference_path} ...", file=sys.stderr)
    inspect_df = read_inspect(str(reference_path))

    # 3) Mark host taxa if --organism was provided and taxonomy is available
    print("Resolving host organism (if provided)...", file=sys.stderr)
    df = mark_host_taxa(df, args.organism, reference_path)

    # 4) Drop strain-level / subspecies taxa BEFORE any stats or denominators
    print(
        "Filtering to species and higher (dropping strains/subspecies)...",
        file=sys.stderr,
    )
    df = filter_to_best_kraken_taxa(df)

    # 5) Compute sample-level minimizer proportions (host-excluded denominator)
    print(
        "Computing minimizer proportions per sample (host excluded from denominator)...",
        file=sys.stderr,
    )
    df = compute_minimizer_proportions(df)

    # 5b) Apply negative control, if provided
    if args.negative_control:
        print(
            f"Applying negative control background from sample '{args.negative_control}'...",
            file=sys.stderr,
        )
        df = apply_negative_control(df, args.negative_control)

    # 6) Merge DB minimizer information (inspect.txt) by taxid
    print("Merging DB minimizer data from inspect.txt...", file=sys.stderr)
    df = df.merge(
        inspect_df[
            [
                "taxid",
                "db_pct",
                "db_clade_minimizers",
                "db_taxon_minimizers",
                "db_rank_code",
                "db_name",
            ]
        ],
        on="taxid",
        how="left",
    )

    # 7) Compute coverage-like estimate using DB minimizers
    print("Computing coverage-like estimates...", file=sys.stderr)
    db_min = df["db_taxon_minimizers"].replace(0, np.nan)
    df["genome_coverage_estimate"] = df["distinct_minimizers"] / db_min
    df["genome_coverage_estimate"] = df["genome_coverage_estimate"]

    # 8) Filter by min_proportion if requested (on non-host denominator-based proportion)
    if args.min_prop > 0.0:
        df = df[df["minimizer_proportion_sample"] >= args.min_prop].copy()

    # 9) Compute z-scores, p-values, FDR on non-host taxa only
    print("Computing z-scores, p-values and FDR (non-host taxa only)...", file=sys.stderr)
    df = compute_statistics(df)

    # 9b) If negative control is present, compute NC-corrected stats too
    if args.negative_control:
        print(
            f"Computing statistics with negative-control correction (negative control = '{args.negative_control}')...",
            file=sys.stderr,
        )
        df = compute_statistics_neg_control(df)

    # 10) Attach metadata if provided
    if args.metadata is not None:
        print("Merging sample metadata...", file=sys.stderr)
        df = attach_metadata(df, args)

    # 11) Sort results
    df = df.sort_values(
        ["sample", "is_host", "fdr", "p_value", "minimizer_proportion_sample"],
        ascending=[True, True, True, True, False],
    )

    # 12) Write results
    out_path = outdir / f"{args.prefix}.tsv"
    print(f"Writing results to {out_path}", file=sys.stderr)
    df.to_csv(out_path, sep="\t", index=False)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
