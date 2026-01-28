import marimo

__generated_with = "0.19.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import polars as pl

    return (pl,)


@app.cell
def _():
    KRAKEN = "/mnt/san/microbio/apps/bonobo/intercalibration/cub.tsv"
    BAM = "/mnt/san/microbio/apps/bonobo/intercalibration/bam.tsv"
    return BAM, KRAKEN


@app.cell
def _(KRAKEN, pl):
    kraken = pl.read_csv(KRAKEN, separator="\t")
    kraken_sel = (
        kraken.filter(pl.col("fragments_direct") != 0)
        .select(
            "sample",
            "taxid",
            "name",
            "fragments_direct",
            "genome_coverage_estimate",
            "z_score",
            "p_value",
            "n_samples_with_taxon",
        )
        .with_columns(pl.col("taxid").cast(pl.Int64, strict=False))
        .rename(
            {
                "z_score": "kraken_z",
                "p_value": "kraken_p",
            }
        )
        # .filter(pl.col("sample") != "MG-Run35-Ech-12")
    )
    return (kraken_sel,)


@app.cell
def _(BAM, pl):
    bam = pl.read_csv(BAM, separator="\t")
    bam_sel = (
        bam.select("sample", "ncbi_taxid", "gtdb_taxonomy", "read_count", "z_score", "p_value")
        .with_columns(pl.col("ncbi_taxid").cast(pl.Int64, strict=False).alias("taxid"))
        .rename(
            {
                "z_score": "bam_z",
                "p_value": "bam_p",
            }
        )
        .drop("ncbi_taxid")
    )
    return (bam_sel,)


@app.cell
def _(bam_sel, kraken_sel, pl):
    merged = (
        bam_sel.join(kraken_sel, on=["sample", "taxid"], how="full")
        .with_columns(
            # Ensure keys are always present
            pl.coalesce([pl.col("sample"), pl.col("sample_right")]).alias("sample"),
            pl.coalesce([pl.col("taxid"), pl.col("taxid_right")]).alias("taxid"),
        )
        .select(
            "sample",
            "taxid",
            "gtdb_taxonomy",
            "name",
            "bam_z",
            "bam_p",
            "read_count",
            "fragments_direct",
            "kraken_z",
            "kraken_p",
            "n_samples_with_taxon",
            "genome_coverage_estimate",
        )
    )
    merged = (
        merged.with_columns(
            pl.col("bam_z").round(2),
            pl.col("bam_p").round(2),
            pl.col("kraken_p").round(2),
            pl.col("kraken_z").round(2),
            pl.col("genome_coverage_estimate").round(6),
        )
        .rename(
            {
                "read_count": "bam_reads",
                "fragments_direct": "kraken_reads",
            }
        )
        .filter(pl.col("name") != "Homo sapiens")
        .sort(by=["sample", "bam_reads", "kraken_reads"], descending=[False, True, True])
    )
    return (merged,)


@app.cell
def _(merged, pl):
    merged.filter(pl.col("sample").str.contains("5_"))
    return


@app.cell
def _(merged, pl):
    N = 30  # max hits per sample – adjust as you like

    best_hits = (
        merged
        # Make sure null reads become 0 for logic below
        .with_columns(
            pl.col("bam_reads").fill_null(0).alias("bam_reads"),
            pl.col("kraken_reads").fill_null(0).alias("kraken_reads"),
        )
        # Flags: present in BAM? present in Kraken?
        .with_columns(
            (pl.col("bam_reads") > 0).alias("has_bam"),
            (pl.col("kraken_reads") > 0).alias("has_kraken"),
        )
        # Evidence category + combined score
        .with_columns(
            # 2 = both, 1 = kraken-only (virus/fungus), 0 = bam-only
            pl.when(pl.col("has_bam") & pl.col("has_kraken"))
            .then(2)
            .when(~pl.col("has_bam") & pl.col("has_kraken"))
            .then(1)
            .otherwise(0)
            .alias("evidence_cat"),
            # combined strength (you can change this if you like)
            pl.when(pl.col("has_bam") & pl.col("has_kraken"))
            .then(pl.col("bam_z").fill_null(0) + pl.col("kraken_z").fill_null(0))
            .when(pl.col("has_bam"))
            .then(pl.col("bam_z").fill_null(0))
            .otherwise(pl.col("kraken_z").fill_null(0))
            .alias("combined_score"),
        )
        # Keep only reasonably strong hits:
        #   - both sets: strong in BAM & Kraken
        #   - kraken-only: strong in Kraken (for virus/fungi)
        .filter(
            # strong BAM + Kraken
            (
                pl.col("has_bam")
                & pl.col("has_kraken")
                & (pl.col("bam_z") >= 1)
                & (pl.col("kraken_z") >= 1)
                & (pl.col("bam_p") <= 0.1)
                & (pl.col("kraken_p") <= 0.1)
            )
            |
            # Kraken-only (e.g. viruses/fungi): strong Kraken signal
            (~pl.col("has_bam") & pl.col("has_kraken") & (pl.col("kraken_z") >= 1) & (pl.col("kraken_p") <= 0.1))
        )
        # Rank inside each sample:
        #   1. evidence_cat (both > kraken-only > bam-only)
        #   2. combined_score (stronger first)
        .sort(
            by=["sample", "evidence_cat", "combined_score"],
            descending=[False, True, True],
        )
        # Take top N rows per sample
        .group_by("sample", maintain_order=True)
        .head(N)
        # Drop helper columns
        .drop(["has_bam", "has_kraken", "evidence_cat", "combined_score"])
        .unique()
    )
    return (best_hits,)


@app.cell
def _(best_hits, pl):
    best_hits.with_columns(pl.col("sample").str.extract(r"(\d+)_").cast(pl.Int32).alias("id")).sort(by=["id", "gtdb_taxonomy"])
    return


if __name__ == "__main__":
    app.run()
