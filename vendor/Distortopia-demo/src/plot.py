#!/usr/bin/env python3

"""
Distortopia plotting: genome-wide crossover distribution with chromosome labels.

This script reads a TSV produced by `disto infer` and generates a single PDF:
- x-axis: concatenated genome coordinate (Mb), labeled by chromosome/scaffold
- y-axis: crossover frequency (per site per read)
- vertical dashed lines mark chromosome boundaries
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import toyplot
import toyplot.pdf
from loguru import logger


def _parse_cx_from_read(read_name: str) -> float:
    """
    Extract crossover position encoded in a read name like "...cx=12345:..."
    Returns np.nan if parsing fails.
    """
    s = str(read_name)
    if "cx=" not in s:
        return np.nan
    try:
        return float(s.split("cx=")[1].split(":")[0])
    except Exception:
        return np.nan


def run_plot(tsv: Path, outdir: Path, prefix: str) -> None:
    """
    Create a genome-wide crossover distribution plot.

    Parameters
    ----------
    tsv : Path
        TSV produced by `disto infer`
    outdir : Path
        Output directory for plots (created if missing)
    prefix : str
        Output filename prefix (PDF will be <prefix>.pdf)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{prefix}.pdf"

    # --- user-tunable parameters (keep these simple) ---
    bins_per_chr = 100
    nreads_total = 1e7   # used only for scaling to "per read"; keep consistent with your pipeline
    read_length = 100_000  # mask first/last read_length bp of each chromosome

    # --- load data ---
    data = pd.read_csv(tsv, sep=r"\s+", engine="python")

    required_cols = {"scaff", "start", "end", "read"}
    missing = required_cols - set(data.columns)
    if missing:
        raise ValueError(f"TSV is missing required columns: {sorted(missing)}")

    # --- infer chromosome/scaffold lengths from TSV (max end per scaff) ---
    chr_lengths = (
        data.groupby("scaff")["end"]
        .max()
        .sort_values(ascending=False)
    )
    chroms = list(chr_lengths.index)

    # --- compute cumulative genome coordinates (bp) ---
    chr_starts_bp: dict[str, int] = {}
    running = 0
    for c in chroms:
        chr_starts_bp[str(c)] = running
        running += int(chr_lengths[c])
    genome_length_bp = running

    # --- build genome-wide x/y arrays ---
    all_x_mb = []
    all_y = []

    # for labeling
    chr_midpoints_mb = []
    chr_labels = []
    chr_boundaries_mb = []

    for chrom in chroms:
        chrom = str(chrom)
        sub = data.loc[data["scaff"] == chrom].copy()
        if sub.empty:
            continue

        # mask ends (your original intent)
        sub = sub.loc[sub["start"] > read_length]
        if sub.empty:
            continue
        sub = sub.loc[sub["end"] < (sub["end"].max() - read_length)]
        if sub.empty:
            continue

        # crossover positions within chromosome
        cxs = sub["read"].apply(_parse_cx_from_read).dropna().values
        if cxs.size == 0:
            # no detectable cx info for this chrom
            # still add labels/boundary so plot is complete
            pass
        else:
            # histogram within chromosome (bp)
            mags, edges = np.histogram(cxs, bins=bins_per_chr)

            # scale to "per site per read" like your original code
            mags = mags / nreads_total

            # convert edges to midpoints (bp)
            mids_bp = (edges[:-1] + edges[1:]) / 2.0

            # shift into genome coordinate and convert to Mb
            mids_bp_shifted = mids_bp + chr_starts_bp[chrom]
            all_x_mb.append(mids_bp_shifted / 1e6)
            all_y.append(mags)

        # boundaries + labels
        start_bp = chr_starts_bp[chrom]
        end_bp = chr_starts_bp[chrom] + int(chr_lengths[chrom])

        chr_boundaries_mb.append(start_bp / 1e6)
        chr_midpoints_mb.append(((start_bp + end_bp) / 2.0) / 1e6)

        # label: keep scaffold name as-is
        chr_labels.append(chrom)

    # add final boundary at genome end
    chr_boundaries_mb.append(genome_length_bp / 1e6)

    if len(all_x_mb) == 0:
        logger.warning("No crossover positions found after filtering; plot will be empty/flat.")
        # still create an empty plot with chromosome labels
        x = np.array([0.0, genome_length_bp / 1e6])
        y = np.array([0.0, 0.0])
    else:
        x = np.concatenate(all_x_mb)
        y = np.concatenate(all_y)

        # sort by x for a clean line
        order = np.argsort(x)
        x = x[order]
        y = y[order]

    # --- plot ---
    canvas = toyplot.Canvas(width=1100, height=380)
    axes = canvas.cartesian(
        xlabel="Chromosome",
        ylabel="Crossover frequency (per site per read)",
        margin=70,
    )

    # x-axis ticks: explicit chromosome labels at midpoints
    axes.x.ticks.locator = toyplot.locator.Explicit(
        locations=chr_midpoints_mb,
        labels=chr_labels,
    )
    axes.x.ticks.angle = -45
    axes.x.ticks.labels.offset = 10
    axes.x.label.offset = 32
    axes.x.label.style["font-size"] = 14

    axes.y.ticks.labels.offset = 10
    axes.y.label.offset = 32
    axes.y.label.style["font-size"] = 14

    # chromosome boundary lines
    for b in chr_boundaries_mb:
        axes.vlines(b, style={"stroke": "lightgray", "stroke-dasharray": "2,2"})

    # line plot (bars removed for clarity on genome-wide plot)
    line = axes.plot(x, y, opacity=0.9)
    # keep your line color tweak
    line._stroke.color = "#7b6de2"

    toyplot.pdf.render(canvas, str(out))
    logger.info(f"wrote {out}")

if __name__ == "__main__":
    main()
