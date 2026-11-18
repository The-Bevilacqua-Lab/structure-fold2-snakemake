##########################################################################
# Pairwise RT-stop correlation (A and C only), plotnine.
# Reads the CSV produced by StructureFold2 rtsc_correlation.py
# (columns: transcript, position, base, <rtsc_path_1>, <rtsc_path_2>, ...)
# and draws one 2-D density panel per replicate pair of log10(RT-stops + 1),
# with Pearson r of the raw counts. Counts are binned with numpy first, so the
# tens of millions of positions never become individual plot points. The bands
# near zero are real: counts are integers, so log10(count + 1) only takes the
# values 0, 0.30, 0.48, 0.60, ... at low counts.
##########################################################################

import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd
from plotnine import (aes, element_blank, facet_wrap, geom_abline, geom_text, geom_tile, ggplot, labs,
                      scale_fill_cmap, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


parser = argparse.ArgumentParser(description="Pairwise RT-stop correlation from rtsc_correlation.py CSV")
parser.add_argument("--input", required=True, help="CSV from rtsc_correlation.py (with -spec AC)")
parser.add_argument("--output", required=True, help="Path to save the output PNG")
parser.add_argument("--bins", type=int, default=80, help="Bins per axis for the 2-D density")
args = parser.parse_args()

header = pd.read_csv(args.input, nrows=0).columns
data_cols = [c for c in header if c not in {"transcript", "position", "base"}]
df = pd.read_csv(args.input, usecols=data_cols, dtype="float32")


def id_from_col(col):
    """.../se/<ID>/combined_plus -> <ID>"""
    parts = col.replace("\\", "/").split("/")
    return parts[-2] if len(parts) >= 2 else col


labels = {c: id_from_col(c) for c in data_cols}
logv = {c: np.log10(df[c].to_numpy(dtype=np.float64) + 1) for c in data_cols}
hi = max(float(np.nanmax(v)) for v in logv.values())
edges = np.linspace(0, hi * 1.0001, args.bins + 1)
width = edges[1] - edges[0]
centers = (edges[:-1] + edges[1:]) / 2

tiles, notes = [], []
pairs = list(itertools.combinations(data_cols, 2))
for a, b in pairs:
    ok = np.isfinite(logv[a]) & np.isfinite(logv[b])
    h, _, _ = np.histogram2d(logv[a][ok], logv[b][ok], bins=[edges, edges])
    xi, yi = np.nonzero(h)
    name = f"{labels[a]} vs {labels[b]}"
    tiles.append(pd.DataFrame({"x": centers[xi], "y": centers[yi], "positions": h[xi, yi], "pair": name}))
    r = np.corrcoef(df[a].to_numpy()[ok], df[b].to_numpy()[ok])[0, 1]   # Pearson r of the raw counts
    notes.append({"pair": name, "x": 0.04 * hi, "y": 0.95 * hi, "label": f"r = {r:.3f}"})
tiles = pd.concat(tiles, ignore_index=True)
notes = pd.DataFrame(notes)
order = [f"{labels[a]} vs {labels[b]}" for a, b in pairs]
tiles["pair"] = pd.Categorical(tiles["pair"], categories=order)
notes["pair"] = pd.Categorical(notes["pair"], categories=order)

ncol = min(3, len(pairs))
nrow = int(np.ceil(len(pairs) / ncol))
p = (ggplot(tiles, aes("x", "y", fill="positions"))
     + geom_tile(width=width, height=width)
     + geom_abline(slope=1, intercept=0, color="red", linetype="dashed", alpha=0.6)
     + geom_text(notes, aes("x", "y", label="label"), inherit_aes=False, ha="left", va="top", size=9)
     + facet_wrap("pair", ncol=ncol)
     + scale_fill_cmap("viridis", trans="log10", name="positions")
     + labs(x="log10(RT-stops + 1), first replicate", y="log10(RT-stops + 1), second replicate")
     + theme_bw(base_size=10)
     # square panels: size the figure to fit them, or plotnine pads the sides with blank space
     + theme(figure_size=(3.5 * ncol + 1.7, 3.6 * nrow + 0.6), aspect_ratio=1, panel_grid_minor=element_blank()))
p.save(args.output, dpi=300, verbose=False)
