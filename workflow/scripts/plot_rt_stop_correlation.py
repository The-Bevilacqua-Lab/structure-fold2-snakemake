##########################################################################
# Pairwise hexbin correlation of log10 RT-stop counts (A and C only).
# Reads the CSV produced by StructureFold2 rtsc_correlation.py
# (columns: transcript, position, base, <rtsc_path_1>, <rtsc_path_2>, ...)
# and plots one hexbin panel per replicate pair with Pearson r.
##########################################################################

import argparse
import itertools
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(
    description="Pairwise hexbin RT-stop correlation from rtsc_correlation.py CSV"
)
parser.add_argument("--input", required=True,
                    help="CSV from rtsc_correlation.py (with -spec AC)")
parser.add_argument("--output", required=True,
                    help="Path to save the output PNG")
args = parser.parse_args()

df = pd.read_csv(args.input)

# The data columns are everything after transcript / position / base
meta_cols = {"transcript", "position", "base"}
data_cols = [c for c in df.columns if c not in meta_cols]

# Derive short display label from path: .../se/<ID>/combined_plus -> <ID>
def id_from_col(col):
    parts = col.replace("\\", "/").split("/")
    # second-to-last component is the sample ID directory
    if len(parts) >= 2:
        return parts[-2]
    return col

labels = {col: id_from_col(col) for col in data_cols}

pairs = list(itertools.combinations(data_cols, 2))
n_pairs = len(pairs)

nrows = min(3, n_pairs)
ncols = math.ceil(n_pairs / nrows)

# no-op: kept so the set_font(...) call sites below don't need touching
def set_font(text_obj):
    pass

fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows), squeeze=False)

for ax, (col_a, col_b) in zip(axes.flat, pairs):
    sub = df[[col_a, col_b]].dropna()
    x = np.log10(sub[col_a].values + 1)
    y = np.log10(sub[col_b].values + 1)

    r, _ = pearsonr(sub[col_a].values, sub[col_b].values)

    hb = ax.hexbin(x, y, gridsize=60, mincnt=1, cmap="viridis", bins="log")
    
    # Colorbar
    cb = fig.colorbar(hb, ax=ax)
    set_font(cb.set_label("log10(count)") or cb.ax.yaxis.label)
    for tick_label in cb.ax.get_yticklabels():
        set_font(tick_label)

    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, "r--", linewidth=0.8, alpha=0.6)

    # Axis labels and title
    set_font(ax.set_xlabel(f"log10(RT-stops + 1)  [{labels[col_a]}]") or ax.xaxis.label)
    set_font(ax.set_ylabel(f"log10(RT-stops + 1)  [{labels[col_b]}]") or ax.yaxis.label)
    #set_font(ax.set_title(f"{labels[col_a]} vs {labels[col_b]}  |  r = {r:.3f}") or ax.title)

    # Annotation
    txt = ax.text(0.05, 0.92, f"r = {r:.3f}", transform=ax.transAxes,
                  fontsize=8, color="black", va="top")
    set_font(txt)

    # Tick labels
    for tick_label in ax.get_xticklabels():
        set_font(tick_label)
    for tick_label in ax.get_yticklabels():
        set_font(tick_label)

for ax in axes.flat[n_pairs:]:
    ax.set_visible(False)

plt.tight_layout()
plt.savefig(args.output, dpi=300, bbox_inches="tight")
plt.close()