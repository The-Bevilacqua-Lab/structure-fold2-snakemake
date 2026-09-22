#!/usr/bin/env python3
"""
Overlay the per-sample metagene profiles from metagene_profile.py on one plot.
One line per sample; colour identifies the sample, line style the DMS condition
(solid = plus / +DMS, dashed = minus / -DMS).
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profiles", nargs="+", required=True, help="metagene_profile.py TSVs")
parser.add_argument("--samples", nargs="+", required=True, help="Sample name for each profile, same order")
parser.add_argument("--conditions", nargs="+", required=True, help="Condition (plus/minus) for each profile")
parser.add_argument("--title", default="Metagene coverage")
parser.add_argument("--output", required=True)
args = parser.parse_args()

if not (len(args.profiles) == len(args.samples) == len(args.conditions)):
    parser.error("--profiles, --samples and --conditions must have the same length")

cmap = plt.get_cmap("tab20" if len(args.profiles) > 10 else "tab10")
fig, ax = plt.subplots(figsize=(8, 5))
for i, (path, sample, cond) in enumerate(zip(args.profiles, args.samples, args.conditions)):
    df = pd.read_csv(path, sep="\t")
    ax.plot(df["position_pct"], df["mean_norm_coverage"], color=cmap(i % cmap.N),
            linestyle="--" if cond == "minus" else "-", linewidth=1.6,
            label=f"{sample} ({'-' if cond == 'minus' else '+'}DMS)")

ax.axhline(1.0, color="grey", linewidth=0.8, linestyle=":")
ax.set_xlim(0, 100)
ax.set_xlabel("Relative position along transcript (5' to 3', %)")
ax.set_ylabel("ln(coverage + 1) / transcript mean")
ax.set_title(args.title)
ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
fig.savefig(args.output, dpi=200)
