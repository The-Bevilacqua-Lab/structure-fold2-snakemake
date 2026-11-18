#!/usr/bin/env python3
"""
Overlay the per-sample metagene profiles from metagene_profile.py on one plot
(plotnine). One line per sample; colour identifies the sample, line type the
DMS condition (solid = plus / +DMS, dashed = minus / -DMS).
"""

import argparse
import os
import sys

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_text,
    geom_hline,
    geom_line,
    ggplot,
    guides,
    guide_legend,
    labs,
    scale_color_manual,
    scale_linetype_manual,
    scale_x_continuous,
    theme,
    theme_bw,
)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


# matplotlib's tab10 / tab20, so colours match the earlier matplotlib version
TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
         "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
TAB20 = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
         "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
         "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
         "#17becf", "#9edae5"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profiles", nargs="+", required=True, help="metagene_profile.py TSVs")
parser.add_argument("--samples", nargs="+", required=True, help="Sample name for each profile, same order")
parser.add_argument("--conditions", nargs="+", required=True, help="Condition (plus/minus) for each profile")
parser.add_argument("--title", default="Metagene coverage")
parser.add_argument("--output", required=True)
args = parser.parse_args()

if not (len(args.profiles) == len(args.samples) == len(args.conditions)):
    parser.error("--profiles, --samples and --conditions must have the same length")

frames = []
labels = []
for path, sample, cond in zip(args.profiles, args.samples, args.conditions):
    label = f"{sample} ({'-' if cond == 'minus' else '+'}DMS)"
    labels.append(label)
    df = pd.read_csv(path, sep="\t")
    df["sample"] = label
    df["condition"] = "-DMS" if cond == "minus" else "+DMS"
    frames.append(df)

data = pd.concat(frames, ignore_index=True)
# Keep samples (legend and colours) in samplesheet order
data["sample"] = pd.Categorical(data["sample"], categories=labels, ordered=True)
data["condition"] = pd.Categorical(data["condition"], categories=["+DMS", "-DMS"])

palette = TAB20 if len(labels) > 10 else TAB10
colours = [palette[i % len(palette)] for i in range(len(labels))]
# Draw each sample's legend key with its own line type (dashed for -DMS)
key_linetypes = ["dashed" if c == "minus" else "solid" for c in args.conditions]

plot = (
    ggplot(data, aes("position_pct", "mean_norm_coverage", colour="sample", linetype="condition"))
    + geom_hline(yintercept=1.0, colour="grey", size=0.4, linetype="dotted")
    + geom_line(size=0.8)
    + scale_color_manual(values=colours, name="Sample")
    + scale_linetype_manual(values={"+DMS": "solid", "-DMS": "dashed"}, name="Condition", drop=False)
    + scale_x_continuous(limits=(0, 100), breaks=range(0, 101, 20), expand=(0, 0))
    + guides(
        colour=guide_legend(ncol=1, order=1, override_aes={"linetype": key_linetypes}),
        linetype=guide_legend(order=2),
    )
    + labs(
        x="Relative position along transcript (5' to 3', %)",
        y="ln(coverage + 1) / transcript mean",
        title=args.title,
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(9, 5),
        panel_grid_minor=element_blank(),
        legend_text=element_text(size=7),
        legend_title=element_text(size=8),
        plot_title=element_text(size=11),
    )
)

plot.save(args.output, dpi=200, verbose=False)
