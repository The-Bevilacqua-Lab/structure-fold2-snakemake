"""
Per-sample RT-stop base specificity as horizontal stacked bars (plotnine).
Inputs are rtsc_specificity.py CSVs (base, count, specificity); -DMS samples are listed first,
then +DMS, each block in the order given.
"""
import argparse
import os
import sys

import pandas as pd
from plotnine import (aes, coord_flip, element_blank, geom_col, ggplot, guide_legend, guides, labs,
                      scale_fill_manual, scale_y_continuous, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


parser = argparse.ArgumentParser(description="Plot per-sample base specificity as horizontal stacked bar graph")
parser.add_argument("--inputs", nargs="+", required=True)
parser.add_argument("--samples", nargs="+", required=True)
parser.add_argument("--conditions", nargs="+", required=True)
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()

frames = []
for path, sample, condition in zip(args.inputs, args.samples, args.conditions):
    df = pd.read_csv(path)
    tmp = df[[df.columns[0], df.columns[2]]].copy()
    tmp.columns = ["base", "specificity"]
    tmp["sample"] = sample
    tmp["condition"] = condition
    frames.append(tmp)
data = pd.concat(frames, ignore_index=True)

# -DMS samples first, then +DMS, each in input order; coord_flip draws the first level at the bottom,
# so reverse to read top-to-bottom
order = ([s for s, c in zip(args.samples, args.conditions) if c == "minus"]
         + [s for s, c in zip(args.samples, args.conditions) if c == "plus"])
condition_of = dict(zip(args.samples, args.conditions))
label_of = {s: f"{s} ({'−' if condition_of[s] == 'minus' else '+'}DMS)" for s in order}
data["label"] = pd.Categorical(data["sample"].map(label_of), categories=[label_of[s] for s in reversed(order)])
data["base"] = data["base"].replace({"T": "U"})
# position_stack puts the first level at the end of the bar: A right-most, then C, G, U from the left
data["base"] = pd.Categorical(data["base"], categories=["A", "C", "G", "U"])
colors = {"A": "#4C72B0", "U": "#DD8452", "G": "#55A868", "C": "#C44E52"}

p = (ggplot(data, aes("label", "specificity", fill="base"))
     + geom_col(width=0.8, color="white", size=0.3)
     + coord_flip(ylim=(0, 1))   # zoom, never drop bars that sum to just over 1
     + scale_fill_manual(values=colors, name="Base")
     + scale_y_continuous(expand=(0, 0))
     + guides(fill=guide_legend(reverse=True))
     + labs(x="Sample", y="Proportion of stops")
     + theme_bw(base_size=10)
     + theme(figure_size=(5.2, 0.38 * len(order) + 1.3), panel_grid_major_y=element_blank(),
             panel_grid_minor=element_blank()))
p.save(args.output, dpi=300, verbose=False)
print(f"Saved plot to {args.output}")
