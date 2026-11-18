"""Pie chart of where each sample's reads go: mRNA, rRNA, lncRNA, other, unmapped
(plot_read_composition_pie rule, report.smk).

Input: read_composition_summary.tsv (pie_<slice> columns: every read entering Bowtie2 in exactly one
slice, assigned by the priority in config read_composition_priority; see read_composition.awk).
One pie per sample plus one for all samples pooled. Slices with no reads in any sample are left out.
Two figures: all reads (png, with an unmapped slice) and mapped reads only (png_mapped, unmapped left
out; antisense-only reads are mapped and count as other). The TSV has the counts, % of all reads and %
of mapped reads.
"""
import os
import sys

import numpy as np
import pandas as pd
from plotnine import (aes, coord_fixed, element_text, facet_wrap, geom_polygon, geom_text, ggplot, labs,
                      scale_fill_manual, scale_x_continuous, scale_y_continuous, theme, theme_void)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")

SLICES = list(sm.params.slices)                       # priority types, then other, unmapped
LABEL = {"mRNA": "mRNA", "rRNA": "rRNA", "lncRNA": "lncRNA", "ncRNA": "other ncRNA",
         "positive_control": "positive control", "other": "other", "unmapped": "unmapped"}
COLOR = {"mRNA": "#2a6fb0", "rRNA": "#d1495b", "lncRNA": "#b5487a", "ncRNA": "#e08a1e",
         "positive_control": "#3a9e5c", "other": "#8c8c8c", "unmapped": "#d9d9d9"}

comp = pd.read_csv(sm.input[0], sep="\t").set_index("sample")
counts = comp[[f"pie_{s}" for s in SLICES]].rename(columns=lambda c: c[4:])
counts.loc["All samples"] = counts.sum()
counts = counts.loc[:, counts.sum() > 0]
mapped = counts.drop(columns="unmapped", errors="ignore")
pct = 100 * counts.div(counts.sum(axis=1), axis=0)
pct_m = 100 * mapped.div(mapped.sum(axis=1), axis=0)
out = pd.concat({"reads": counts, "pct": pct.round(2), "pct_of_mapped": pct_m.round(2)}, axis=1)
out.columns = [f"{s}_{m}" for m, s in out.columns]
out.insert(0, "mapped_reads", mapped.sum(axis=1))
out.insert(0, "total_reads", counts.sum(axis=1))
out.to_csv(sm.output.tsv, sep="\t")


def pie(counts, what, path):
    """One pie per row of counts (wedges as polygons: plotnine has no polar coordinates), clockwise
    from 12 o'clock; % labels on slices >= 3% (outside the pie below 10%)."""
    slices = list(counts.columns)
    poly, lab = [], []
    panels = list(counts.index)
    for smp in panels:
        frac = (counts.loc[smp] / counts.loc[smp].sum()).to_numpy()
        edges = np.concatenate([[0], np.cumsum(frac)])
        for s, a0, a1 in zip(slices, edges[:-1], edges[1:]):
            if a1 - a0 <= 0:
                continue
            th = np.linspace(a0, a1, max(2, int(np.ceil((a1 - a0) * 200)) + 1)) * 2 * np.pi
            x = np.concatenate([[0], np.sin(th), [0]]); y = np.concatenate([[0], np.cos(th), [0]])
            poly.append(pd.DataFrame({"x": x, "y": y, "slice": s, "sample": smp, "grp": f"{smp}|{s}"}))
            if a1 - a0 >= 0.03:
                m = (a0 + a1) * np.pi
                r = 0.62 if a1 - a0 >= 0.10 else 1.2         # small slices: label just outside the pie
                lab.append({"x": r * np.sin(m), "y": r * np.cos(m), "sample": smp, "slice": s,
                            "label": f"{100 * (a1 - a0):.0f}%"})
    poly = pd.concat(poly, ignore_index=True)
    lab = pd.DataFrame(lab)
    title = {smp: f"{smp}\n{int(counts.loc[smp].sum()):,} {what}" for smp in panels}
    for d in (poly, lab):
        d["sample"] = pd.Categorical(d["sample"].map(title), categories=[title[p] for p in panels])
        d["slice"] = pd.Categorical(d["slice"].map(LABEL), categories=[LABEL[s] for s in slices])
    ncol = min(4, len(panels))
    nrow = int(np.ceil(len(panels) / ncol))
    p = (ggplot(poly, aes("x", "y", group="grp", fill="slice"))
         + geom_polygon(color="white", size=0.3)
         + geom_text(lab, aes("x", "y", label="label"), inherit_aes=False, size=8, color="black")
         + facet_wrap("sample", ncol=ncol)
         + scale_fill_manual(values={LABEL[s]: COLOR[s] for s in slices})
         + scale_x_continuous(limits=(-1.38, 1.38), expand=(0, 0))
         + scale_y_continuous(limits=(-1.38, 1.38), expand=(0, 0))
         + coord_fixed()
         + labs(fill="")
         + theme_void(base_size=10)
         + theme(figure_size=(3.0 * ncol + 1.6, 3.1 * nrow), legend_position="right",
                 strip_text=element_text(size=9, margin={"b": 4, "units": "pt"}),
                 panel_spacing=0.03, plot_margin=0.04))
    p.save(path, dpi=200, verbose=False)


pie(counts, "reads", sm.output.png)
pie(mapped, "mapped reads", sm.output.png_mapped)
