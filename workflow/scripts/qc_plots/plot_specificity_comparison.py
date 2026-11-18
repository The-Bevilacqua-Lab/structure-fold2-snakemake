"""
RT-stop base specificity compared across several pipeline runs (plotnine).

Each --run NAME OUTPUT_DIR SAMPLESHEET names one pipeline run; its samples' specificity CSVs
(OUTPUT_DIR/specificity_<sample>.csv, from rtsc_specificity.py: base, count, specificity) are read
import sys
from the samplesheet's sample/condition (and ID, if present) columns. One panel per run, in the order
given, each with its own rows ("<ID> (±DMS)", -DMS above +DMS; runs may name their IDs differently).
Bars use the same base colours as plot_specificity.py.

Also writes --manifest: one row per sample (file, group, condition, replicate).
"""
import argparse
import os

import pandas as pd
from plotnine import (aes, coord_flip, element_blank, element_text, facet_wrap, geom_col, ggplot,
                      guide_legend, guides, labs, scale_fill_manual, scale_y_continuous, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--run", nargs=3, action="append", required=True, metavar=("NAME", "OUTPUT_DIR", "SAMPLESHEET"))
parser.add_argument("--manifest", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

rows = []
for name, out_dir, sheet_path in args.run:
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str)
    sheet = sheet.apply(lambda c: c.str.strip()).drop_duplicates(subset="sample", keep="first")  # one row per run
    for _, r in sheet.iterrows():
        rid = r["ID"] if "ID" in sheet.columns and pd.notna(r["ID"]) else r["sample"]
        rows.append({"file": os.path.join(out_dir, f"specificity_{r['sample']}.csv"), "group": name,
                     "condition": r["condition"], "replicate": rid})
manifest = pd.DataFrame(rows)
manifest.to_csv(args.manifest, sep="\t", index=False)

frames = []
for m in manifest.itertuples(index=False):
    df = pd.read_csv(m.file)
    tmp = df[[df.columns[0], df.columns[2]]].copy()
    tmp.columns = ["base", "specificity"]
    frames.append(tmp.assign(group=m.group, condition=m.condition, replicate=m.replicate))
data = pd.concat(frames, ignore_index=True)

data["label"] = data["replicate"] + data["condition"].map({"minus": " (−DMS)", "plus": " (+DMS)"})
replicates = list(dict.fromkeys(manifest.sort_values("replicate")["replicate"]))
order = [f"{r} (−DMS)" for r in replicates] + [f"{r} (+DMS)" for r in replicates]
order = [lab for lab in order if lab in set(data["label"])]
data["label"] = pd.Categorical(data["label"], categories=list(reversed(order)))  # coord_flip: first at top
data["group"] = pd.Categorical(data["group"], categories=[r[0] for r in args.run])
data["base"] = pd.Categorical(data["base"].replace({"T": "U"}), categories=["A", "C", "G", "U"])
colors = {"A": "#4C72B0", "U": "#DD8452", "G": "#55A868", "C": "#C44E52"}

n_groups = len(args.run)
n_rows = int(data.groupby("group", observed=True)["label"].nunique().max())
p = (ggplot(data, aes("label", "specificity", fill="base"))
     + geom_col(width=0.8, color="white", size=0.3)
     + coord_flip(ylim=(0, 1))
     + facet_wrap("group", nrow=1, scales="free_y")   # runs may use different ID names
     + scale_fill_manual(values=colors, name="Base")
     + scale_y_continuous(expand=(0, 0), breaks=[0, 0.25, 0.5, 0.75, 1], labels=["0", ".25", ".5", ".75", "1"])
     + guides(fill=guide_legend(reverse=True))
     + labs(x="Sample", y="Proportion of stops")
     + theme_bw(base_size=10)
     + theme(figure_size=(3.0 * n_groups + 2.2, max(3.0, 0.38 * n_rows + 1.3)),
             panel_grid_major_y=element_blank(), panel_grid_minor=element_blank(),
             strip_text=element_text(weight="bold"), panel_spacing_x=0.04))
p.save(args.output, dpi=300, verbose=False)
print(f"{len(manifest)} samples from {n_groups} runs; saved plot to {args.output}")
