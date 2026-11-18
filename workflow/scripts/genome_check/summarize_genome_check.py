#!/usr/bin/env python3
"""
Summarise where the reads that did not align to the transcriptome land on the genome, over all
samples (inputs from classify_genome_hits.py and HISAT2's --summary-file).

Writes:
  --summary  one row per sample: total reads (Bowtie2 log), reads unaligned to the transcriptome,
             of those aligned to the genome (unique = one genome locus, NH = 1), and per category
             reads and % of the transcriptome-unaligned reads
  --loci     the --top loci with the most reads over all samples: read 5' ends on one chromosome
             strand no more than --gap nt apart; per-sample reads, unique reads, and the locus's
             main category / gene (the one with most reads)
  --plot     per sample, the transcriptome-unaligned reads by category (% stacked bars)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from plotnine import (aes, coord_flip, element_blank, geom_col, ggplot, guide_legend, guides, labs,
                      scale_fill_manual, scale_y_continuous, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


FIXED_FIRST = ["transcriptome exon"]
FIXED_LAST = ["intron", "antisense", "intergenic"]
NOT_ALIGNED = "not aligned to genome"
MAX_EXON_CATEGORIES = 7   # beyond this, smaller "<biotype> exon" categories become "other annotated exon"


def hisat2_counts(path):
    txt = open(path).read()
    total = int(re.search(r"(\d+) reads; of these", txt).group(1))
    zero = re.search(r"(\d+) \([\d.]+%\) aligned 0 times", txt)
    one = re.search(r"(\d+) \([\d.]+%\) aligned exactly 1 time", txt)
    many = re.search(r"(\d+) \([\d.]+%\) aligned >1 times", txt)
    get = lambda m: int(m.group(1)) if m else 0
    return total, get(zero), get(one), get(many)


def bowtie2_total(path):
    m = re.search(r"(\d+) reads; of these", open(path).read())
    return int(m.group(1)) if m else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--positions", nargs="+", required=True)
    ap.add_argument("--hisat2-summaries", nargs="+", required=True)
    ap.add_argument("--bowtie2-logs", nargs="+", required=True)
    ap.add_argument("--gap", type=int, default=50)
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--plot", required=True)
    args = ap.parse_args()
    S = args.samples

    counts = pd.concat([pd.read_csv(f, sep="\t") for f in args.counts], ignore_index=True)
    # category order: transcriptome exon, other exon categories by reads, intron/antisense/intergenic, organelles
    tot = counts.groupby("category")["reads"].sum().sort_values(ascending=False)
    exon_cats = [c for c in tot.index if c.endswith(" exon") and c not in FIXED_FIRST]
    lumped = exon_cats[MAX_EXON_CATEGORIES:]
    organelles = [c for c in tot.index if c.startswith("organelle")]
    order = ([c for c in FIXED_FIRST if c in tot.index] + exon_cats[:MAX_EXON_CATEGORIES]
             + (["other annotated exon"] if lumped else []) + [c for c in FIXED_LAST if c in tot.index] + organelles)

    rows = []
    for s, hs, bl in zip(S, args.hisat2_summaries, args.bowtie2_logs):
        total, zero, one, many = hisat2_counts(hs)
        c = counts[counts["sample"] == s].set_index("category")
        r = {"sample": s, "total_reads": bowtie2_total(bl), "transcriptome_unaligned": total,
             "genome_aligned": one + many, "genome_aligned_unique": one, "genome_aligned_multi": many}
        r["pct_of_total_unaligned"] = 100 * total / r["total_reads"] if r["total_reads"] else np.nan
        r["pct_unaligned_on_genome"] = 100 * (one + many) / total if total else np.nan
        for cat in order:
            members = lumped if cat == "other annotated exon" else [cat]
            n = int(c.reindex(members)["reads"].fillna(0).sum()) if len(c) else 0
            r[f"{cat}_reads"] = n
            r[f"{cat}_pct"] = 100 * n / total if total else np.nan
        r[f"{NOT_ALIGNED}_pct"] = 100 * zero / total if total else np.nan
        rows.append(r)
    summary = pd.DataFrame(rows)
    summary.round(3).to_csv(args.summary, sep="\t", index=False)

    # ---------------------------------------------------------------- top loci
    pos = []
    for s, f in zip(S, args.positions):
        p = pd.read_csv(f, sep="\t", dtype={"chrom": str, "gene_id": str, "gene_name": str, "biotype": str},
                        keep_default_na=False)
        pos.append(p.assign(sample=s))
    pos = pd.concat(pos, ignore_index=True)
    if len(pos):
        pos = pos.sort_values(["chrom", "strand", "pos5"], kind="stable").reset_index(drop=True)
        new = ((pos["chrom"] != pos["chrom"].shift()) | (pos["strand"] != pos["strand"].shift())
               | (pos["pos5"] - pos["pos5"].shift() > args.gap))
        pos["locus"] = new.cumsum()
        loc = pos.groupby("locus").agg(chrom=("chrom", "first"), strand=("strand", "first"), start=("pos5", "min"),
                                       end=("pos5", "max"), reads=("reads", "sum"), unique_reads=("unique_reads", "sum"))
        main = (pos.groupby(["locus", "category", "gene_id", "gene_name", "biotype"])["reads"].sum().reset_index()
                .sort_values(["locus", "reads"], ascending=[True, False]).drop_duplicates("locus").set_index("locus"))
        loc = loc.join(main[["category", "gene_id", "gene_name", "biotype"]])
        top = loc.sort_values("reads", ascending=False).head(args.top)
        per = (pos[pos["locus"].isin(top.index)].pivot_table(index="locus", columns="sample", values="reads",
                                                               aggfunc="sum", fill_value=0)
               .reindex(columns=S, fill_value=0).add_prefix("reads_"))
        top = top.join(per)
        top.insert(top.columns.get_loc("unique_reads") + 1, "pct_unique", (100 * top["unique_reads"] / top["reads"]).round(1))
        top.insert(0, "rank", np.arange(1, len(top) + 1))
    else:
        top = pd.DataFrame(columns=["rank", "chrom", "strand", "start", "end", "reads", "unique_reads", "pct_unique",
                                    "category", "gene_id", "gene_name", "biotype"])
    top.to_csv(args.loci, sep="\t", index=False)

    # ---------------------------------------------------------------- plot
    cats = order + [NOT_ALIGNED]
    long = summary.melt(id_vars="sample", value_vars=[f"{c}_pct" for c in cats], var_name="category", value_name="pct")
    long["category"] = long["category"].str[:-4]
    long["category"] = pd.Categorical(long["category"], categories=list(reversed(cats)))
    long["sample"] = pd.Categorical(long["sample"], categories=list(reversed(S)))
    base = ["#4c72b0", "#e3a21a", "#8c564b", "#9b59b6", "#17becf", "#bcbd22", "#e377c2", "#7f7f7f", "#1a3e7a"]
    exonish = [c for c in order if c.endswith(" exon") and c != "transcriptome exon"]
    pal = {"transcriptome exon": "#d42a1e", "intron": "#f2b5a0", "antisense": "#a1c9f4", "intergenic": "#555555",
           NOT_ALIGNED: "#dddddd"}
    pal.update({c: base[i % len(base)] for i, c in enumerate(exonish)})
    org_cols = ["#2ca02c", "#98df8a", "#1f5f1f"]
    pal.update({c: org_cols[i % len(org_cols)] for i, c in enumerate(organelles)})
    p = (ggplot(long, aes("sample", "pct", fill="category"))
         + geom_col(width=0.75)
         + coord_flip(ylim=(0, 100))   # zoom, never drop bars
         + scale_fill_manual(values=pal)
         + scale_y_continuous(expand=(0, 0))
         + guides(fill=guide_legend(ncol=4, reverse=True))   # legend in category order
         + labs(x="", y="% of reads not aligned to the transcriptome", fill="",
                title="Where the transcriptome-unaligned reads align on the genome")
         + theme_bw(base_size=10)
         + theme(figure_size=(10, 0.32 * len(S) + 1.6 + 0.25 * -(-len(cats) // 4)), panel_grid_major_y=element_blank(),
                 panel_grid_minor=element_blank(), legend_position="bottom"))
    p.save(args.plot, dpi=300, verbose=False)
    print(summary[["sample", "transcriptome_unaligned", "genome_aligned", "pct_unaligned_on_genome"]].to_string(index=False))


if __name__ == "__main__":
    main()
