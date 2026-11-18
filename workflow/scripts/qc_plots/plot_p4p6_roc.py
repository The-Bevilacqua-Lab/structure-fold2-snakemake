#!/usr/bin/env python3
"""
ROC curves of p4p6 DMS reactivity against the known P4-P6 secondary structure (plotnine).

Classes come from the XRNA drawing's base pairs, over A/C residues (the only ones DMS reports on):
  positive  unpaired A/C
  negative  A/C in a Watson-Crick-Franklin pair (A-U or C-G)
A/C in any other pair (G-A, A-A, ...) are left out. High reactivity predicts "unpaired", so each curve
sweeps a reactivity cut-off and plots the true-positive rate (unpaired A/C called unpaired) against
the false-positive rate (WCF-paired A/C called unpaired). The AUC is the chance that a random unpaired
A/C is more reactive than a random WCF-paired one (ties count half).

Reactivity position p of the p4p6 transcript is XRNA nucleotide p - --offset. --trim3 leaves out the
last N transcript positions (default 50; poorly covered by RT). One curve per --inputs file; writes the plot and
a CSV of AUC and class counts per label.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from plotnine import (aes, coord_equal, element_blank, element_text, geom_abline, geom_path, ggplot, labs,
                      scale_color_manual, scale_x_continuous, scale_y_continuous, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


WCF = {"AU", "UA", "CG", "GC"}
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
           "#bcbd22", "#17becf"]


def read_xrna(path):
    """Sequence and base pairs (1-based (i, j)) from an XRNA file."""
    text = open(path).read()
    block = text.split("DataType='NucChar.XPos.YPos'>")[1].split("<")[0]
    seq = "".join(line.split()[0].upper().replace("T", "U") for line in block.strip().splitlines() if line.strip())
    pairs = []
    for a, n, b in re.findall(r"<BasePairs nucID='(\d+)' length='(\d+)' bpNucID='(\d+)'", text):
        a, n, b = int(a), int(n), int(b)
        pairs += [(a + k, b - k) for k in range(n)]
    return seq, pairs


def structure_classes(seq, pairs):
    """XRNA nucleotide -> 'wcf' / 'unpaired' / 'other' (other = in a non-WCF pair)."""
    cls = {i: "unpaired" for i in range(1, len(seq) + 1)}
    for i, j in pairs:
        kind = "wcf" if seq[i - 1] + seq[j - 1] in WCF else "other"
        cls[i] = cls[j] = kind
    return cls


def roc_points(score, positive):
    """(FPR, TPR) at every distinct score threshold, highest score first; ties move diagonally."""
    order = np.argsort(-score, kind="stable")
    s, y = score[order], positive[order]
    last = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]       # last index of each tie group
    tp, fp = np.cumsum(y)[last], np.cumsum(~y)[last]
    return np.r_[0, fp / (~y).sum()], np.r_[0, tp / y.sum()]


def auc(score, positive):
    ranks = pd.Series(score).rank(method="average").to_numpy()
    n_pos, n_neg = positive.sum(), (~positive).sum()
    return (ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True, help="reactivity.csv (transcript,position,base,reactivity)")
    ap.add_argument("--labels", nargs="+", required=True, help="One curve label per --inputs")
    ap.add_argument("--xrna", required=True, help="P4-P6 XRNA structure drawing")
    ap.add_argument("--transcript", default="p4p6")
    ap.add_argument("--offset", type=int, default=2, help="[default = 2] reactivity position minus XRNA number")
    ap.add_argument("--trim3", type=int, default=50, help="[default = 50] leave out the last N transcript positions")
    ap.add_argument("--title", default="p4p6: reactivity vs known structure")
    ap.add_argument("--plot", required=True)
    ap.add_argument("--auc-csv", required=True)
    args = ap.parse_args()
    if len(args.inputs) != len(args.labels):
        ap.error("--inputs and --labels must be the same length")

    seq, pairs = read_xrna(args.xrna)
    cls = structure_classes(seq, pairs)
    curves, rows = [], []
    for path, label in zip(args.inputs, args.labels):
        df = pd.read_csv(path)
        df.columns = ["transcript", "position", "base", "reactivity"]
        df = df[df["transcript"] == args.transcript].copy()
        if df.empty:
            sys.exit(f"{path}: no rows for transcript {args.transcript!r}")
        last = df["position"].max()
        df["nt"] = df["position"] - args.offset
        df["reactivity"] = pd.to_numeric(df["reactivity"], errors="coerce")
        df = df[(df["nt"] >= 1) & (df["nt"] <= len(seq)) & (df["position"] <= last - args.trim3)
                & df["reactivity"].notna() & df["base"].str.upper().isin(["A", "C"])]
        df["class"] = df["nt"].map(cls)
        df = df[df["class"].isin(["wcf", "unpaired"])]
        positive = (df["class"] == "unpaired").to_numpy()
        score = df["reactivity"].to_numpy()                       # high reactivity -> called unpaired
        a = auc(score, positive)
        fpr, tpr = roc_points(score, positive)
        name = f"{label} (AUC {a:.3f})"
        curves.append(pd.DataFrame({"fpr": fpr, "tpr": tpr, "curve": name}))
        rows.append({"label": label, "auc": round(a, 4), "n_unpaired": int(positive.sum()),
                     "n_wcf_paired": int((~positive).sum()), "trim3": args.trim3})
        print(f"{label}: AUC {a:.4f} ({positive.sum()} unpaired, {(~positive).sum()} WCF-paired A/C)")

    pd.DataFrame(rows).to_csv(args.auc_csv, index=False)
    data = pd.concat(curves, ignore_index=True)
    names = [c["curve"].iloc[0] for c in curves]
    data["curve"] = pd.Categorical(data["curve"], categories=names)
    sub = f"positive = unpaired A/C (n = {rows[0]['n_unpaired']}), negative = WCF-paired A/C (n = {rows[0]['n_wcf_paired']})"
    if args.trim3:
        sub += f"; last {args.trim3} nt left out"
    p = (ggplot(data, aes("fpr", "tpr", color="curve"))
         + geom_abline(slope=1, intercept=0, color="grey", linetype="dashed")
         + geom_path(size=1)
         + scale_color_manual(values=[PALETTE[i % len(PALETTE)] for i in range(len(names))], name="")
         + scale_x_continuous(limits=(0, 1), expand=(0.01, 0.01))
         + scale_y_continuous(limits=(0, 1), expand=(0.01, 0.01))
         + coord_equal()
         + labs(x="False-positive rate (WCF-paired called unpaired)", y="True-positive rate (unpaired called unpaired)",
                title=args.title, subtitle=sub)
         + theme_bw(base_size=10)
         + theme(figure_size=(6.4, 5), panel_grid_minor=element_blank(), legend_position="right",
                 legend_title=element_blank(),
                 plot_subtitle=element_text(size=8, color="#555555")))
    p.save(args.plot, dpi=300, verbose=False)


if __name__ == "__main__":
    main()
