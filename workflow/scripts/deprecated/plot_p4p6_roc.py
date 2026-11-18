#!/usr/bin/env python3
"""
P4-P6 ROC curves and AUC, before and after depth matching.

Score = reactivity, positive class = unpaired, at the A/C positions of the
positive control that have a 0/1 pairing label and a reactivity, leaving out
the last --trim3 nt (poor RT coverage at the 3' end). The pairing file's
index is the transcript position minus --offset.

One panel per set (--before, --after), one curve per library, AUC in the legend;
the AUC table has one row per set x library.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_PAIRING = Path(__file__).resolve().parent.parent.parent / "resources" / "p4p6" / "p4p6_pairing.txt"


def load(path, name):
    parts = [c[c["transcript"] == name] for c in pd.read_csv(path, chunksize=1_000_000)]
    return pd.concat(parts).set_index("position")["reactivity"]


def roc(score, unpaired):
    order = np.argsort(-score, kind="stable")
    y = unpaired[order]
    tpr = np.concatenate([[0], np.cumsum(y) / y.sum()])
    fpr = np.concatenate([[0], np.cumsum(~y) / (~y).sum()])
    return fpr, tpr


def auc(score, unpaired):
    ranks = pd.Series(score).rank(method="average").to_numpy()
    n_pos, n_neg = unpaired.sum(), (~unpaired).sum()
    return (ranks[unpaired].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--before", nargs="*", default=[], help="label=reactivity.csv, before depth matching")
    ap.add_argument("--after", nargs="+", required=True, help="label=reactivity.csv, depth matched")
    ap.add_argument("--pairing", default=str(DEFAULT_PAIRING),
                    help="index,seq,pairing CSV (1 = paired, 0 = unpaired); default: resources/p4p6/p4p6_pairing.txt")
    ap.add_argument("--name", default="p4p6")
    ap.add_argument("--offset", type=int, default=2)
    ap.add_argument("--trim3", type=int, default=50)
    ap.add_argument("--title", default="")
    ap.add_argument("--plot", required=True)
    ap.add_argument("--auc-csv", required=True)
    args = ap.parse_args()

    pairing = pd.read_csv(args.pairing, dtype=str).dropna()
    pairing = pairing[pairing["pairing"].isin(["0", "1"])]
    paired = {int(i) + args.offset: int(p) for i, p in pairing[["index", "pairing"]].values}

    sets = [("Before depth matching", args.before), ("After depth matching", args.after)]
    sets = [(t, s) for t, s in sets if s]
    fig, axes = plt.subplots(1, len(sets), figsize=(4.6 * len(sets), 4.4), squeeze=False)
    colors = plt.get_cmap("tab10").colors
    rows = []
    for ax, (title, specs) in zip(axes[0], sets):
        for k, spec in enumerate(specs):
            label, path = spec.split("=", 1)
            react = load(path, args.name)
            last = react.index.max()
            pos = [p for p in react.index if p in paired and p <= last - args.trim3 and react[p] == react[p]]
            score = react[pos].to_numpy(dtype=float)
            unpaired = np.array([paired[p] == 0 for p in pos])
            a = auc(score, unpaired)
            rows.append(dict(set=title, library=label, n_positions=len(pos), n_unpaired=int(unpaired.sum()), auc=a))
            fpr, tpr = roc(score, unpaired)
            ax.plot(fpr, tpr, color=colors[k % 10], lw=1.6, label=f"{label} (AUC {a:.3f})")
        ax.plot([0, 1], [0, 1], color="0.6", lw=0.8, ls="--")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("False positive rate (paired)")
        ax.set_ylabel("True positive rate (unpaired)")
        ax.set_aspect("equal")
        ax.legend(frameon=False, fontsize=7.5, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{args.name} ROC, A/C positions, last {args.trim3} nt excluded"
                 + (f"\n{args.title}" if args.title else ""), fontsize=10.5)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    pd.DataFrame(rows).to_csv(args.auc_csv, index=False, float_format="%.4f")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
