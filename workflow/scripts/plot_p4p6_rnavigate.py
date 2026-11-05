#!/usr/bin/env python3
"""
Plot p4p6 positive-control DMS reactivity onto the Tetrahymena P4-P6 domain
secondary structure with RNAvigate (https://rnavigate.readthedocs.io).

The structure and its drawing coordinates come from an XRNA file
(resources/p4p6/p4p6.xrna). Reactivity position p of the p4p6 transcript maps
to XRNA nucleotide p - --offset (the p4p6 construct carries 2 extra 5'
nucleotides that the XRNA drawing doesn't). One panel per --inputs file,
labelled by --labels. Each residue is a circle filled by its reactivity bin;
residues without a reactivity (G/U, NA) are grey.
"""

import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import rnavigate as rnav
from matplotlib.lines import Line2D

BIN_COLORS = ["black", "orange", "red"]


def xrna_sequence(path):
    """Nucleotide letters of the XRNA NucListData block, in order."""
    seq, inside = [], False
    with open(path) as f:
        for line in f:
            if line.startswith("<NucListData"):
                inside = True
            elif line.startswith("</NucListData"):
                break
            elif inside and line.strip():
                seq.append(line.split()[0].upper().replace("T", "U"))
    return "".join(seq)


def profile_for(csv, transcript, sequence, offset, thresholds):
    df = pd.read_csv(csv)
    df.columns = ["transcript", "position", "base", "reactivity"]
    df = df[df["transcript"] == transcript]
    if df.empty:
        sys.exit(f"{csv}: no rows for transcript {transcript!r}")
    df = df.assign(Nucleotide=df["position"] - offset)
    df = df[(df["Nucleotide"] >= 1) & (df["Nucleotide"] <= len(sequence))].set_index("Nucleotide")

    xrna_bases = pd.Series(list(sequence), index=range(1, len(sequence) + 1))
    shared = xrna_bases.index.intersection(df.index)
    mismatched = [i for i in shared if df.at[i, "base"].upper().replace("T", "U") != xrna_bases[i]]
    if len(mismatched) > 2:
        print(f"WARNING: {csv}: {len(mismatched)} bases differ from the XRNA sequence "
              f"(XRNA positions {mismatched[:10]}) -- check --offset")

    table = pd.DataFrame({
        "Nucleotide": xrna_bases.index,
        "Sequence": xrna_bases.values,
        "Profile": pd.to_numeric(df["reactivity"], errors="coerce").reindex(xrna_bases.index).to_numpy(),
    })
    return rnav.data.Profile(
        input_data=table,
        metric="DMS",
        metric_defaults={
            "DMS": {
                "metric_column": "Profile",
                "error_column": None,
                "cmap": BIN_COLORS,
                "normalization": "bins",
                "values": thresholds,
                "title": "DMS reactivity",
                "extend": "max",
            }
        },
    ), int(np.isfinite(table["Profile"]).sum())


WATSON_CRICK = {"AU", "UA", "GC", "CG"}


def draw_basepairs(ax, structure):
    """Watson-Crick pairs as a short dash, every other pair (G-U, G-A, ...) as
    dots -- drawn here because RNAvigate's bp_style can't split them this way."""
    x, y, seq = structure.xcoordinates, structure.ycoordinates, structure.sequence.upper()
    for i, j in structure.get_pairs():
        (x0, x1), (y0, y1) = x[[i - 1, j - 1]], y[[i - 1, j - 1]]
        # Trim each end by the same 1/3 unit RNAvigate uses, clearing the circles.
        length = np.hypot(x1 - x0, y1 - y0)
        ux, uy = (x1 - x0) / length / 3, (y1 - y0) / length / 3
        xs, ys = [x0 + ux, x1 - ux], [y0 + uy, y1 - uy]
        if seq[i - 1] + seq[j - 1] in WATSON_CRICK:
            ax.plot(xs, ys, color="grey", linewidth=1.2, solid_capstyle="butt", zorder=0)
        else:
            ax.scatter(np.linspace(*xs, 5), np.linspace(*ys, 5), c="grey", marker=".", s=4, zorder=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inputs", nargs="+", required=True, help="reactivity.csv (transcript,position,base,reactivity)")
    parser.add_argument("--labels", nargs="+", required=True, help="One panel label per --inputs")
    parser.add_argument("--xrna", required=True, help="P4-P6 XRNA structure drawing")
    parser.add_argument("--transcript", default="p4p6", help="[default = p4p6] Transcript name in reactivity.csv")
    parser.add_argument("--offset", type=int, default=2,
                        help="[default = 2] reactivity position minus XRNA nucleotide number")
    parser.add_argument("--thresholds", default="0.3,0.7",
                        help="[default = 0.3,0.7] low/medium/high reactivity bin edges (black/orange/red)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if len(args.inputs) != len(args.labels):
        parser.error("--inputs and --labels must be the same length")

    thresholds = [float(v) for v in args.thresholds.split(",")]
    if len(thresholds) != 2:
        parser.error("--thresholds needs exactly two values")
    sequence = xrna_sequence(args.xrna)

    samples = []
    for csv, label in zip(args.inputs, args.labels):
        profile, n = profile_for(csv, args.transcript, sequence, args.offset, thresholds)
        print(f"{label}: {n} of {len(sequence)} P4-P6 nucleotides have a reactivity")
        samples.append(rnav.Sample(sample=label, ss=args.xrna, dms=profile))

    # Residues as outlined circles filled by reactivity bin, letters in a
    # contrasting color on top.
    rnav.styles.settings["ss"]["nucleotides"].update({"edgecolor": "black", "linewidths": 0.4})
    plot = rnav.plot_ss(
        samples=samples,
        structure="ss",
        profile="dms",
        colors={"sequence": "contrast", "nucleotides": "profile", "structure": "grey", "basepairs": None},
        nt_ticks=None,
        colorbars=False,
        plot_kwargs={"cols": min(len(samples), 3)},
    )
    structure = rnav.data.SecondaryStructure(input_data=args.xrna)
    for ax in plot.axes.flat[: len(samples)]:
        draw_basepairs(ax, structure)

    # RNAvigate draws its colorbar as a separate figure, so put the bin legend
    # on the structure figure itself.
    low, high = thresholds
    labels = [f"< {low:g}", f"{low:g} - {high:g}", f">= {high:g}"]
    handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor=c, markeredgecolor="black",
                      markeredgewidth=0.4, markersize=9, label=l)
               for c, l in zip(BIN_COLORS + ["grey"], labels + ["no data (G/U, NA)"])]
    plot.fig.legend(handles=handles, title="DMS reactivity", loc="lower center",
                    ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, -0.02))
    plot.fig.savefig(args.output, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
