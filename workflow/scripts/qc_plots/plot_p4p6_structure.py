#!/usr/bin/env python3
"""
Plot p4p6 positive-control DMS reactivity onto the Tetrahymena P4-P6 domain secondary structure.

Structure and drawing coordinates come from an XRNA file (resources/p4p6/p4p6.xrna). Reactivity
position p of the p4p6 transcript maps to XRNA nucleotide p - --offset (the construct carries 2
extra 5' nucleotides the drawing doesn't). One panel per --inputs file, labelled by --labels.

Drawing style: each residue is a filled circle (no outline) coloured by reactivity bin -- black
(low), gold (medium), brick red (high), grey (no reactivity: G/U, NA) -- with its letter in white on
dark fills and dark on light fills. A thin dark line traces the backbone; base pairs are drawn as
two lines (G-C), one line (A-U), a filled dot (G-U) or an open dot (any other pair). Every 10th
residue gets a tick and its number (transcript position, i.e. XRNA number + --offset).

Drawn with matplotlib directly: circles are sized in structure units so neighbours sit edge to edge
at any figure size, which plotnine can't do.
"""

import argparse
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402  Arial in every figure


BIN_COLORS = ["#000000", "#f6cc4e", "#8a2a24"]   # low / medium / high reactivity
NO_DATA_COLOR = "#969696"                        # no reactivity (G/U, NA)
LINE_COLOR = "#222222"                           # backbone and base-pair symbols
NUMBER_COLOR = "#555555"
LETTER_LUMINANCE_CUTOFF = 0.2   # white letters on fills darker than this (black, red), dark on lighter
RADIUS = 0.47                   # circle radius, in units of the median backbone step
INCHES_PER_UNIT = 0.026         # drawing scale: inches per XRNA coordinate unit
NUMBER_EVERY = 10
LEGEND_TOP_LABEL = ">1.0"       # label above the top (high) segment of the colour bar


# ----------------------------------------------------------------------------- input
def read_xrna(path):
    """Sequence, x/y coordinates and base pairs (1-based (i, j), i < j) from an XRNA file."""
    text = open(path).read()
    block = text.split("DataType='NucChar.XPos.YPos'>")[1].split("<")[0]
    rows = [line.split() for line in block.strip().splitlines() if line.strip()]
    seq = "".join(r[0].upper().replace("T", "U") for r in rows)
    xy = np.array([[float(r[1]), float(r[2])] for r in rows])
    pairs = []
    for a, n, b in re.findall(r"<BasePairs nucID='(\d+)' length='(\d+)' bpNucID='(\d+)'", text):
        a, n, b = int(a), int(n), int(b)
        pairs += [(a + k, b - k) for k in range(n)]
    return seq, xy, pairs


def reactivity_for(csv, transcript, sequence, offset):
    """Reactivity per XRNA nucleotide (NaN where missing)."""
    df = pd.read_csv(csv)
    df.columns = ["transcript", "position", "base", "reactivity"]
    df = df[df["transcript"] == transcript]
    if df.empty:
        sys.exit(f"{csv}: no rows for transcript {transcript!r}")
    df = df.assign(nt=df["position"] - offset)
    df = df[(df["nt"] >= 1) & (df["nt"] <= len(sequence))].set_index("nt")
    mismatched = [i for i in df.index if df.at[i, "base"].upper().replace("T", "U") != sequence[i - 1]]
    if len(mismatched) > 2:
        print(f"WARNING: {csv}: {len(mismatched)} bases differ from the XRNA sequence "
              f"(XRNA positions {mismatched[:10]}) -- check --offset")
    return pd.to_numeric(df["reactivity"], errors="coerce").reindex(range(1, len(sequence) + 1)).to_numpy()


# ----------------------------------------------------------------------------- styling
def luminance(color):
    r, g, b = (c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(color))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def fill_colors(values, thresholds):
    low, high = thresholds
    out = []
    for v in values:
        if not np.isfinite(v):
            out.append(NO_DATA_COLOR)
        else:
            out.append(BIN_COLORS[0] if v < low else BIN_COLORS[1] if v < high else BIN_COLORS[2])
    return out


def letter_color(fill):
    return "white" if luminance(fill) < LETTER_LUMINANCE_CUTOFF else "#1a1a1a"


# ----------------------------------------------------------------------------- drawing
def draw_panel(ax, seq, xy, pairs, fills, r, pt_per_unit, title, offset):
    n = len(seq)
    partner = {i: j for i, j in pairs} | {j: i for i, j in pairs}
    # backbone
    ax.plot(xy[:, 0], xy[:, 1], color=LINE_COLOR, linewidth=0.8, zorder=1, solid_capstyle="round")
    # base pairs, trimmed to the circle edges
    for i, j in pairs:
        p0, p1 = xy[i - 1], xy[j - 1]
        d = p1 - p0
        length = np.hypot(*d)
        u = d / length
        a, b = p0 + u * r, p1 - u * r
        kind = "".join(sorted(seq[i - 1] + seq[j - 1]))
        if kind == "CG":                                   # G-C: two lines
            perp = np.array([-u[1], u[0]]) * r * 0.16
            for s in (-1, 1):
                ax.plot([a[0] + s * perp[0], b[0] + s * perp[0]], [a[1] + s * perp[1], b[1] + s * perp[1]],
                        color=LINE_COLOR, linewidth=0.9, zorder=1, solid_capstyle="butt")
        elif kind == "AU":                                 # A-U: one line
            ax.plot([a[0], b[0]], [a[1], b[1]], color=LINE_COLOR, linewidth=0.9, zorder=1, solid_capstyle="butt")
        else:                                              # G-U filled dot, anything else open dot
            mid = (p0 + p1) / 2
            ax.plot([a[0], b[0]], [a[1], b[1]], color=LINE_COLOR, linewidth=0.9, zorder=1, solid_capstyle="butt")
            ax.add_patch(Circle(mid, r * 0.28, facecolor=LINE_COLOR if kind == "GU" else "white",
                                edgecolor=LINE_COLOR, linewidth=0.9, zorder=2))
    # residues
    font = 1.15 * r * pt_per_unit
    for k in range(n):
        ax.add_patch(Circle(xy[k], r, facecolor=fills[k], edgecolor="none", zorder=3))
        ax.text(xy[k, 0], xy[k, 1], seq[k], ha="center", va="center_baseline", fontsize=font,
                color=letter_color(fills[k]), zorder=4, family=plot_style.FONT)
    # numbering every NUMBER_EVERY transcript positions: a short tick pointing away from the structure
    # locally, then the number anchored at the tick's end; moved out (or flipped) until clear of circles
    step = np.median(np.hypot(*np.diff(xy, axis=0).T))
    placed = []
    for k in range(n):
        number = k + 1 + offset
        if number % NUMBER_EVERY:
            continue
        if (k + 1) in partner:                             # paired: away from the partner
            out = xy[k] - xy[partner[k + 1] - 1]
        else:                                              # unpaired: away from its close neighbours
            nb = [xy[m] for m in (k - 1, k + 1) if 0 <= m < n and np.hypot(*(xy[m] - xy[k])) < 2.5 * step]
            out = xy[k] - np.mean(nb, axis=0) if nb else xy[k] - xy.mean(axis=0)
            if len(nb) == 2 and np.hypot(*out) < 0.3 * step:
                # on a straight run: perpendicular to the backbone, facing away from the centre
                tangent = nb[1] - nb[0]
                out = np.array([-tangent[1], tangent[0]])
                if np.dot(out, xy[k] - xy.mean(axis=0)) < 0:
                    out = -out
        out = out / np.hypot(*out)

        def clear(point):
            d_res = np.hypot(*(xy - point).T)
            d_res[k] = np.inf
            d_lab = [np.hypot(*(q - point)) for q in placed]
            return d_res.min() > 2.4 * r and (not d_lab or min(d_lab) > 4 * r)

        choice = None
        for direction in (out, -out):
            for reach in (2.3, 3.3, 4.5):
                tip = xy[k] + direction * r * reach
                label = tip + direction * r * 0.5
                if clear(label + direction * r * 1.2):     # label centre sits ~1 char beyond the anchor
                    choice = (direction, tip, label)
                    break
            if choice:
                break
        direction, tip, label = choice or (out, xy[k] + out * r * 2.3, xy[k] + out * r * 2.8)
        placed.append(label + direction * r * 1.2)
        t0 = xy[k] + direction * r * 1.2
        ax.plot([t0[0], tip[0]], [t0[1], tip[1]], color=NUMBER_COLOR, linewidth=0.8, zorder=1)
        ha = "left" if direction[0] > 0.35 else "right" if direction[0] < -0.35 else "center"
        va = "bottom" if direction[1] > 0.35 else "top" if direction[1] < -0.35 else "center"
        ax.text(label[0], label[1], str(number), ha=ha, va=va, fontsize=font * 1.45, color=NUMBER_COLOR, zorder=4)
    ax.set_title(title, fontsize=12, pad=4)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_legend(ax, thresholds):
    """Vertical colour bar: three equal segments (low / medium / high) labelled 0.0, the two cutoffs and
    LEGEND_TOP_LABEL, then a grey square for residues without a reactivity."""
    low, high = thresholds
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    edge = dict(edgecolor="#333333", linewidth=0.8)
    x, w, seg, top = 0.18, 0.28, 0.13, 0.80
    ax.text(0.0, 0.96, "Reactivity", ha="left", va="top", fontsize=15)
    ax.text(x + w / 2, top + 0.02, LEGEND_TOP_LABEL, ha="center", va="bottom", fontsize=13)
    for i, color in enumerate(reversed(BIN_COLORS)):          # high (top) to low (bottom)
        ax.add_patch(Rectangle((x, top - (i + 1) * seg), w, seg, facecolor=color, **edge))
    for i, value in enumerate((high, low)):                    # cutoffs at the inner boundaries
        ax.text(x + w + 0.08, top - (i + 1) * seg, f"{value:g}", ha="left", va="center", fontsize=13)
    bottom = top - 3 * seg
    ax.text(x + w / 2, bottom - 0.02, "0.0", ha="center", va="top", fontsize=13)
    square_y = bottom - 0.27                                # clear gap below the "0.0" label
    ax.add_patch(Rectangle((x, square_y), w, seg, facecolor=NO_DATA_COLOR, **edge))
    ax.text(x + w + 0.08, square_y + seg / 2, "G/U", ha="left", va="center", fontsize=13)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inputs", nargs="+", required=True, help="reactivity.csv (transcript,position,base,reactivity)")
    parser.add_argument("--labels", nargs="+", required=True, help="One panel label per --inputs")
    parser.add_argument("--xrna", required=True, help="P4-P6 XRNA structure drawing")
    parser.add_argument("--transcript", default="p4p6", help="[default = p4p6] Transcript name in reactivity.csv")
    parser.add_argument("--offset", type=int, default=2,
                        help="[default = 2] reactivity position minus XRNA nucleotide number")
    parser.add_argument("--thresholds", default="0.3,0.7",
                        help="[default = 0.3,0.7] low/medium/high reactivity bin edges (black/gold/red)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if len(args.inputs) != len(args.labels):
        parser.error("--inputs and --labels must be the same length")
    thresholds = [float(v) for v in args.thresholds.split(",")]
    if len(thresholds) != 2:
        parser.error("--thresholds needs exactly two values")

    seq, xy, pairs = read_xrna(args.xrna)
    step = np.median(np.hypot(*np.diff(xy, axis=0).T))
    r = RADIUS * step
    pad = r * 5                                            # room for circles and numbers
    (x0, y0), (x1, y1) = xy.min(axis=0) - pad, xy.max(axis=0) + pad
    w_in, h_in = (x1 - x0) * INCHES_PER_UNIT, (y1 - y0) * INCHES_PER_UNIT
    pt_per_unit = INCHES_PER_UNIT * 72

    cols = min(len(args.inputs), 3)
    rows = int(np.ceil(len(args.inputs) / cols))
    gap, title_h, legend_w = 0.3, 0.45, 1.6                # inches
    fig_w = cols * w_in + (cols + 1) * gap + legend_w
    fig_h = rows * (h_in + title_h) + gap
    fig = plt.figure(figsize=(fig_w, fig_h))
    for idx, (csv, label) in enumerate(zip(args.inputs, args.labels)):
        values = reactivity_for(csv, args.transcript, seq, args.offset)
        print(f"{label}: {int(np.isfinite(values).sum())} of {len(seq)} P4-P6 nucleotides have a reactivity")
        rr, cc = divmod(idx, cols)
        left = (gap + cc * (w_in + gap)) / fig_w
        bottom = (fig_h - (rr + 1) * (h_in + title_h)) / fig_h
        ax = fig.add_axes([left, bottom, w_in / fig_w, h_in / fig_h])
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        draw_panel(ax, seq, xy, pairs, fill_colors(values, thresholds), r, pt_per_unit, label, args.offset)

    legend_h = 3.4                                         # inches, centred on the right
    legend_ax = fig.add_axes([(fig_w - legend_w) / fig_w, (fig_h - legend_h) / 2 / fig_h,
                              (legend_w - 0.1) / fig_w, legend_h / fig_h])
    draw_legend(legend_ax, thresholds)
    fig.savefig(args.output, dpi=300, facecolor="white")


if __name__ == "__main__":
    main()
