#!/usr/bin/env python3
"""
Compare reactivity between two genotypes over matched transcripts (plotnine figure).

Inputs are the final reactivity.csv of every reactivity ID (transcript, position, base,
reactivity), each ID's genotype, and match_transcripts.py's position map. Per genotype, reactivity
is the mean over its IDs at each position (where defined). Per matched pair, positions are mapped
through the gap-free aligned blocks and kept where both genotypes have a reactivity (A/C in both).

Writes:
  --per-transcript  one row per pair with >= --min-positions shared positions: n_positions,
                    pearson, spearman, mean_abs_delta, mean reactivity per genotype, identity
  --per-position    (gzip TSV) every shared position: both names/positions/bases, the two genotype
                    means and their difference (b - a)
  --summary         median per-transcript correlation for each comparison below
  --plot            distribution of per-transcript Pearson r for: genotype means (a vs b); and,
                    as a like-for-like baseline, single replicates -- the first two IDs within each
                    genotype, and the first ID of a vs the first ID of b
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from plotnine import (aes, element_blank, geom_density, geom_vline, ggplot, labs, scale_color_manual, theme,
                      theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure



def read_reactivity(path, keep):
    parts = []
    for chunk in pd.read_csv(path, header=0, names=["transcript", "position", "base", "reactivity"],
                             dtype={"transcript": str, "base": str}, na_values=["NA"], chunksize=5_000_000):
        chunk = chunk[chunk["reactivity"].notna() & chunk["transcript"].isin(keep)]
        if len(chunk):
            parts.append(chunk)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["transcript", "position", "base", "reactivity"])
    out = {}
    for name, g in df.groupby("transcript", sort=False):
        out[name] = (g["position"].to_numpy(np.int64) - 1, g["reactivity"].to_numpy(np.float64), g["base"].to_numpy())
    return out


def dense(entry, length):
    """Reactivity vector (NaN where undefined) and bases over 0-based positions."""
    v = np.full(length, np.nan)
    b = np.full(length, "", dtype=object)
    if entry is not None:
        pos, val, base = entry
        ok = (pos >= 0) & (pos < length)
        v[pos[ok]], b[pos[ok]] = val[ok], base[ok]
    return v, b


def pearson(x, y):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    return pearson(pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy())


def block_index(blocks):
    ia, ib = [], []
    for block in str(blocks).split(";"):
        if not block or block == "nan":
            continue
        a0, b0, n = map(int, block.split(","))
        ia.append(np.arange(a0, a0 + n))
        ib.append(np.arange(b0, b0 + n))
    if not ia:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(ia), np.concatenate(ib)


def main():
    warnings.filterwarnings("ignore", message="Mean of empty slice")   # positions with no reactivity in any ID
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--labels", nargs="+", required=True)
    ap.add_argument("--genotypes", nargs="+", required=True, help="Genotype of each --inputs")
    ap.add_argument("--genotype-a", required=True)
    ap.add_argument("--genotype-b", required=True)
    ap.add_argument("--liftover", required=True)
    ap.add_argument("--min-positions", type=int, default=20)
    ap.add_argument("--per-transcript", required=True)
    ap.add_argument("--per-position", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--plot", required=True)
    args = ap.parse_args()
    if not len(args.inputs) == len(args.labels) == len(args.genotypes):
        ap.error("--inputs, --labels and --genotypes must be the same length")
    ga, gb = args.genotype_a, args.genotype_b

    lift = pd.read_csv(args.liftover, sep="\t", keep_default_na=False)
    lift = lift[lift["method"] != "too_long"]
    keep = {ga: set(lift["transcript_a"]), gb: set(lift["transcript_b"])}
    react = {}
    for path, label, g in zip(args.inputs, args.labels, args.genotypes):
        react[label] = read_reactivity(path, keep[g])
        print(f"{label} ({g}): {len(react[label]):,} matched transcripts with reactivity")
    ids = {g: [lab for lab, gg in zip(args.labels, args.genotypes) if gg == g] for g in (ga, gb)}
    if not ids[ga] or not ids[gb]:
        sys.exit(f"need reactivity IDs from both genotypes, got {ids}")

    rows, pos_rows, dist = [], [], []
    for r in lift.itertuples(index=False):
        ta, tb, la, lb = r.transcript_a, r.transcript_b, int(r.length_a), int(r.length_b)
        if not any(ta in react[i] for i in ids[ga]) or not any(tb in react[i] for i in ids[gb]):
            continue
        va = {i: dense(react[i].get(ta), la) for i in ids[ga]}
        vb = {i: dense(react[i].get(tb), lb) for i in ids[gb]}
        with np.errstate(all="ignore"):
            mean_a = np.nanmean(np.vstack([v for v, _ in va.values()]), axis=0)
            mean_b = np.nanmean(np.vstack([v for v, _ in vb.values()]), axis=0)
        base_a = va[ids[ga][0]][1]
        base_b = vb[ids[gb][0]][1]
        for i in ids[ga][1:]:
            base_a = np.where(base_a == "", va[i][1], base_a)
        for i in ids[gb][1:]:
            base_b = np.where(base_b == "", vb[i][1], base_b)
        ia, ib = block_index(r.blocks)
        ok = (ia < la) & (ib < lb)
        ia, ib = ia[ok], ib[ok]
        both = ~np.isnan(mean_a[ia]) & ~np.isnan(mean_b[ib])
        ia, ib = ia[both], ib[both]
        if len(ia) < args.min_positions:
            continue
        x, y = mean_a[ia], mean_b[ib]
        rows.append({"transcript_a": ta, "transcript_b": tb, "identity": r.identity, "n_positions": len(ia),
                     "pearson": pearson(x, y), "spearman": spearman(x, y), "mean_abs_delta": float(np.mean(np.abs(y - x))),
                     f"mean_reactivity_{ga}": float(np.mean(x)), f"mean_reactivity_{gb}": float(np.mean(y))})
        pos_rows.append(pd.DataFrame({"transcript_a": ta, "position_a": ia + 1, "base_a": base_a[ia],
                                      "transcript_b": tb, "position_b": ib + 1, "base_b": base_b[ib],
                                      f"reactivity_{ga}": np.round(x, 4), f"reactivity_{gb}": np.round(y, 4),
                                      "delta": np.round(y - x, 4)}))
        dist.append(("genotype means: " + f"{ga} vs {gb}", rows[-1]["pearson"]))
        # single-replicate baselines on the same transcripts
        a1, b1 = va[ids[ga][0]][0], vb[ids[gb][0]][0]
        m = ~np.isnan(a1[ia]) & ~np.isnan(b1[ib])
        if m.sum() >= args.min_positions:
            dist.append((f"one replicate each: {ids[ga][0]} vs {ids[gb][0]}", pearson(a1[ia][m], b1[ib][m])))
        for g, vv, idx in ((ga, va, ia), (gb, vb, ib)):
            if len(ids[g]) >= 2:
                r1, r2 = vv[ids[g][0]][0][idx], vv[ids[g][1]][0][idx]
                m = ~np.isnan(r1) & ~np.isnan(r2)
                if m.sum() >= args.min_positions:
                    dist.append((f"within {g}: {ids[g][0]} vs {ids[g][1]}", pearson(r1[m], r2[m])))

    per_tx = pd.DataFrame(rows)
    per_tx.to_csv(args.per_transcript, sep="\t", index=False, float_format="%.4f")
    (pd.concat(pos_rows, ignore_index=True) if pos_rows else pd.DataFrame()).to_csv(
        args.per_position, sep="\t", index=False, compression="gzip")
    d = pd.DataFrame(dist, columns=["comparison", "pearson"]).dropna()
    summary = d.groupby("comparison", sort=False)["pearson"].agg(transcripts="size", median_pearson="median").reset_index()
    summary.to_csv(args.summary, sep="\t", index=False, float_format="%.4f")
    print(f"{len(per_tx):,} matched transcripts with >= {args.min_positions} shared positions")
    print(summary.to_string(index=False))
    if d.empty:
        sys.exit("no matched transcript has enough shared positions to compare")

    order = summary["comparison"].tolist()
    d["comparison"] = pd.Categorical(d["comparison"], categories=order)
    palette = ["#14823c", "#d42a1e", "#4c72b0", "#8c8c8c", "#e3a21a"]
    p = (ggplot(d, aes("pearson", color="comparison"))
         + geom_density(size=1)
         + geom_vline(data=summary.assign(comparison=pd.Categorical(summary["comparison"], categories=order)),
                      mapping=aes(xintercept="median_pearson", color="comparison"), linetype="dashed")
         + scale_color_manual(values=palette[: len(order)], name="")
         + labs(x=f"Per-transcript Pearson r of reactivity (A/C positions, >= {args.min_positions} per transcript)",
                y="Density", title=f"Reactivity agreement between genotypes ({len(per_tx):,} matched transcripts)")
         + theme_bw(base_size=10)
         + theme(figure_size=(8, 4.5), panel_grid_minor=element_blank(), legend_position="bottom",
                 legend_direction="vertical", legend_title=element_blank()))
    p.save(args.plot, dpi=300, verbose=False)


if __name__ == "__main__":
    main()
