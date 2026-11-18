#!/usr/bin/env python3
"""
Pool -DMS RT-stops across the two genotypes, position by position.

Input is the -DMS pooled within each genotype (one .rtsc holding both genotypes' transcripts under
their "<genotype>__<id>" names) and match_transcripts.py's pair list. For every matched pair, each
transcript gets its own stops plus the other genotype's stops at the matching positions, so both
members of the pair end up with the same pooled background, each on its own coordinates. Positions
with no partner (inside an indel) keep only their own stops.

Where the matching positions come from:

  alignment (default)  match_transcripts.py's gap-free aligned blocks. Pairs below --min-identity and
                       unaligned (too_long) pairs are written unchanged. Equal-length pairs are taken
                       as 1:1 without aligning, so a pair whose indels cancel in length is offset
                       between them.

  --coordinate-map     A variant-derived position map (TSV: transcript_id, ref_pos, alt_pos; 1-based;
                       -1 = no partner), e.g. the pseudo-transcriptome's coordinate_matching.tsv.gz.
                       ref_pos is in --map-ref-genotype's coordinates. Every matched pair whose two
                       members share a transcript ID is pooled, with no identity cutoff: transcripts
                       listed in the map use its rows, unlisted ones are 1:1 (and must be the same
                       length). A pair whose map positions fall outside either transcript is written
                       unchanged and counted.

Unmatched transcripts are written unchanged.
"""
import argparse

import numpy as np
import pandas as pd

SEP = "__"


def read_rtsc(path):
    out, order = {}, []
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    i = 0
    while i < len(lines):
        name = lines[i].strip()
        if not name:
            i += 1
            continue
        out[name] = np.array(lines[i + 1].split("\t"), dtype=np.float64) if lines[i + 1] else np.zeros(0)
        order.append(name)
        i += 3
    return out, order


def fmt(v):
    return str(int(v)) if float(v).is_integer() else f"{v:.6g}"


def split_name(name):
    genotype, _, tid = name.partition(SEP)
    return genotype, tid


def read_coordinate_map(path):
    """{transcript_id: (ref_idx, alt_idx)} as 0-based int arrays, rows with a partner on both sides."""
    m = pd.read_csv(path, sep="\t", usecols=["transcript_id", "ref_pos", "alt_pos"],
                    dtype={"transcript_id": "category", "ref_pos": np.int64, "alt_pos": np.int64})
    m = m[(m.ref_pos > 0) & (m.alt_pos > 0)]
    codes = m.transcript_id.cat.codes.to_numpy()
    ref, alt = m.ref_pos.to_numpy() - 1, m.alt_pos.to_numpy() - 1
    order = np.argsort(codes, kind="stable")
    codes, ref, alt = codes[order], ref[order], alt[order]
    bounds = np.flatnonzero(np.diff(codes)) + 1
    starts = np.r_[0, bounds]
    ends = np.r_[bounds, len(codes)]
    cats = m.transcript_id.cat.categories
    return {cats[codes[s]]: (ref[s:e], alt[s:e]) for s, e in zip(starts, ends)}


def pool_alignment(lift, stops, pooled, min_identity):
    n = {"pooled": 0, "below_identity_or_unaligned": 0, "missing_minus": 0}
    for row in lift.itertuples(index=False):
        a, b = row.transcript_a, row.transcript_b
        if a not in stops or b not in stops:
            n["missing_minus"] += 1
            continue
        identity = float(row.identity) if row.identity != "" else float("nan")
        if row.method == "too_long" or not identity >= min_identity:
            n["below_identity_or_unaligned"] += 1
            continue
        for block in row.blocks.split(";"):
            a0, b0, k = map(int, block.split(","))
            k = min(k, len(stops[a]) - a0, len(stops[b]) - b0)
            if k <= 0:
                continue
            pooled[a][a0:a0 + k] += stops[b][b0:b0 + k]
            pooled[b][b0:b0 + k] += stops[a][a0:a0 + k]
        n["pooled"] += 1
    return n


def pool_coordinate_map(lift, stops, pooled, cmap, ref_genotype):
    n = {"pooled_from_map": 0, "pooled_1to1": 0, "missing_minus": 0, "different_ids": 0,
         "unlisted_length_mismatch": 0, "map_out_of_range": 0}
    for row in lift.itertuples(index=False):
        a, b = row.transcript_a, row.transcript_b
        if a not in stops or b not in stops:
            n["missing_minus"] += 1
            continue
        (ga, ta), (gb, tb) = split_name(a), split_name(b)
        if ta != tb:
            n["different_ids"] += 1
            continue
        ref, alt = (a, b) if ga == ref_genotype else (b, a)
        if ta in cmap:
            ri, ai = cmap[ta]
            if len(ri) and (ri.max() >= len(stops[ref]) or ai.max() >= len(stops[alt])):
                n["map_out_of_range"] += 1
                continue
            n["pooled_from_map"] += 1
        else:
            if len(stops[ref]) != len(stops[alt]):
                n["unlisted_length_mismatch"] += 1
                continue
            ri = ai = np.arange(len(stops[ref]))
            n["pooled_1to1"] += 1
        pooled[ref][ri] += stops[alt][ai]
        pooled[alt][ai] += stops[ref][ri]
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--minus", required=True, help="-DMS .rtsc pooled within each genotype")
    ap.add_argument("--liftover", required=True, help="match_transcripts.py output (the matched pairs)")
    ap.add_argument("--min-identity", type=float, default=0.9, help="alignment mode only")
    ap.add_argument("--coordinate-map", help="variant-derived position map; replaces the alignment blocks")
    ap.add_argument("--map-ref-genotype", help="genotype whose coordinates are the map's ref_pos")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if args.coordinate_map and not args.map_ref_genotype:
        ap.error("--coordinate-map needs --map-ref-genotype")

    stops, order = read_rtsc(args.minus)
    pooled = {name: v.copy() for name, v in stops.items()}
    lift = pd.read_csv(args.liftover, sep="\t", keep_default_na=False)
    if args.coordinate_map:
        cmap = read_coordinate_map(args.coordinate_map)
        print(f"coordinate map: {args.coordinate_map} ({len(cmap):,} transcripts listed; "
              f"ref = {args.map_ref_genotype})")
        n = pool_coordinate_map(lift, stops, pooled, cmap, args.map_ref_genotype)
    else:
        print(f"alignment blocks: {args.liftover} (min identity {args.min_identity})")
        n = pool_alignment(lift, stops, pooled, args.min_identity)
    with open(args.output, "w") as out:
        for name in order:
            out.write(name + "\n" + "\t".join(map(fmt, pooled[name])) + "\n\n")
    print("pairs: " + "; ".join(f"{k.replace('_', ' ')}: {v:,}" for k, v in n.items()))


if __name__ == "__main__":
    main()
