#!/usr/bin/env python3
"""
Depth-matched reactivity: every library thinned, per transcript, to the same
number of +DMS RT-stops, and reactivity averaged over many thinnings.

--shard i --n-shards N processes only the i-th (0-based) of N contiguous blocks of the sorted
transcript list. Each transcript's random draws are seeded by its index in the full list, so
concatenating the N shards' outputs in shard order (merge_depth_matched_shards.py) gives exactly
the output of one unsharded run.

--normalization sets how each channel's ln(stops + 1) is put on a common scale before -DMS is subtracted:
  all  divided by its mean over all positions of the transcript (StructureFold2; default)
  gu   divided by its mean over the transcript's G/U positions, which DMS does not modify on the Watson-Crick
       face, so the mean measures background only (an all-position +DMS mean also contains the A/C DMS signal,
       which scales the +DMS background below the -DMS background in proportion to the signal). A library has no
       scale for a draw when the transcript has fewer than --min-gu G/U positions or either channel has no stops
       at G/U; libraries without a scale in more than half of the draws are not reported, as with "all".
"""
import argparse
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

BASES = b"AC"  # bases that get a reactivity
GU_BASES = b"GT"  # background bases for --normalization gu
CAP = 7.0  # reactivity ceiling
MIN_DEFINED = 0.5  # report a transcript only if at least this fraction of draws had a scale

# Filled in main(); read by the worker processes (fork).
PLUS = []  # per library: {transcript: +DMS stops}
MINUS = []  # per distinct -DMS file: {transcript: stops}
MINUS_OF = []  # per library: index into MINUS
SEQS = {}
NAME_OF = []  # per library: {matching key: transcript name in that library} (--key-map)
ARGS = None


def read_rtsc(path, keep=None):
    out = {}
    with open(path) as f:
        while True:
            name = f.readline()
            if not name:
                break
            stops = f.readline()
            f.readline()
            name = name.strip()
            if keep is None or name in keep:
                out[name] = np.array(stops.split("\t"), dtype=np.float32)
    return out


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf).upper()
    return seqs


def round_counts(counts, rng):
    """Integer counts with total round(counts.sum()), unbiased per position."""
    floor = np.floor(counts)
    frac = counts - floor
    out = floor.astype(np.int64)
    extra = int(round(counts.sum())) - int(out.sum())
    if extra > 0:
        nz = np.flatnonzero(frac > 0)
        extra = min(extra, len(nz))
        p = frac[nz] / frac[nz].sum()
        out[rng.choice(nz, size=extra, replace=False, p=p)] += 1
    return out


def thin(counts, total, rng):
    if counts.sum() <= total:
        return counts
    return rng.multivariate_hypergeometric(counts, total)


def base_mask(seq, n, bases=BASES):
    """True at positions 1..n whose base is in `bases` (positions past the sequence end count as N)."""
    return np.isin(np.frombuffer(seq[:n].ljust(n, "N").encode(), dtype=np.uint8),
                   np.frombuffer(bases, dtype=np.uint8))


def reactivity(plus, minus, mask, gu_mask=None, min_gu=20):
    """Capped reactivity at positions 1..len(mask), 2-8% scaled over `mask`; None without a scale.
    gu_mask (positions 1..len(mask)): G/U normalization -- each channel's ln(stops + 1) divided by its mean over
    the G/U positions instead of its mean over all positions."""
    L = len(plus)
    short = min(L, len(minus))
    p, m = np.log1p(plus), np.log1p(minus)
    if gu_mask is None:
        sp, sm = p.sum(), m.sum()
        if sp == 0 or sm == 0:
            return None
        d = np.maximum(0.0, p[:short] / sp * L - m[:short] / sm * L)[1:short]  # positions 1..short-1
    else:
        if gu_mask.sum() < min_gu:
            return None
        gp, gm = p[1:short][gu_mask].mean(), m[1:short][gu_mask].mean()
        if gp == 0 or gm == 0:
            return None
        d = np.maximum(0.0, p[1:short] / gp - m[1:short] / gm)  # positions 1..short-1
    vals = d[mask]
    first, last = int(len(vals) * 0.02), int(len(vals) * 0.1)
    if last <= first:
        return None
    scale = np.sort(vals)[::-1][first:last].mean()
    if scale <= 0:
        return None
    return np.minimum(d / scale, CAP)


def process(task):
    idx, key = task
    rng = np.random.default_rng([ARGS.seed, idx])
    real = [names.get(key, key) for names in NAME_OF]   # each library's own name for this transcript
    plus = [lib[key].astype(np.float64) for lib in PLUS]
    # one background per distinct (-DMS file, transcript name), in file order
    bgs = list(dict.fromkeys((f, real[i]) for f in range(len(MINUS)) for i in range(len(PLUS)) if MINUS_OF[i] == f))
    bg_of = [bgs.index((MINUS_OF[i], real[i])) for i in range(len(PLUS))]
    minus = [MINUS[f][r].astype(np.float64) for f, r in bgs]
    totals = [int(round(c.sum())) for c in plus]
    target = min(totals)
    if target <= 0:
        return None
    minus_target = min(int(round(c.sum())) for c in minus)
    masks = [base_mask(SEQS[real[i]], min(len(c), len(minus[bg_of[i]])) - 1) for i, c in enumerate(plus)]
    gu_masks = ([base_mask(SEQS[real[i]], len(m), GU_BASES) for i, m in enumerate(masks)]
                if ARGS.normalization == "gu" else [None] * len(plus))

    sums = [np.zeros(len(m)) for m in masks]
    ndef = [0] * len(plus)
    draw1 = None
    for draw in range(ARGS.draws):
        bg = [thin(round_counts(c, rng), minus_target, rng) for c in minus] if ARGS.match_minus else minus
        matched = [thin(round_counts(c, rng), target, rng) for c in plus]
        if draw == 0:
            draw1 = matched
        for i, c in enumerate(matched):
            r = reactivity(c.astype(np.float64), bg[bg_of[i]], masks[i], gu_masks[i], ARGS.min_gu)
            if r is not None:
                sums[i] += r
                ndef[i] += 1
    means = []
    for s, n, m in zip(sums, ndef, masks):
        if n >= MIN_DEFINED * ARGS.draws:
            r = s / n
            r[~m] = np.nan
            means.append(r)
        else:
            means.append(None)
    return key, real, target, totals, ndef, means, draw1


def main():
    global ARGS, SEQS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", nargs="+", required=True, help="One name per library")
    ap.add_argument("--plus", nargs="+", required=True, help="+DMS .rtsc per library (3' bias corrected)")
    ap.add_argument("--minus", nargs="+", required=True, help="-DMS .rtsc per library, or one for all")
    ap.add_argument("--fasta", required=True, help="Transcriptome FASTA")
    ap.add_argument("--outdirs", nargs="+", required=True, help="Output directory per library")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--draws", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20261014)
    ap.add_argument("--match-minus", action="store_true", help="Also match the -DMS backgrounds")
    ap.add_argument("--key-map", help="TSV (name, key): transcripts sharing a key are matched across libraries "
                                      "(e.g. orthologs of two genotypes); default: by name")
    ap.add_argument("--normalization", choices=("all", "gu"), default="all",
                    help="scale of ln(stops + 1) before the -DMS subtraction: mean over all positions (default) "
                         "or over G/U positions")
    ap.add_argument("--min-gu", type=int, default=20, help="--normalization gu: fewest G/U positions for a scale")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0, help="0-based shard index (with --n-shards)")
    ap.add_argument("--n-shards", type=int, default=1, help="split the transcripts into this many blocks")
    ARGS = ap.parse_args()
    if not 0 <= ARGS.shard < ARGS.n_shards:
        raise SystemExit("--shard must be in 0 .. --n-shards - 1")
    n = len(ARGS.labels)
    if len(ARGS.minus) == 1:
        ARGS.minus *= n
    for key in ("plus", "minus", "outdirs"):
        if len(getattr(ARGS, key)) != n:
            raise SystemExit(f"--{key} needs one entry per label" + (", or one for all" if key == "minus" else ""))

    key_of = {}
    if ARGS.key_map:
        km = pd.read_csv(ARGS.key_map, sep="\t", dtype=str)
        key_of = dict(zip(km["name"], km["key"]))

    # Transcripts (by matching key) in every library, with +DMS stops in each.
    names = None
    for path in ARGS.plus:
        d, own = {}, {}
        for k, v in read_rtsc(path).items():
            if v.sum() > 0:
                d[key_of.get(k, k)] = v
                own[key_of.get(k, k)] = k
        names = set(d) if names is None else names & set(d)
        PLUS.append(d)
        NAME_OF.append(own)
    minus_paths = sorted(set(ARGS.minus))
    MINUS_OF.extend(minus_paths.index(path) for path in ARGS.minus)
    needed = {own[k] for own in NAME_OF for k in names}
    for path in minus_paths:
        MINUS.append(read_rtsc(path, needed))
    for i, f in enumerate(MINUS_OF):
        names = {k for k in names if NAME_OF[i][k] in MINUS[f]}
    SEQS = read_fasta(ARGS.fasta, {own[k] for own in NAME_OF for k in names})
    names = sorted(k for k in names if all(own[k] in SEQS for own in NAME_OF))
    print(f"{len(names):,} transcripts in all {n} libraries, {ARGS.draws} draws each; normalization: "
          f"{ARGS.normalization}" + (f" (>= {ARGS.min_gu} G/U positions)" if ARGS.normalization == "gu" else ""), flush=True)

    tasks = list(enumerate(names))
    if ARGS.n_shards > 1:
        lo = len(tasks) * ARGS.shard // ARGS.n_shards
        hi = len(tasks) * (ARGS.shard + 1) // ARGS.n_shards
        tasks = tasks[lo:hi]
        print(f"shard {ARGS.shard + 1}/{ARGS.n_shards}: transcripts {lo + 1:,}-{hi:,} ({len(tasks):,})", flush=True)
    if ARGS.threads > 1:
        with Pool(ARGS.threads) as pool:
            results = [r for r in pool.imap(process, tasks, chunksize=20) if r is not None]
    else:
        results = [r for r in map(process, tasks) if r is not None]

    handles = []
    for outdir in ARGS.outdirs:
        os.makedirs(outdir, exist_ok=True)
        react = open(os.path.join(outdir, "reactivity.csv"), "w")
        react.write("transcript,position,base,reactivity\n")
        handles.append((react, open(os.path.join(outdir, "matched_plus_draw1.rtsc"), "w")))
    rows = []
    reported = dict.fromkeys(ARGS.labels, 0)
    for key, real, target, totals, ndef, means, draw1 in results:
        row = dict(transcript=key, target_plus=target)
        for i, label in enumerate(ARGS.labels):
            name, seq = real[i], SEQS[real[i]]
            row[f"{label}_plus_total"] = totals[i]
            row[f"{label}_draws_defined"] = ndef[i]
            react, rtsc = handles[i]
            rtsc.write(name + "\n" + "\t".join(map(str, draw1[i])) + "\n\n")
            if means[i] is None:
                continue
            reported[label] += 1
            text = ["NA" if v != v else f"{v:.3f}" for v in means[i]] + ["NA"]
            react.write("".join(f"{name},{p},{seq[p - 1] if p <= len(seq) else 'N'},{t}\n"
                                for p, t in enumerate(text, start=1)))
        rows.append(row)
    for react, rtsc in handles:
        react.close()
        rtsc.close()
    columns = ["transcript", "target_plus"] + [f"{label}_{c}" for label in ARGS.labels
                                               for c in ("plus_total", "draws_defined")]
    pd.DataFrame(rows, columns=columns).to_csv(ARGS.summary, sep="\t", index=False)
    print(f"{len(rows):,} transcripts matched; reported per library: {reported}")
    unscaled = sum(1 for *_, ndef, means, _ in results if any(m is None for m in means))
    print(f"{unscaled:,} transcripts not reported in at least one library (a scale in < {MIN_DEFINED:.0%} of draws)")


if __name__ == "__main__":
    main()
