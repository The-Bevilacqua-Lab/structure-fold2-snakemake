#!/usr/bin/env python3
"""
Depth-matched reactivity: every library thinned, per transcript, to the same
number of +DMS RT-stops, and reactivity averaged over many thinnings.
"""
import argparse
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

BASES = b"AC"  # bases that get a reactivity
CAP = 7.0  # reactivity ceiling
MIN_DEFINED = 0.5  # report a transcript only if at least this fraction of draws had a scale

# Filled in main(); read by the worker processes (fork).
PLUS = []  # per library: {transcript: +DMS stops}
MINUS = []  # per distinct -DMS file: {transcript: stops}
MINUS_OF = []  # per library: index into MINUS
SEQS = {}
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


def base_mask(seq, n):
    """True at positions 1..n whose base is in BASES (positions past the sequence end count as N)."""
    return np.isin(np.frombuffer(seq[:n].ljust(n, "N").encode(), dtype=np.uint8),
                   np.frombuffer(BASES, dtype=np.uint8))


def reactivity(plus, minus, mask):
    """Capped reactivity at positions 1..len(mask), 2-8% scaled over `mask`; None without a scale."""
    L = len(plus)
    short = min(L, len(minus))
    p, m = np.log1p(plus), np.log1p(minus)
    sp, sm = p.sum(), m.sum()
    if sp == 0 or sm == 0:
        return None
    d = np.maximum(0.0, p[:short] / sp * L - m[:short] / sm * L)[1:short]  # positions 1..short-1
    vals = d[mask]
    first, last = int(len(vals) * 0.02), int(len(vals) * 0.1)
    if last <= first:
        return None
    scale = np.sort(vals)[::-1][first:last].mean()
    if scale <= 0:
        return None
    return np.minimum(d / scale, CAP)


def process(task):
    idx, name = task
    rng = np.random.default_rng([ARGS.seed, idx])
    plus = [lib[name].astype(np.float64) for lib in PLUS]
    minus = [bg[name].astype(np.float64) for bg in MINUS]
    totals = [int(round(c.sum())) for c in plus]
    target = min(totals)
    if target <= 0:
        return None
    minus_target = min(int(round(c.sum())) for c in minus)
    masks = [base_mask(SEQS[name], min(len(c), len(minus[j])) - 1) for c, j in zip(plus, MINUS_OF)]

    sums = [np.zeros(len(m)) for m in masks]
    ndef = [0] * len(plus)
    draw1 = None
    for draw in range(ARGS.draws):
        bg = [thin(round_counts(c, rng), minus_target, rng) for c in minus] if ARGS.match_minus else minus
        matched = [thin(round_counts(c, rng), target, rng) for c in plus]
        if draw == 0:
            draw1 = matched
        for i, c in enumerate(matched):
            r = reactivity(c.astype(np.float64), bg[MINUS_OF[i]], masks[i])
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
    return name, target, totals, ndef, means, draw1


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
    ap.add_argument("--threads", type=int, default=1)
    ARGS = ap.parse_args()
    n = len(ARGS.labels)
    if len(ARGS.minus) == 1:
        ARGS.minus *= n
    for key in ("plus", "minus", "outdirs"):
        if len(getattr(ARGS, key)) != n:
            raise SystemExit(f"--{key} needs one entry per label" + (", or one for all" if key == "minus" else ""))

    # Transcripts in every library, with +DMS stops in each.
    names = None
    for path in ARGS.plus:
        d = {k: v for k, v in read_rtsc(path, names).items() if v.sum() > 0}
        names = set(d) if names is None else names & set(d)
        PLUS.append(d)
    minus_paths = sorted(set(ARGS.minus))
    for path in minus_paths:
        MINUS.append(read_rtsc(path, names))
        names &= set(MINUS[-1])
    MINUS_OF.extend(minus_paths.index(path) for path in ARGS.minus)
    SEQS = read_fasta(ARGS.fasta, names)
    names = sorted(names & set(SEQS))
    print(f"{len(names):,} transcripts in all {n} libraries, {ARGS.draws} draws each", flush=True)

    tasks = list(enumerate(names))
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
    for name, target, totals, ndef, means, draw1 in results:
        row = dict(transcript=name, target_plus=target)
        seq = SEQS[name]
        for i, label in enumerate(ARGS.labels):
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
    pd.DataFrame(rows).to_csv(ARGS.summary, sep="\t", index=False)
    print(f"{len(rows):,} transcripts matched; reported per library: {reported}")


if __name__ == "__main__":
    main()
