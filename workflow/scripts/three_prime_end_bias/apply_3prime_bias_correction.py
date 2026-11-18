#!/usr/bin/env python3
"""
Apply a fitted 3' bias model (fit_3prime_bias_model.py) to +DMS RT-stops.

Every transcript in the FASTA (mRNA, ncRNA and the positive control alike) has
the count at .rtsc index i divided by the model's fitted bias exp(prediction),
with position p = i and the transcript's FASTA length L (linear or spline form,
whichever model.tsv holds). The model is trained on protein-coding transcripts
only, but once 3' bias correction is on, every reactivity downstream is built
from corrected stops. Transcripts missing from the FASTA are written unchanged.

The corrected counts are then multiplied by one global constant so the
corrected total equals the raw one. This keeps the +DMS counts
on their own library scale -- the model intercept mostly reflects the
+DMS/-DMS library-size ratio -- which matters because rtsc_to_react takes
ln(count + 1) before normalizing. It doesn't change any position's count
relative to another's.

Also writes a before/after RT-stop metagene (training transcripts, per-
transcript mean-normalized, .rtsc index 0 dropped) against the -DMS reference.
"""

import argparse

import numpy as np
import pandas as pd

from bias_model import metagene, predict_log, read_model
from plot_bias import plot_rtsc


def read_rtsc(path):
    with open(path) as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        yield from zip(lines, lines)


def fasta_lengths(path):
    lengths, name, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def fit_factor(coefs, n, length):
    return np.exp(predict_log(coefs, np.arange(n), length))


def fmt(v):
    return "0" if v == 0 else f"{v:.6g}"


def long_table(vectors, lengths):
    """{transcript: counts} -> arrays for metagene(), dropping index 0."""
    names = list(vectors)
    t = np.concatenate([np.full(len(vectors[n]) - 1, n, dtype=object) for n in names])
    p = np.concatenate([np.arange(1, len(vectors[n])) for n in names])
    L = np.concatenate([np.full(len(vectors[n]) - 1, lengths[n]) for n in names])
    v = np.concatenate([vectors[n][1:] for n in names])
    return t, p, L, v


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True)
    parser.add_argument("--plus-rtsc", required=True)
    parser.add_argument("--minus-rtsc", required=True, help="-DMS reference .rtsc (plot only)")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--training", required=True)
    parser.add_argument("--bins", type=int, default=50)
    parser.add_argument("--title", default="")
    parser.add_argument("--output", required=True, help="Corrected +DMS .rtsc")
    parser.add_argument("--binned", required=True)
    parser.add_argument("--plot", required=True)
    args = parser.parse_args()

    coefs, model = read_model(args.model)
    lengths = fasta_lengths(args.fasta)
    training = set(pd.read_csv(args.training, sep="\t", dtype={"transcript": str})["transcript"])

    # Pass 1: totals for the library-scale constant.
    raw_total = corrected_total = 0.0
    for name, stops in read_rtsc(args.plus_rtsc):
        if name in lengths:
            counts = np.array(stops.split("\t"), dtype=float)
            raw_total += counts.sum()
            corrected_total += (counts / fit_factor(coefs, len(counts), lengths[name])).sum()
    scale = raw_total / corrected_total if corrected_total > 0 else 1.0
    print(f"Library-scale constant: {scale:.6g} (raw total {raw_total:,.0f})")

    # Pass 2: write, keeping the training transcripts for the plot.
    raw_train, corr_train, n_corrected, n_passed = {}, {}, 0, 0
    with open(args.output, "w") as out:
        for name, stops in read_rtsc(args.plus_rtsc):
            if name in lengths:
                counts = np.array(stops.split("\t"), dtype=float)
                corrected = counts / fit_factor(coefs, len(counts), lengths[name]) * scale
                stops = "\t".join(map(fmt, corrected))
                n_corrected += 1
                if name in training:
                    raw_train[name], corr_train[name] = counts, corrected
            else:
                n_passed += 1
            out.write(f"{name}\n{stops}\n\n")
    print(f"Corrected {n_corrected:,} transcripts; {n_passed:,} not in the FASTA written unchanged")

    minus_train = {name: np.array(stops.split("\t"), dtype=float)
                   for name, stops in read_rtsc(args.minus_rtsc) if name in raw_train}
    names = [n for n in raw_train if n in minus_train]
    meta = {}
    for key, vectors in [("plus_raw", raw_train), ("plus_corrected", corr_train), ("minus", minus_train)]:
        t, p, L, v = long_table({n: vectors[n] for n in names}, lengths)
        meta[key] = metagene(t, p, L, v, args.bins)

    x = np.arange(1, args.bins + 1)
    binned = pd.DataFrame(index=pd.Index(x, name="bin"))
    for key, (mean, sem) in meta.items():
        binned[f"rtsc_{key}_mean"], binned[f"rtsc_{key}_sem"] = mean, sem
    binned.to_csv(args.binned, sep="\t")
    plot_rtsc(binned, model, args.title, args.plot, n_training=len(names))

if __name__ == "__main__":
    main()
