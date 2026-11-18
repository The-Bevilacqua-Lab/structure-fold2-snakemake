#!/usr/bin/env python3
"""
Metagene coverage profile for ONE sample.

For every transcript, per-position coverage c is transformed to
    ln(c + 1) / mean_over_transcript(ln(c + 1))
(transcripts whose mean is 0, i.e. no coverage at all, are skipped), then each
transcript is split into --bins equal-width bins of relative position (5'->3',
0-100% of its length) and averaged within each bin. Finally the bins are
averaged across transcripts, so every transcript carries equal weight
regardless of its length.

Coverage comes from one of two sources (--mode):
  depth  `samtools depth -a` output (chrom, pos, depth), read from --input
         (default "-" = stdin). Requires -a so every position of every
         transcript is present.
  rtsc   a StructureFold2 <.rtsc> file: per transcript, an ID line then a line of
         tab-separated RT-stop counts (one per position), then a blank line.

Output TSV columns: bin, position_pct (bin centre, 0-100), mean_norm_coverage,
n_transcripts.
"""

import argparse
import sys

import numpy as np
import pandas as pd

CHUNK_ROWS = 5_000_000


class Metagene:
    def __init__(self, nbins):
        self.nbins = nbins
        self.total = np.zeros(nbins)
        self.count = np.zeros(nbins, dtype=np.int64)

    def add(self, values, starts):
        """values: concatenated per-position coverage of consecutive
        transcripts; starts: index in values where each transcript begins."""
        starts = np.asarray(starts, dtype=np.int64)
        lengths = np.diff(np.append(starts, len(values)))
        n_tx = len(starts)

        ln_cov = np.log1p(values)
        means = np.add.reduceat(ln_cov, starts) / lengths
        valid = means > 0
        if not valid.any():
            return

        tx_idx = np.repeat(np.arange(n_tx), lengths)
        keep = valid[tx_idx]
        norm = ln_cov[keep] / means[tx_idx[keep]]

        rel_idx = np.arange(len(values)) - np.repeat(starts, lengths)
        bins = np.minimum(rel_idx * self.nbins // np.repeat(lengths, lengths), self.nbins - 1)

        key = tx_idx[keep] * self.nbins + bins[keep]
        size = n_tx * self.nbins
        sums = np.bincount(key, weights=norm, minlength=size).reshape(n_tx, self.nbins)
        cnts = np.bincount(key, minlength=size).reshape(n_tx, self.nbins)

        filled = cnts > 0
        per_tx_bin_mean = np.where(filled, sums / np.maximum(cnts, 1), 0.0)
        self.total += per_tx_bin_mean.sum(axis=0)
        self.count += filled.sum(axis=0)


def run_depth(handle, meta):
    """Stream `samtools depth -a` output; transcripts are contiguous rows."""
    carry_name, carry_vals = None, np.empty(0)
    reader = pd.read_csv(handle, sep="\t", header=None, usecols=[0, 2],
                         names=["contig", "_pos", "depth"], dtype={"contig": str, "depth": float},
                         chunksize=CHUNK_ROWS)
    for chunk in reader:
        names = chunk["contig"].to_numpy()
        vals = chunk["depth"].to_numpy()
        if carry_name is not None:
            names = np.concatenate([np.full(len(carry_vals), carry_name, dtype=object), names])
            vals = np.concatenate([carry_vals, vals])
        starts = np.append(0, np.flatnonzero(names[1:] != names[:-1]) + 1)
        last = starts[-1]
        carry_name, carry_vals = names[last], vals[last:]
        if last > 0:
            meta.add(vals[:last], starts[:-1])
    if carry_name is not None and len(carry_vals):
        meta.add(carry_vals, [0])


def run_rtsc(path, meta):
    """Read alternating ID / counts lines, batching transcripts for speed."""
    batch, batch_rows = [], 0

    def flush():
        nonlocal batch, batch_rows
        if batch:
            lengths = [len(b) for b in batch]
            starts = np.cumsum([0] + lengths[:-1])
            meta.add(np.concatenate(batch), starts)
        batch, batch_rows = [], 0

    # Records are: transcript ID line, counts line, blank line -- so read the
    # non-blank lines in (ID, counts) pairs.
    with open(path) as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        for _header, counts in zip(lines, lines):
            arr = np.array(counts.split("\t"), dtype=float)
            batch.append(arr)
            batch_rows += len(arr)
            if batch_rows >= CHUNK_ROWS:
                flush()
    flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", required=True, choices=["depth", "rtsc"])
    parser.add_argument("--input", default="-", help="samtools depth output or .rtsc file ('-' = stdin, depth mode only)")
    parser.add_argument("--bins", type=int, default=100, help="[default = 100] Number of relative-position bins")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    meta = Metagene(args.bins)
    if args.mode == "depth":
        handle = sys.stdin if args.input == "-" else args.input
        run_depth(handle, meta)
    else:
        if args.input == "-":
            parser.error("--mode rtsc needs a file path for --input")
        run_rtsc(args.input, meta)

    out = pd.DataFrame({
        "bin": np.arange(1, args.bins + 1),
        "position_pct": (np.arange(args.bins) + 0.5) * 100.0 / args.bins,
        "mean_norm_coverage": np.where(meta.count > 0, meta.total / np.maximum(meta.count, 1), np.nan),
        "n_transcripts": meta.count,
    })
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
