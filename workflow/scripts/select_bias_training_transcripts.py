#!/usr/bin/env python3
"""
Pick the training set for the 3' bias model, shared by every reactivity ID:
the top --n protein-coding transcripts by +DMS RT-stop coverage, ranked by
each transcript's LOWEST coverage across the given +DMS .rtsc files (one per
ID), so every selected transcript is well covered in every replicate.

Coverage is computed exactly as StructureFold's rtsc_coverage.py (and so
coverage.csv) does it: RT-stops on A/C bases divided by the number of A/C
bases, with .rtsc index i+1 paired with base i (the last base has no stop).
Only transcripts in --protein-coding are considered.

Outputs a TSV (transcript, length, min_coverage, coverage_<label>...) and a
BED of whole transcripts (for `samtools depth -b`).
"""

import argparse
import sys

import numpy as np

SPECIFIC = np.zeros(256, dtype=bool)
SPECIFIC[[ord("A"), ord("C")]] = True


def read_fasta(path, keep):
    seqs, name, parts = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name in keep:
        seqs[name] = "".join(parts)
    return seqs


def read_rtsc(path):
    with open(path) as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        yield from zip(lines, lines)


def coverages(rtsc, seqs):
    out = {}
    for name, stops in read_rtsc(rtsc):
        seq = seqs.get(name)
        if seq is None:
            continue
        mask = SPECIFIC[np.frombuffer(seq[:-1].encode("ascii", "replace"), dtype=np.uint8)]
        n_specific = mask.sum()
        if n_specific == 0:
            continue
        counts = np.array(stops.split("\t")[1:], dtype=float)
        k = min(len(counts), len(mask))
        out[name] = counts[:k][mask[:k]].sum() / n_specific
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rtsc", nargs="+", required=True, help="+DMS .rtsc, one per ID")
    parser.add_argument("--labels", nargs="+", required=True, help="One label per --rtsc")
    parser.add_argument("--fasta", required=True, help="Transcriptome FASTA (uppercase)")
    parser.add_argument("--protein-coding", required=True, help="Protein-coding transcript IDs, one per line")
    parser.add_argument("--n", type=int, required=True, help="Number of transcripts to keep")
    parser.add_argument("--output", required=True, help="TSV: transcript, length, min_coverage, per-label coverage")
    parser.add_argument("--bed", required=True)
    args = parser.parse_args()
    if len(args.rtsc) != len(args.labels):
        parser.error("--rtsc and --labels must be the same length")

    with open(args.protein_coding) as f:
        pc = {line.strip() for line in f if line.strip()}
    seqs = read_fasta(args.fasta, pc)

    # The same .rtsc can back several IDs (e.g. pool_replicates: plus) -- read it once.
    by_file = {path: coverages(path, seqs) for path in dict.fromkeys(args.rtsc)}
    per_label = [by_file[path] for path in args.rtsc]

    rows = []
    for name in seqs:
        covs = [c.get(name, 0.0) for c in per_label]
        if min(covs) > 0:
            rows.append((name, len(seqs[name]), min(covs), covs))
    rows.sort(key=lambda r: -r[2])
    top = rows[: args.n]
    if not top:
        sys.exit("No protein-coding transcript has +DMS RT-stop coverage > 0 in every ID -- nothing to train on")
    not_pc = [name for name, *_ in top if name not in pc]
    if not_pc:
        sys.exit(f"BUG: {len(not_pc)} selected transcripts are not protein-coding, e.g. {not_pc[:5]}")
    if len(top) < args.n:
        print(f"WARNING: only {len(top):,} protein-coding transcripts are covered in every ID (asked for {args.n:,})")

    with open(args.output, "w") as out, open(args.bed, "w") as bed:
        out.write("\t".join(["transcript", "length", "min_coverage"] + [f"coverage_{l}" for l in args.labels]) + "\n")
        for name, length, low, covs in top:
            out.write("\t".join([name, str(length), f"{low:.6g}"] + [f"{c:.6g}" for c in covs]) + "\n")
            bed.write(f"{name}\t0\t{length}\n")
    print(f"Selected {len(top):,} of {len(rows):,} protein-coding transcripts covered in all of "
          f"{', '.join(args.labels)}; lowest selected min-coverage = {top[-1][2]:.4g}")


if __name__ == "__main__":
    main()
