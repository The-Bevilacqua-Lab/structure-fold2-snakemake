#!/usr/bin/env python3
"""
Pick the training set for the 3' bias model, shared by every reactivity ID.
Each protein-coding transcript is scored by its LOWEST +DMS RT-stop coverage
across the given +DMS .rtsc files (one per ID), then (--mode):

  top_n        the top --n transcripts by that lowest coverage (default), so
               every selected transcript is well covered in every replicate;
  all_covered  every transcript whose coverage is >= --min-coverage in ALL of
               the .rtsc files (no cap on the number).

Coverage is computed exactly as StructureFold's rtsc_coverage.py (and so
coverage.csv) does it: RT-stops on A/C bases divided by the number of A/C
bases, with .rtsc index i+1 paired with base i (the last base has no stop).
Only transcripts in --protein-coding are considered.

With --groups (one per label; genotype mode), a transcript named "<group><sep><id>" is scored
only over the labels of its own group, and top_n keeps the top --n of each group.

Outputs a TSV (transcript, length, min_coverage, coverage_<label>...; NA for labels of another
group) and a BED of whole transcripts (for `samtools depth -b`).
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
    parser.add_argument("--mode", choices=["top_n", "all_covered"], default="top_n",
                        help="top_n: top --n by lowest coverage; all_covered: every transcript "
                             "with coverage >= --min-coverage in every .rtsc")
    parser.add_argument("--n", type=int, help="Number of transcripts to keep (--mode top_n)")
    parser.add_argument("--min-coverage", type=float, default=1.0,
                        help="Coverage every .rtsc must reach (--mode all_covered; default 1)")
    parser.add_argument("--groups", nargs="+", help="Group (genotype) of each --rtsc label")
    parser.add_argument("--group-sep", default="__", help="Separator between group and transcript ID")
    parser.add_argument("--output", required=True, help="TSV: transcript, length, min_coverage, per-label coverage")
    parser.add_argument("--bed", required=True)
    args = parser.parse_args()
    if len(args.rtsc) != len(args.labels):
        parser.error("--rtsc and --labels must be the same length")
    if args.mode == "top_n" and (args.n is None or args.n < 1):
        parser.error("--mode top_n needs --n >= 1")

    with open(args.protein_coding) as f:
        pc = {line.strip() for line in f if line.strip()}
    seqs = read_fasta(args.fasta, pc)

    # The same .rtsc can back several IDs (e.g. pool_replicates: plus) -- read it once.
    by_file = {path: coverages(path, seqs) for path in dict.fromkeys(args.rtsc)}
    per_label = [by_file[path] for path in args.rtsc]

    if args.groups and len(args.groups) != len(args.labels):
        parser.error("--groups and --labels must be the same length")

    def group_of(name):
        return name.split(args.group_sep, 1)[0] if args.groups else None

    rows = []
    for name in seqs:
        idx = [i for i in range(len(per_label)) if not args.groups or args.groups[i] == group_of(name)]
        if not idx:
            continue
        covs = [per_label[i].get(name, 0.0) if i in idx else None for i in range(len(per_label))]
        low = min(covs[i] for i in idx)
        if low > 0:
            rows.append((name, len(seqs[name]), low, covs))
    rows.sort(key=lambda r: -r[2])
    if args.mode == "top_n":
        top, taken = [], {}
        for r in rows:                      # top --n overall, or of each group
            g = group_of(r[0])
            if taken.get(g, 0) < args.n:
                top.append(r)
                taken[g] = taken.get(g, 0) + 1
    else:
        top = [r for r in rows if r[2] >= args.min_coverage]
    if not top:
        need = "> 0" if args.mode == "top_n" else f">= {args.min_coverage:g}"
        sys.exit(f"No protein-coding transcript has +DMS RT-stop coverage {need} in every ID -- nothing to train on")
    not_pc = [name for name, *_ in top if name not in pc]
    if not_pc:
        sys.exit(f"BUG: {len(not_pc)} selected transcripts are not protein-coding, e.g. {not_pc[:5]}")
    if args.mode == "top_n" and not args.groups and len(top) < args.n:
        print(f"WARNING: only {len(top):,} protein-coding transcripts are covered in every ID (asked for {args.n:,})")

    with open(args.output, "w") as out, open(args.bed, "w") as bed:
        out.write("\t".join(["transcript", "length", "min_coverage"] + [f"coverage_{l}" for l in args.labels]) + "\n")
        for name, length, low, covs in top:
            out.write("\t".join([name, str(length), f"{low:.6g}"] + ["NA" if c is None else f"{c:.6g}" for c in covs]) + "\n")
            bed.write(f"{name}\t0\t{length}\n")
    rule_txt = f"top {args.n:,}" if args.mode == "top_n" else f"all with coverage >= {args.min_coverage:g}"
    print(f"Selected {len(top):,} ({rule_txt}) of {len(rows):,} protein-coding transcripts covered in all of "
          f"{', '.join(args.labels)}; lowest selected min-coverage = {top[-1][2]:.4g}")


if __name__ == "__main__":
    main()
