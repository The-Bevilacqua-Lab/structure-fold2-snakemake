#!/usr/bin/env python3
"""
Match transcripts between two genotypes and map positions between each matched pair.

Pairs: transcripts with the same ID in both genotypes, or the rows of --matching (a TSV with one
column per genotype name, e.g. "nipponbare<TAB>indica"; extra columns are ignored, a transcript may
appear in only one row). FASTA records are named "<genotype><sep><transcript_id>".

Positions: a pair of equal length whose ungapped identity is at least --ungapped-min-identity is
mapped 1:1 (the usual case for individuals of one species, where only SNPs differ). Any other pair
is aligned (Biopython PairwiseAligner; global, end gaps free so UTR length differences cost
nothing). Pairs whose sequence lengths multiply to more than --max-cells are left unaligned.

Writes one row per matched pair: transcript_a, transcript_b (full prefixed names), length_a,
length_b, method (ungapped / aligned / too_long), aligned_columns, identity (matches / aligned
columns), and blocks -- gap-free aligned runs as "start_a,start_b,length;..." with 0-based starts
(the same indexing as an .rtsc line).
"""
import argparse
import sys
from multiprocessing import Pool

import pandas as pd
from Bio import Align

ALIGNER = None


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def make_aligner():
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 1.0, -1.0
    a.open_gap_score, a.extend_gap_score = -3.0, -1.0
    a.target_end_gap_score = a.query_end_gap_score = 0.0
    return a


def init_worker():
    global ALIGNER
    ALIGNER = make_aligner()


def map_pair(job):
    name_a, seq_a, name_b, seq_b, min_ungapped, max_cells = job
    la, lb = len(seq_a), len(seq_b)
    if la == lb and la > 0:
        same = sum(x == y for x, y in zip(seq_a, seq_b))
        if same / la >= min_ungapped:
            return name_a, name_b, la, lb, "ungapped", la, same / la, f"0,0,{la}"
    if la * lb > max_cells:
        return name_a, name_b, la, lb, "too_long", 0, float("nan"), ""
    aln = ALIGNER.align(seq_a, seq_b)[0]
    blocks, cols, same = [], 0, 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        n = a1 - a0
        blocks.append(f"{a0},{b0},{n}")
        cols += n
        same += sum(seq_a[a0 + k] == seq_b[b0 + k] for k in range(n))
    return name_a, name_b, la, lb, "aligned", int(cols), (float(same) / cols if cols else float("nan")), ";".join(blocks)


def pairs_from_ids(seqs_a, seqs_b, ga, gb, sep):
    ids_a = {n[len(ga) + len(sep):]: n for n in seqs_a}
    ids_b = {n[len(gb) + len(sep):]: n for n in seqs_b}
    return [(ids_a[t], ids_b[t]) for t in ids_a if t in ids_b]


def pairs_from_table(path, seqs_a, seqs_b, ga, gb, sep):
    table = pd.read_csv(path, sep="\t", dtype=str)
    missing = [g for g in (ga, gb) if g not in table.columns]
    if missing:
        sys.exit(f"{path}: needs one column per genotype name {ga!r} and {gb!r}; missing {missing}")
    table = table[[ga, gb]].dropna()
    pairs, used_a, used_b, absent, dup = [], set(), set(), 0, 0
    for ta, tb in zip(table[ga].str.strip(), table[gb].str.strip()):
        na, nb = f"{ga}{sep}{ta}", f"{gb}{sep}{tb}"
        if na not in seqs_a or nb not in seqs_b:
            absent += 1
            continue
        if na in used_a or nb in used_b:
            dup += 1
            continue
        used_a.add(na)
        used_b.add(nb)
        pairs.append((na, nb))
    if absent:
        print(f"{absent:,} matching rows name a transcript missing from its genotype's transcriptome (skipped)")
    if dup:
        print(f"{dup:,} matching rows reuse an already-matched transcript (first row kept, rest skipped)")
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-a", required=True)
    ap.add_argument("--fasta-b", required=True)
    ap.add_argument("--genotype-a", required=True)
    ap.add_argument("--genotype-b", required=True)
    ap.add_argument("--sep", default="__")
    ap.add_argument("--matching", help="TSV with one column per genotype name (default: match identical IDs)")
    ap.add_argument("--ungapped-min-identity", type=float, default=0.8)
    ap.add_argument("--max-cells", type=float, default=1e8, help="[default = 1e8] skip aligning longer pairs")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    ga, gb, sep = args.genotype_a, args.genotype_b, args.sep
    seqs_a = {n: s for n, s in read_fasta(args.fasta_a).items() if n.startswith(ga + sep)}
    seqs_b = {n: s for n, s in read_fasta(args.fasta_b).items() if n.startswith(gb + sep)}
    pairs = (pairs_from_table(args.matching, seqs_a, seqs_b, ga, gb, sep) if args.matching
             else pairs_from_ids(seqs_a, seqs_b, ga, gb, sep))
    print(f"{ga}: {len(seqs_a):,} transcripts; {gb}: {len(seqs_b):,}; matched pairs: {len(pairs):,} "
          f"({'transcript_matching file' if args.matching else 'identical IDs'})")
    if not pairs:
        sys.exit("no transcripts matched between the genotypes -- check the IDs or provide transcript_matching")

    jobs = ((a, seqs_a[a], b, seqs_b[b], args.ungapped_min_identity, args.max_cells) for a, b in pairs)
    with Pool(args.threads, initializer=init_worker) as pool:
        rows = list(pool.imap(map_pair, jobs, chunksize=64))
    out = pd.DataFrame(rows, columns=["transcript_a", "transcript_b", "length_a", "length_b", "method",
                                      "aligned_columns", "identity", "blocks"])
    out.to_csv(args.output, sep="\t", index=False, float_format="%.4f")
    print(out["method"].value_counts().to_string())
    print(f"median identity {out['identity'].median():.4f}")


if __name__ == "__main__":
    main()
