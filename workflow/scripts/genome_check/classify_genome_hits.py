#!/usr/bin/env python3
"""
Place one sample's genome alignments (reads Bowtie2 could not align to the transcriptome) on the
genome annotation.

Input --alignments: gzip TSV of primary alignments, no header: read, flag, chrom, pos, cigar, NH.
Reads are sense to the RNA (the pipeline keeps forward-strand transcriptome alignments), so a read on
the genome's + strand comes from a + strand transcript. Each read is placed by its 5' end (the base
after the RT stop): pos for + strand reads, the last aligned base for - strand reads.

Category (first that applies):
  organelle (<chrom>)        chromosome matches --organelles (e.g. Mt, Pt)
  transcriptome exon         exon of a transcript that is in the transcriptome (should have aligned:
                             usually an unannotated splice junction, indels, or a mismatch at base 1
                             tolerated by HISAT2 but not Bowtie2)
  <biotype> exon             exon of an annotated transcript left out of the transcriptome
                             (other isoforms, ncRNAs dropped by cds_only, rRNA/tRNA, ...)
  intron                     inside a gene on the same strand, not in any of its exons
  antisense                  inside a gene (exon or intron) on the opposite strand only
  intergenic                 none of the above

Writes --positions (gzip TSV): one row per (chrom, strand, 5' end, category, gene): reads,
unique_reads (NH = 1); and --counts: one row per category: reads, unique_reads.
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd

CIGAR_REF = re.compile(r"(\d+)([MDN=X])")


def ref_span(cigar):
    return sum(int(n) for n, _ in CIGAR_REF.findall(cigar))


class Painter:
    """Disjoint segments of one chromosome strand, each labelled with its highest-priority interval."""

    def __init__(self, starts, ends, prio, label):
        bp = np.unique(np.concatenate([starts, ends + 1]))
        seg = np.full(len(bp), -1, dtype=np.int64)
        for i in np.argsort(prio, kind="stable"):
            a, b = np.searchsorted(bp, starts[i]), np.searchsorted(bp, ends[i] + 1)
            s = seg[a:b]
            s[s < 0] = label[i]
        self.bp, self.seg = bp, seg

    def lookup(self, pos):
        k = np.searchsorted(self.bp, pos, side="right") - 1
        out = np.full(len(pos), -1, dtype=np.int64)
        ok = k >= 0
        out[ok] = self.seg[k[ok]]
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignments", required=True)
    ap.add_argument("--intervals", required=True, help="annotation_intervals.py output")
    ap.add_argument("--organelles", default="")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--positions", required=True)
    ap.add_argument("--counts", required=True)
    args = ap.parse_args()

    aln = pd.read_csv(args.alignments, sep="\t", header=None, names=["read", "flag", "chrom", "pos", "cigar", "NH"],
                      dtype={"read": str, "chrom": str, "cigar": str})
    rev = (aln["flag"].to_numpy() & 16) > 0
    pos = aln["pos"].to_numpy(np.int64)
    if len(aln):
        span = np.fromiter((ref_span(c) for c in aln["cigar"]), dtype=np.int64, count=len(aln))
        pos = np.where(rev, pos + span - 1, pos)
    aln["strand"] = np.where(rev, "-", "+")
    aln["pos5"] = pos
    aln["unique"] = (aln["NH"].fillna(1).astype(int) == 1).astype(int)

    iv = pd.read_csv(args.intervals, sep="\t", dtype={"chrom": str, "gene_id": str, "gene_name": str,
                                                       "transcript_id": str, "biotype": str},
                     keep_default_na=False)
    if len(iv):
        iv["in_transcriptome"] = iv["in_transcriptome"].astype(int)
        iv["prio"] = np.select([(iv["feature"] == "exon") & (iv["in_transcriptome"] == 1), iv["feature"] == "exon"],
                               [0, 1], 2)
        iv["category"] = np.select([iv["prio"] == 0, iv["prio"] == 1],
                                   ["transcriptome exon", iv["biotype"] + " exon"], "intron")
    labels = iv[["category", "gene_id", "gene_name", "biotype"]].reset_index(drop=True) if len(iv) else None

    n = len(aln)
    sense = np.full(n, -1, dtype=np.int64)
    anti = np.full(n, -1, dtype=np.int64)
    if len(iv):
        for (chrom, strand), g in iv.groupby(["chrom", "strand"], sort=False):
            painter = Painter(g["start"].to_numpy(np.int64), g["end"].to_numpy(np.int64), g["prio"].to_numpy(),
                              g.index.to_numpy(np.int64))
            same = np.flatnonzero((aln["chrom"].to_numpy() == chrom) & (aln["strand"].to_numpy() == strand))
            other = np.flatnonzero((aln["chrom"].to_numpy() == chrom) & (aln["strand"].to_numpy() != strand))
            if len(same):
                sense[same] = painter.lookup(pos[same])
            if len(other):
                hit = painter.lookup(pos[other])
                anti[other] = np.where(anti[other] >= 0, anti[other], hit)
    cat = np.full(n, "intergenic", dtype=object)
    gene = np.full(n, "", dtype=object)
    gname = np.full(n, "", dtype=object)
    btype = np.full(n, "", dtype=object)
    for idx, name in ((anti, "antisense"), (sense, None)):    # sense overrides antisense
        m = idx >= 0
        if m.any():
            lab = labels.iloc[idx[m]]
            cat[m] = name if name else lab["category"].to_numpy()
            gene[m] = lab["gene_id"].to_numpy()
            gname[m] = lab["gene_name"].to_numpy()
            btype[m] = lab["biotype"].to_numpy()
    if args.organelles:
        org = aln["chrom"].str.fullmatch(args.organelles, case=False).to_numpy(bool)
        cat[org] = "organelle (" + aln.loc[org, "chrom"] + ")"
    aln["category"], aln["gene_id"], aln["gene_name"], aln["biotype"] = cat, gene, gname, btype

    keys = ["chrom", "strand", "pos5", "category", "gene_id", "gene_name", "biotype"]
    positions = aln.groupby(keys, sort=False).agg(reads=("read", "size"), unique_reads=("unique", "sum")).reset_index()
    positions.sort_values(["chrom", "strand", "pos5"]).to_csv(args.positions, sep="\t", index=False, compression="gzip")
    counts = aln.groupby("category").agg(reads=("read", "size"), unique_reads=("unique", "sum")).reset_index()
    counts.insert(0, "sample", args.sample)
    counts.to_csv(args.counts, sep="\t", index=False)
    print(f"{args.sample}: {n:,} genome-aligned reads", file=sys.stderr)
    print(counts.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
