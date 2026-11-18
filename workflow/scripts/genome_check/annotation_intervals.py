#!/usr/bin/env python3
"""
Genome annotation (GTF or GFF3) -> the exon and gene intervals used to place genome-aligned reads.

Writes a gzip TSV, one row per interval (1-based, inclusive):
  chrom start end strand feature gene_id gene_name transcript_id biotype in_transcriptome
feature is "exon" (one row per exon of every transcript) or "gene" (the span of a gene: its gene
feature, else the span of its transcripts). biotype is the transcript's biotype (transcript_biotype
/ biotype / the GFF3 feature type, e.g. tRNA), else its gene's. in_transcriptome = 1 for exons of
transcripts whose ID is a header of --transcriptome (after removing --prefix and any
"transcript:" prefix), i.e. the transcripts the reads were aligned to.

GFF3 transcripts without exon children use their CDS/UTR pieces as exons (as gffread does).
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

import pandas as pd

TX_FEATURES_NOT = {"gene", "exon", "CDS", "five_prime_UTR", "three_prime_UTR", "five_prime_utr", "three_prime_utr",
                   "start_codon", "stop_codon", "UTR", "chromosome", "region", "biological_region", "intron",
                   "pseudogene", "ncRNA_gene", "Selenocysteine", "scaffold", "contig", "supercontig"}
GENE_FEATURES = {"gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"}
PIECE_FEATURES = {"CDS", "five_prime_UTR", "three_prime_UTR", "five_prime_utr", "three_prime_utr", "UTR"}


def strip_id(x):
    return re.sub(r"^(transcript|gene|rna|mrna):", "", x, flags=re.IGNORECASE) if x else x


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_attrs(field, gtf):
    if gtf:
        return {k: v for k, v in re.findall(r'(\S+)\s+"([^"]*)"', field)}
    out = {}
    for part in field.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def read_annotation(path):
    gtf = path.lower().rstrip(".gz").endswith((".gtf", ".gtf2"))
    genes = {}                         # gene_id -> dict(chrom, start, end, strand, name, biotype)
    txs = {}                           # transcript_id -> dict(gene, biotype, chrom, strand)
    exons = defaultdict(list)          # transcript_id -> [(start, end)]
    pieces = defaultdict(list)         # GFF3: transcript_id -> CDS/UTR pieces
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            chrom, ftype, start, end, strand = c[0], c[2], int(c[3]), int(c[4]), c[6]
            a = parse_attrs(c[8], gtf)
            if gtf:
                gid, tid = strip_id(a.get("gene_id")), strip_id(a.get("transcript_id"))
                if ftype == "gene" and gid:
                    genes[gid] = dict(chrom=chrom, start=start, end=end, strand=strand,
                                      name=a.get("gene_name", ""), biotype=a.get("gene_biotype", a.get("gene_type", "")))
                if tid:
                    t = txs.setdefault(tid, dict(gene=gid, biotype="", chrom=chrom, strand=strand))
                    t["biotype"] = t["biotype"] or a.get("transcript_biotype", a.get("transcript_type", ""))
                    if gid and gid not in genes:
                        genes[gid] = dict(chrom=chrom, start=start, end=end, strand=strand,
                                          name=a.get("gene_name", ""), biotype=a.get("gene_biotype", a.get("gene_type", "")),
                                          from_tx=True)
                    elif gid and genes[gid].get("from_tx"):
                        genes[gid]["start"] = min(genes[gid]["start"], start)
                        genes[gid]["end"] = max(genes[gid]["end"], end)
                    if ftype == "exon":
                        exons[tid].append((start, end))
                continue
            # GFF3
            fid = strip_id(a.get("ID", ""))
            parents = [strip_id(p) for p in a.get("Parent", "").split(",") if p]
            if ftype in GENE_FEATURES:
                genes[fid] = dict(chrom=chrom, start=start, end=end, strand=strand, name=a.get("Name", ""),
                                  biotype=a.get("biotype", a.get("gene_biotype", "pseudogene" if ftype == "pseudogene" else "")))
            elif ftype == "exon":
                for p in parents:
                    exons[p].append((start, end))
            elif ftype in PIECE_FEATURES:
                for p in parents:
                    pieces[p].append((start, end))
            elif ftype not in TX_FEATURES_NOT and fid and "Parent" in a:
                tid = strip_id(a.get("transcript_id", fid))
                bt = a.get("biotype", a.get("transcript_biotype", ""))
                if not bt and ftype not in ("mRNA", "transcript"):
                    bt = ftype
                elif not bt and ftype == "mRNA":
                    bt = "protein_coding"
                txs[fid] = dict(gene=parents[0] if parents else "", biotype=bt, chrom=chrom, strand=strand, tid=tid,
                                start=start, end=end)
    if not gtf:
        # GFF3: exons hang off the feature ID; report the transcript_id attribute when given
        for fid, t in txs.items():
            if not exons.get(fid) and pieces.get(fid):
                exons[fid] = merge(pieces[fid])
            if not exons.get(fid):
                exons[fid] = [(t["start"], t["end"])]
        exons = {txs[fid].get("tid", fid): v for fid, v in exons.items() if fid in txs}
        txs = {t.get("tid", fid): t for fid, t in txs.items()}
        for t in txs.values():   # gene span from transcripts when the gene line is missing
            g = t["gene"]
            if g and g not in genes:
                genes[g] = dict(chrom=t["chrom"], start=t["start"], end=t["end"], strand=t["strand"], name="", biotype="")
    return genes, txs, exons


def merge(iv):
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation", default="", help="GTF or GFF3 (empty: write an empty table)")
    ap.add_argument("--transcriptome", required=True, help="FASTA the reads were aligned to")
    ap.add_argument("--prefix", default="", help="Remove this prefix from transcriptome IDs (genotype mode)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    cols = ["chrom", "start", "end", "strand", "feature", "gene_id", "gene_name", "transcript_id", "biotype",
            "in_transcriptome"]
    if not args.annotation:
        pd.DataFrame(columns=cols).to_csv(args.output, sep="\t", index=False, compression="gzip")
        print("no annotation: wrote an empty interval table", file=sys.stderr)
        return

    in_tx = set()
    with open_text(args.transcriptome) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if args.prefix:
                    if not name.startswith(args.prefix):
                        continue
                    name = name[len(args.prefix):]
                in_tx.add(strip_id(name))

    genes, txs, exons = read_annotation(args.annotation)
    rows = []
    for tid, ivs in exons.items():
        t = txs.get(tid)
        if t is None:
            continue
        g = genes.get(t["gene"], {})
        bt = t["biotype"] or g.get("biotype", "") or "unknown"
        flag = int(tid in in_tx)
        for s, e in ivs:
            rows.append((t["chrom"], s, e, t["strand"], "exon", t["gene"], g.get("name", ""), tid, bt, flag))
    for gid, g in genes.items():
        rows.append((g["chrom"], g["start"], g["end"], g["strand"], "gene", gid, g.get("name", ""), "",
                     g.get("biotype", "") or "unknown", 0))
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(args.output, sep="\t", index=False, compression="gzip")
    n_hit = len(in_tx & set(exons))
    print(f"{len(genes):,} genes, {len(exons):,} transcripts, {(df['feature'] == 'exon').sum():,} exons; "
          f"{n_hit:,} of {len(in_tx):,} transcriptome IDs found in the annotation", file=sys.stderr)
    if in_tx and n_hit == 0:
        print("WARNING: no transcriptome ID matches an annotation transcript ID; every hit will be "
              "labelled 'not in transcriptome'", file=sys.stderr)


if __name__ == "__main__":
    main()
