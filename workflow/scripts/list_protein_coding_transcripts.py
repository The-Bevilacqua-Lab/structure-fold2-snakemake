#!/usr/bin/env python3
"""
List the protein-coding transcripts of a GTF or GFF3 that are present in the
transcriptome FASTA (one ID per line).

A transcript counts as protein-coding when its transcript-level line carries
  GTF:  transcript_biotype/transcript_type "protein_coding" (Ensembl/GENCODE),
        falling back to gene_biotype/gene_type on that line;
  GFF3: biotype/transcript_biotype/transcript_type=protein_coding, or -- when
        the line has no biotype attribute at all -- feature type mRNA.
GFF3 IDs are matched both as written and with a "transcript:" prefix removed
(Ensembl GFF3), so either form of FASTA header works.
"""

import argparse
import re
import sys

GTF_ATTR = re.compile(r'(\S+) "([^"]*)"')
TRANSCRIPT_TYPES = {"transcript", "mRNA"}
BIOTYPE_KEYS = ("transcript_biotype", "transcript_type", "biotype", "gene_biotype", "gene_type")


def parse_attributes(text):
    if "=" in text and '"' not in text.split(";")[0]:
        return dict(kv.split("=", 1) for kv in text.strip().strip(";").split(";") if "=" in kv), "gff3"
    return dict(GTF_ATTR.findall(text)), "gtf"


def protein_coding_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in TRANSCRIPT_TYPES:
                continue
            attrs, fmt = parse_attributes(fields[8])
            biotype = next((attrs[k] for k in BIOTYPE_KEYS if k in attrs), None)
            if biotype is None:
                is_pc = fields[2] == "mRNA"
            else:
                is_pc = biotype == "protein_coding"
            if not is_pc:
                continue
            tid = attrs.get("transcript_id") if fmt == "gtf" else attrs.get("ID")
            if tid:
                ids.add(tid)
                if tid.startswith("transcript:"):
                    ids.add(tid[len("transcript:"):])
    return ids


def fasta_ids(path):
    with open(path) as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotation", required=True, help="GTF or GFF3")
    parser.add_argument("--fasta", required=True, help="Transcriptome FASTA")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pc = protein_coding_ids(args.annotation)
    keep = [t for t in fasta_ids(args.fasta) if t in pc]
    if not keep:
        sys.exit(
            f"No protein-coding transcripts from {args.annotation} ({len(pc)} found there) "
            f"match the IDs in {args.fasta} -- check the annotation has biotypes and uses the same transcript IDs."
        )
    with open(args.output, "w") as out:
        out.writelines(t + "\n" for t in keep)
    print(f"{len(pc):,} protein-coding transcripts in annotation, {len(keep):,} present in the transcriptome")


if __name__ == "__main__":
    main()
