#!/usr/bin/env python3
"""
Annotate each position in each transcript as 5UTR, CDS, or 3UTR.
For CDS positions, report the codon position (1, 2, or 3).

Strategy
--------
gffread -x extracts CDS sequences and -w extracts full transcript sequences
from the reference genome.  For each transcript we find the CDS as a
substring of the reference transcript (str.find), which gives us the
transcript-space CDS start and end.  Those coordinates are then applied to
the actual TRANSCRIPTOME used in this pipeline run (which may differ from
the reference at variant positions but preserves the same exon structure).

Output columns (CSV)
--------------------
transcript_id,position,region,codon_position
  position        : 1-based position in the transcript
  region          : 5UTR | CDS | 3UTR | noncoding
  codon_position  : 1, 2, or 3 for CDS bases; NA otherwise
"""

import argparse
import sys
from Bio import SeqIO


def normalize_tx_id(tx_id):
    """
    Strip a GFF3-style "transcript:" ID prefix (the Ensembl convention,
    e.g. ID=transcript:OS01T0100100-01) so transcript ids line up
    regardless of source: gffread's -x/-w FASTA headers preserve the GFF's
    raw ID attribute (prefix and all), while a plain Ensembl cdna FASTA
    (transcript_source: file) does not carry that prefix. Without this,
    every transcript would fail the cds_coords lookup and fall back to
    "noncoding".
    """
    prefix = "transcript:"
    return tx_id[len(prefix):] if tx_id.startswith(prefix) else tx_id


def load_fasta(path):
    """Return dict {normalized record.id: str(sequence).upper()}."""
    return {normalize_tx_id(rec.id): str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref_tx_fasta", required=True,
                    help="Reference transcript sequences (gffread -w output)")
    ap.add_argument("--cds_fasta",    required=True,
                    help="CDS sequences (gffread -x output)")
    ap.add_argument("--tx_fasta",     required=True,
                    help="Transcriptome used in this pipeline run (TRANSCRIPTOME)")
    ap.add_argument("--output",       required=True,
                    help="Output CSV file")
    args = ap.parse_args()

    print("Loading reference transcripts...", file=sys.stderr)
    ref_tx = load_fasta(args.ref_tx_fasta)

    print("Loading CDS sequences...", file=sys.stderr)
    cds_seqs = load_fasta(args.cds_fasta)

    print("Building CDS coordinate table...", file=sys.stderr)
    # cds_coords: tx_id -> (cds_start_0based, cds_end_0based_exclusive)
    cds_coords = {}
    missing = 0
    for tx_id, cds_seq in cds_seqs.items():
        ref_seq = ref_tx.get(tx_id)
        if ref_seq is None:
            missing += 1
            continue
        pos = ref_seq.find(cds_seq)
        if pos == -1:
            print(f"WARNING: CDS not found in reference transcript for {tx_id}",
                  file=sys.stderr)
            missing += 1
            continue
        cds_coords[tx_id] = (pos, pos + len(cds_seq))

    if missing:
        print(f"WARNING: {missing} transcripts had no CDS match; "
              "those positions are marked 'noncoding'.", file=sys.stderr)

    print("Writing per-position annotations...", file=sys.stderr)
    noncoding_total = 0
    with open(args.output, "w") as out:
        out.write("transcript_id,position,region,codon_position\n")

        for rec in SeqIO.parse(args.tx_fasta, "fasta"):
            tx_id  = normalize_tx_id(rec.id)
            tx_len = len(rec.seq)

            if tx_id not in cds_coords:
                # Non-coding or no annotation match
                noncoding_total += 1
                for pos in range(1, tx_len + 1):
                    out.write(f"{tx_id},{pos},noncoding,NA\n")
                continue

            cds_start, cds_end = cds_coords[tx_id]

            # 5'UTR
            for pos in range(1, cds_start + 1):
                out.write(f"{tx_id},{pos},5UTR,NA\n")

            # CDS with codon positions
            for i in range(cds_end - cds_start):
                pos = cds_start + i + 1
                codon_pos = (i % 3) + 1
                out.write(f"{tx_id},{pos},CDS,{codon_pos}\n")

            # 3'UTR
            for pos in range(cds_end + 1, tx_len + 1):
                out.write(f"{tx_id},{pos},3UTR,NA\n")

    print(f"Done. {noncoding_total} transcripts written as noncoding.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
