#!/usr/bin/env python3
"""
Slice a fasta file down to a per-transcript region using 1-based inclusive
start/end coordinates from a coordinates TSV (transcript, start, end) --
see region_coordinates.py. Transcripts not listed in the coordinates file
are dropped (e.g. noncoding transcripts have no CDS/UTR region).
"""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fasta", required=True)
parser.add_argument("--coords", required=True, help="TSV: transcript, start, end (1-based, inclusive)")
parser.add_argument("--output", required=True)
args = parser.parse_args()

coords = {}
with open(args.coords) as f:
    next(f)  # header
    for line in f:
        transcript, start, end = line.rstrip("\n").split("\t")
        coords[transcript] = (int(start), int(end))

written = 0
with open(args.output, "w") as out:
    for record in SeqIO.parse(args.fasta, "fasta"):
        bounds = coords.get(record.id)
        if bounds is None:
            continue
        start, end = bounds
        region_seq = str(record.seq)[start - 1:end]
        out.write(f">{record.id}\n{region_seq}\n")
        written += 1

print(f"Wrote {written} region sequences to {args.output}")
