#!/usr/bin/env python3
"""
For each transcript, find the start/end (1-based, inclusive, transcript
coordinates matching the fasta used to build transcript_position_annotations.csv)
of the requested region (5UTR, CDS, or 3UTR).

A given transcript's 5'UTR/CDS/3'UTR each occupy one unbroken stretch of a
mature mRNA, so the region's bounds are just the min/max matching position.
"""

import argparse
import csv

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--annotations", required=True, help="transcript_position_annotations.csv")
parser.add_argument("--region", required=True, choices=["5UTR", "CDS", "3UTR"])
parser.add_argument("--output", required=True, help="Output TSV: transcript, start, end")
args = parser.parse_args()

bounds = {}
with open(args.annotations) as f:
    for row in csv.DictReader(f):
        if row["region"] != args.region:
            continue
        transcript = row["transcript_id"]
        position = int(row["position"])
        if transcript not in bounds:
            bounds[transcript] = [position, position]
        else:
            entry = bounds[transcript]
            entry[0] = min(entry[0], position)
            entry[1] = max(entry[1], position)

with open(args.output, "w") as out:
    out.write("transcript\tstart\tend\n")
    for transcript, (start, end) in sorted(bounds.items()):
        out.write(f"{transcript}\t{start}\t{end}\n")

print(f"Found {args.region} in {len(bounds)} transcripts")
