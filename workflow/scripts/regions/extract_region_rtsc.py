#!/usr/bin/env python3
"""
Slice a <.rtsc> file's per-transcript RT-stop-count vectors down to the same
per-transcript region used by extract_region_fasta.py, so the sliced counts
line up 1:1 (same length, same base at each index) with the region-sliced
transcriptome fed to rtsc_to_react.py. Transcripts not listed in the
coordinates file are dropped.
"""

import argparse
from itertools import islice

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rtsc", required=True)
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
with open(args.rtsc) as f, open(args.output, "w") as out:
    while True:
        chunk = list(islice(f, 3))
        if not chunk:
            break
        transcript, counts_line, _blank = [n.rstrip("\n") for n in chunk]
        bounds = coords.get(transcript)
        if bounds is None:
            continue
        start, end = bounds
        region_counts = counts_line.split("\t")[start - 1:end]
        out.write(transcript + "\n")
        out.write("\t".join(region_counts) + "\n\n")
        written += 1

print(f"Wrote {written} region RT-stop vectors to {args.output}")
