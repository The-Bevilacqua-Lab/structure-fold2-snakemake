#!/usr/bin/env python3
"""
Sum RT-stop counts per transcript from a .rtsc file, producing raw integer
counts to use as a read-count proxy for differential expression testing.

Run on the untreated (-DMS) channel: unlike the +DMS channel, its RT stops
aren't confounded by structure-dependent DMS reactivity, so their per-
transcript total tracks transcript abundance (same rationale as
StructureFold2's own rtsc_abundances.py). DESeq2 expects raw, un-normalized
counts and does its own library-size normalization, so no RPKM/TPM
conversion happens here.
"""

import argparse
import csv
from itertools import islice


def read_total_stops(rtsc_path):
    totals = {}
    with open(rtsc_path) as f:
        while True:
            record = list(islice(f, 3))
            if not record:
                break
            transcript, stops, _ = (line.strip() for line in record)
            totals[transcript] = sum(int(x) for x in stops.split("\t"))
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rtsc", required=True, help="Input .rtsc file (untreated/-DMS channel)")
    parser.add_argument("--output", required=True, help="Output CSV: transcript,count")
    args = parser.parse_args()

    totals = read_total_stops(args.rtsc)

    with open(args.output, "w", newline="") as g:
        writer = csv.writer(g)
        writer.writerow(["transcript", "count"])
        for transcript, count in totals.items():
            writer.writerow([transcript, count])


if __name__ == "__main__":
    main()
