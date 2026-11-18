#!/usr/bin/env python3
"""
Convert a StructureFold2 .rtsc file to per-nucleotide RT-stop counts.

Output is a long-format CSV with columns:
  transcript, position, nucleotide, rt_stop_count

Position is 1-based.  The same nucleotide trimming convention used by
rtsc_coverage.py is applied: the first stop value (sentinel) and the last
nucleotide of the sequence are excluded, so positions 1..N-1 are reported.
"""

import argparse
import csv
import itertools


def read_fasta(path):
    sequences = {}
    with open(path) as fh:
        current_id = None
        current_seq = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences


def read_rtsc(path):
    data = {}
    with open(path) as fh:
        while True:
            chunk = list(itertools.islice(fh, 3))
            if not chunk:
                break
            transcript, stops, _ = [l.strip() for l in chunk]
            data[transcript] = [int(x) for x in stops.split("\t")]
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Output per-nucleotide RT-stop counts from a .rtsc file"
    )
    parser.add_argument("--rtsc", required=True, help="Input .rtsc file")
    parser.add_argument("--fasta", required=True, help="Reference transcriptome FASTA")
    parser.add_argument("--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    sequences = read_fasta(args.fasta)
    rtsc_data = read_rtsc(args.rtsc)

    with open(args.output, "w", newline="") as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(["transcript", "position", "nucleotide", "rt_stop_count"])
        for transcript in sorted(rtsc_data):
            stops = rtsc_data[transcript]
            seq = sequences.get(transcript, "")
            # Convention from rtsc_coverage.py:
            #   stops[0] is a sentinel; effective stops are stops[1:]
            #   effective sequence drops the last base: seq[:-1]
            effective_stops = stops
            effective_seq = seq.upper() if seq else ""
            for i, (nuc, count) in enumerate(zip(effective_seq, effective_stops), start=1):
                writer.writerow([transcript, i, nuc, count])


if __name__ == "__main__":
    main()
