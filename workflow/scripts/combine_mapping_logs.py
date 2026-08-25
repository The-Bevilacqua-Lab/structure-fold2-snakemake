"""
Combines multiple per-run Bowtie2 stderr mapping logs (one per sequencing
run of the same sample -- see combine_sam_runs) into a single log in the
same textual format Bowtie2 itself produces, by summing the raw read
counts and recomputing the percentages. Written this way (rather than,
say, averaging the runs' percentages) so alignment_stats -- which parses
this exact format -- needs no changes and the combined percentages stay
correct even when runs have very different read counts.
"""

import argparse
import re


def parse_log(path):
    with open(path) as f:
        text = f.read()
    total = int(re.search(r"^(\d+) reads; of these:", text, re.M).group(1))
    unmapped = int(re.search(r"^\s*(\d+) \([\d.]+%\) aligned 0 times", text, re.M).group(1))
    unique = int(re.search(r"^\s*(\d+) \([\d.]+%\) aligned exactly 1 time", text, re.M).group(1))
    multi = int(re.search(r"^\s*(\d+) \([\d.]+%\) aligned >1 times", text, re.M).group(1))
    return total, unmapped, unique, multi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", help="Per-run Bowtie2 mapping log files to combine")
    parser.add_argument("-o", "--output", required=True, help="Combined output log path")
    args = parser.parse_args()

    total = unmapped = unique = multi = 0
    for log in args.logs:
        t, u0, u1, m = parse_log(log)
        total += t
        unmapped += u0
        unique += u1
        multi += m

    aligned = unique + multi
    overall_pct = (aligned / float(total) * 100) if total else 0.0
    unmapped_pct = (unmapped / float(total) * 100) if total else 0.0
    unique_pct = (unique / float(total) * 100) if total else 0.0
    multi_pct = (multi / float(total) * 100) if total else 0.0

    with open(args.output, "w") as g:
        g.write(f"{total} reads; of these:\n")
        g.write(f"  {total} (100.00%) were unpaired; of these:\n")
        g.write(f"    {unmapped} ({unmapped_pct:.2f}%) aligned 0 times\n")
        g.write(f"    {unique} ({unique_pct:.2f}%) aligned exactly 1 time\n")
        g.write(f"    {multi} ({multi_pct:.2f}%) aligned >1 times\n")
        g.write(f"{overall_pct:.2f}% overall alignment rate\n")


if __name__ == "__main__":
    main()
