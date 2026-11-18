##########################################################################
# Write a plain-text list of transcripts that have mean RT-stop coverage
# >= threshold in ALL supplied coverage CSV files (intersection).
##########################################################################

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Intersect covered transcripts across replicates"
)
parser.add_argument("--coverages", nargs="+", required=True,
                    help="coverage.csv files from calculate_stop_coverage, one per replicate")
parser.add_argument("--threshold", type=float, default=1.0,
                    help="Minimum mean coverage to keep a transcript (default: 1.0)")
parser.add_argument("--output", required=True,
                    help="Output plain-text file: one transcript name per line")
args = parser.parse_args()

covered_sets = []
for path in args.coverages:
    df = pd.read_csv(path, header=0)
    df.columns = ["transcript", "coverage"]
    covered_sets.append(set(df.loc[df["coverage"] >= args.threshold, "transcript"]))

common = set.intersection(*covered_sets)

with open(args.output, "w") as fh:
    for transcript in sorted(common):
        fh.write(transcript + "\n")

print(f"Wrote {len(common)} transcripts to {args.output}")
