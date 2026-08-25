##########################################################################
# UpSet plot of transcripts with RT-stop coverage >= threshold, one set per
# input coverage.csv (calculate_stop_coverage). Visualizes every overlap
# combination across the given replicates/IDs -- including the single bar
# for transcripts covered in every one of them, which
# filter_covered_transcripts.py reduces to a plain-text list.
##########################################################################

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_contents


def main():
    parser = argparse.ArgumentParser(
        description="UpSet plot of covered-transcript sets across replicates/IDs"
    )
    parser.add_argument(
        "--coverages", nargs="+", required=True,
        help="coverage.csv files from calculate_stop_coverage, one per replicate/ID",
    )
    parser.add_argument(
        "--labels", nargs="+", required=True,
        help="Set label for each --coverages file, in the same order",
    )
    parser.add_argument(
        "--threshold", type=float, default=1.0,
        help="Minimum mean coverage to count a transcript as covered (default: 1.0)",
    )
    parser.add_argument("--output", required=True, help="Output PNG path")
    args = parser.parse_args()

    if len(args.coverages) != len(args.labels):
        parser.error(
            f"Got {len(args.coverages)} --coverages but {len(args.labels)} --labels; "
            "need exactly one label per coverage file"
        )

    covered_sets = {}
    for path, label in zip(args.coverages, args.labels):
        df = pd.read_csv(path, header=0)
        df.columns = ["transcript", "coverage"]
        covered_sets[label] = set(df.loc[df["coverage"] >= args.threshold, "transcript"])

    data = from_contents(covered_sets)
    fig = plt.figure(figsize=(max(8, 2 + 2 * len(covered_sets)), 6))
    UpSet(data, subset_size="count", show_counts=True, sort_by="cardinality").plot(fig=fig)
    fig.suptitle(f"Transcripts with RT-stop coverage ≥ {args.threshold}")
    fig.savefig(args.output, dpi=150)


if __name__ == "__main__":
    main()
