#!/usr/bin/env python3
"""
NA out the 3'-trimmed tail of every transcript in a <.react> file.

rtsc_to_react.py's -trim3 (and rtsc_to_react_plus_only.py's copy of the same
logic) already excludes the last `trim3` nucleotide positions from 2-8%
normalization-SCALE generation (see generate_normalization_scale's
`stop = max(1, len(reactivities) - trim3)`), but it does NOT exclude those
same positions from the per-position reactivity VALUES it reports -- they
still get a real computed number in the .react file. Since heat correction
(react_heat_correct.py's sum_react) sums every non-NA value across the whole
transcript to compute its scaling factor, those trim3'd positions were being
folded into that sum even though trim3 says to disregard them, and they show
up as ordinary (non-NA) values in reactivity.csv too.

This script re-derives the exact same `stop` cutoff and masks everything at
or past it to 'NA', so a value trim3 already excluded from the scale is also
excluded from the reported reactivity, from heat correction's scaling-factor
calculation (apply_correction leaves 'NA' as 'NA', and its own sum already
skips non-float entries), and from anything reading the resulting .csv.

Run on the .react output of any rtsc_to_react.py/rtsc_to_react_plus_only.py
invocation whose *scale* was (or -- for a shared/externally-supplied scale,
would have been, had it generated its own -- see rtsc_to_react_heat_shared in
workflow/rules/structurefold2.smk) generated with this trim3 value. trim3=0
is a no-op: every .react file already carries a structural trailing 'NA' at
its very last position regardless of trim3 (see calculate_final_reactivity),
so masking with trim3=0 reproduces the input unchanged.
"""

import argparse
from itertools import islice


def read_react(path):
    """Return {transcript: [tab-separated value strings]} -- values are left
    as raw strings ('NA' or a numeric literal) since this script only ever
    rewrites the tail to 'NA' and never needs to parse the rest."""
    info = {}
    with open(path) as f:
        while True:
            lines = list(islice(f, 2))
            if not lines:
                break
            transcript, values = (line.strip() for line in lines)
            info[transcript] = values.split("\t")
    return info


def mask_trim3(values, trim3):
    """
    Mask the trailing (trim3 + 1) entries to 'NA': the same `stop` cutoff
    generate_normalization_scale uses (see module docstring) leaves
    `stop - 1` leading positions in the scale, where
    stop = max(1, len(values) - trim3) -- so everything from index
    `stop - 1` onward (inclusive) is what trim3 already excludes from the
    scale, and now also gets excluded from the reported value. trim3 <= 0
    returns the input unchanged.
    """
    if trim3 <= 0:
        return values
    stop = max(1, len(values) - trim3)
    keep = stop - 1
    return values[:keep] + ["NA"] * (len(values) - keep)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input <.react> file")
    parser.add_argument("--output", required=True, help="Output <.react> file")
    parser.add_argument(
        "--trim3", type=int, required=True,
        help="Bases trimmed from the 3' end when the scale was generated (0 = no-op)",
    )
    args = parser.parse_args()

    data = read_react(args.input)
    with open(args.output, "w") as g:
        for transcript, values in data.items():
            g.write(transcript + "\n")
            g.write("\t".join(mask_trim3(values, args.trim3)) + "\n")


if __name__ == "__main__":
    main()
