"""
Writes an uppercase-only copy of a FASTA file (headers unchanged).

Alignment is case-insensitive, but every base-identity check downstream of
mapping in this pipeline (reactivity/specificity "is this position an A or
C" tests, UTR/CDS annotation) compares against an uppercase specificity set
(e.g. "AC") case-sensitively, so soft-masked (lowercase) bases in the source
FASTA would otherwise be silently skipped -- reported as NA / excluded from
coverage and specificity counts, even though they're a perfectly good A/C.
"""

import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Source FASTA")
    parser.add_argument("--output", required=True, help="Uppercased FASTA to write")
    args = parser.parse_args()

    with open(args.output, "w") as out:
        for record in SeqIO.parse(args.input, "fasta"):
            record.seq = record.seq.upper()
            SeqIO.write(record, out, "fasta")


if __name__ == "__main__":
    main()
