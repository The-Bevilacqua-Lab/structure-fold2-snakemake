"""
Calculate reactivities from a single +DMS <.rtsc> file, with no -DMS
(background) subtraction.

Mirrors StructureFold2's rtsc_to_react.py exactly, except for the one step
it's named for: calculate_raw_reactivity() there produces
max(0, normalized_plus - normalized_minus) per position; here it's just
normalized_plus. Everything else -- the natural-log transform, the
sum-normalize-by-length step, the 2-8% percentile scale generation, and the
final threshold-capped scaling -- is the same formula, same defaults,
applied to one channel instead of two.

Output is the same two-line-per-transcript <.react> format StructureFold2
uses (name line, then tab-separated values with 'NA' at non-specificity
positions and always at the final position), so it's a drop-in input to
react_to_csv.py.
"""

import argparse
import math
import os
from itertools import islice

from Bio import SeqIO


def read_in_rtsc(rtsc_file):
    """Reads a <.rtsc> file into a dict: transcript -> [stop counts]."""
    information = {}
    with open(rtsc_file) as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript, stops, _blank = [n.strip() for n in next_n_lines]
            information[transcript] = [float(x) for x in stops.split("\t")]
    return information


def read_in_fasta(fasta_file):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}


def calculate_raw_reactivity(plus_data, ln_off=False):
    """
    Per-position: optionally ln(count + 1), then divide by that transcript's
    own sum and multiply by length -- same normalization rtsc_to_react.py
    applies to each channel before subtracting one from the other. No
    subtraction here, so no floor at 0 is needed (values are already >= 0).
    """
    data_out = {}
    for transcript, stops in plus_data.items():
        plus_vector = stops if ln_off else [math.log(v + 1, math.e) for v in stops]
        sum_plus = sum(plus_vector)
        length = len(plus_vector)
        if sum_plus != 0:
            data_out[transcript] = [v / sum_plus * length for v in plus_vector]
    return data_out


def generate_normalization_scale(derived_reactivities, transcript_seqs, specificity, trim3=0):
    """2-8% percentile scale -- identical to rtsc_to_react.py."""
    data = {}
    for transcript, reactivities in derived_reactivities.items():
        sequence = transcript_seqs[transcript]
        stop = max(1, len(reactivities) - trim3)
        accepted = sorted(
            [reactivities[k] for k in range(1, stop) if sequence[k - 1] in specificity],
            reverse=True,
        )
        top = accepted[int(len(accepted) * 0.02):int(len(accepted) * 0.1)]
        top_average = sum(top) / len(top) if top else 0
        if top_average > 0:
            data[transcript] = top_average
    return data


def read_normalization_scale(scale_file):
    info = {}
    with open(scale_file) as f:
        for line in f:
            if line.startswith("transcript"):
                continue
            transcript, value = line.strip().split(",")
            info[transcript] = float(value)
    return info


def write_normalization_scale(scale_dict, outfile):
    with open(outfile, "w") as g:
        g.write("transcript,value\n")
        for transcript, value in scale_dict.items():
            g.write(f"{transcript},{value}\n")


def calculate_final_reactivity(derived_reactivities, transcript_sequences, specificity, threshold, nrm_scale, norm_off=False):
    """Divide by scale, cap at threshold, mask non-specificity positions -- identical to rtsc_to_react.py."""
    data_out, missing_transcripts = {}, {}
    for transcript, reactivities in derived_reactivities.items():
        if transcript in nrm_scale:
            normalizer = nrm_scale[transcript] if not norm_off else 1
            sequence = transcript_sequences[transcript]
            normalized_values = [
                f"{min(reactivities[x] / normalizer, threshold):.3f}"
                if sequence[x - 1] in specificity else "NA"
                for x in range(1, len(reactivities))
            ] + ["NA"]
            data_out[transcript] = normalized_values
        else:
            missing_transcripts[transcript] = None
    return data_out, missing_transcripts


def write_out_reactivity_file(info, outfile):
    with open(outfile, "w") as g:
        for name, values in info.items():
            g.write(name + "\n")
            g.write("\t".join(values) + "\n")


def write_out_missing_file(info, outfile):
    with open(outfile, "w") as g:
        for transcript in info:
            g.write(transcript + "\n")


def check_extension(s, extension):
    return s if s.endswith(extension) else s + extension


def filter_dictionary(in_dict, filter_dict):
    for k in list(in_dict.keys()):
        if k not in filter_dict:
            del in_dict[k]


def get_covered_transcripts(coverage_file):
    with open(coverage_file) as f:
        return {line.strip(): None for line in f}


def main():
    parser = argparse.ArgumentParser(
        description="Generate a <.react> file from a single +DMS <.rtsc> file, no -DMS subtraction"
    )
    parser.add_argument("plus", type=str, help="+DMS <.rtsc> file")
    parser.add_argument("fasta", type=str, help="Transcript <.fasta> file")
    parser.add_argument("-threshold", type=float, default=7.0, help="[default = 7.0] Reactivity cap")
    parser.add_argument("-ln_off", action="store_true", help="Do not take the natural log of the stop counts")
    parser.add_argument("-nrm_off", action="store_true", help="Turn off 2-8%% normalization of the derived reactivity")
    parser.add_argument("-save_fails", action="store_true", help="Log transcripts with zero or missing scales")
    parser.add_argument("-scale", type=str, default=None, help="Provide a normalization <.scale> for calculation")
    parser.add_argument("-bases", type=str, default="AC", help="[default = AC] Reaction specificity")
    parser.add_argument("-name", type=str, default=None, help="Change the name of the outfile")
    parser.add_argument("-restrict", default=None, help="Limit analysis to these specific transcripts <.txt>")
    parser.add_argument("-trim3", type=int, default=0, help="[default = 0] Bases to exclude from the 3' end before normalization")
    args = parser.parse_args()

    base_name = os.path.basename(args.plus).replace(".rtsc", "")
    out_name = check_extension(args.name, ".react") if args.name else f"{base_name}_plus_only.react"

    plus_data = read_in_rtsc(args.plus)
    if args.restrict is not None:
        covered = get_covered_transcripts(args.restrict)
        filter_dictionary(plus_data, covered)

    data = calculate_raw_reactivity(plus_data, args.ln_off)

    seqs = read_in_fasta(args.fasta)

    normalization_scale = (
        generate_normalization_scale(data, seqs, args.bases, args.trim3)
        if args.scale is None else read_normalization_scale(args.scale)
    )
    if args.scale is None:
        write_normalization_scale(normalization_scale, out_name.replace(".react", ".scale"))

    out_reactivity, out_missing = calculate_final_reactivity(
        data, seqs, args.bases, args.threshold, normalization_scale, args.nrm_off
    )

    write_out_reactivity_file(out_reactivity, out_name)
    if args.save_fails:
        write_out_missing_file(out_missing, out_name.replace(".react", "_unresolvable_transcripts.txt"))


if __name__ == "__main__":
    main()
