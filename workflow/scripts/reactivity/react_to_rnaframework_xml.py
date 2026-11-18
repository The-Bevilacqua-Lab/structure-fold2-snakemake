"""Write reactivities as RNA Framework XML files (the rf-norm output format), one file per transcript.

https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/

    <?xml version="1.0" encoding="UTF-8"?>
    <data tool="rf-norm" scoring="Ding" norm="2-8%" reactive="AC" ...>
        <transcript id="ID" length="N">
            <sequence>
                ACGU...
            </sequence>
            <reactivity>
                NaN,0.123,...
            </reactivity>
        </transcript>
    </data>

Input is a StructureFold <.react> (two lines per transcript: name, then tab-separated values with NA) or a
reactivity CSV (transcript,position,base,reactivity; rows of a transcript contiguous). NA -> NaN. The
<data> attributes describe how StructureFold computed the values: Ding et al. 2014 scoring
(ln(count + pseudo) +DMS minus -DMS, each scaled by its transcript mean), 2-8% normalization over the
whole transcript (win = offset = 1e9, rf-norm's default whole-transcript window), only the reactive
bases reported, values capped at max. Attributes, their order and the line wrapping (60 nt / 60 values
per line) follow rf-norm 2.9.6. The tool attribute is "rf-norm" so RNA Framework tools (rf-fold,
rf-correlate, rf-compare, ...) read the files as reactivity data.

Transcripts without a single numeric value are not written (rf-norm leaves out transcripts that do not
pass its coverage filters the same way). The sequence is the transcript FASTA's sequence, uppercase, T.
"""
import argparse
import os
import re
import sys
from xml.sax.saxutils import escape


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_react(path):
    with open(path) as f:
        while True:
            name = f.readline()
            if not name:
                return
            values = f.readline()
            yield name.strip(), values.split()


def read_csv(path):
    cur, vals = None, []
    with open(path) as f:
        header = f.readline().rstrip("\r\n").split(",")
        ti, vi = header.index("transcript"), header.index("reactivity")
        for line in f:
            col = line.rstrip("\r\n").split(",")
            if col[ti] != cur:
                if cur is not None:
                    yield cur, vals
                cur, vals = col[ti], []
            vals.append(col[vi])
    if cur is not None:
        yield cur, vals


def fmt(v, decimals):
    if v in ("NA", "NaN", "nan", ""):
        return None
    x = float(v)
    return None if x != x else f"{x:.{decimals}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--react", help="StructureFold <.react> file")
    src.add_argument("--csv", help="reactivity CSV (transcript,position,base,reactivity)")
    ap.add_argument("--fasta", required=True, help="transcript FASTA the reactivities were computed on")
    ap.add_argument("--outdir", required=True, help="output directory (one <transcript>.xml per transcript)")
    ap.add_argument("--scoring", default="Ding")
    ap.add_argument("--norm", default="2-8%", help='normalization ("2-8%%", or "raw")')
    ap.add_argument("--reactive", default="AC")
    ap.add_argument("--max", type=float, default=7.0, help="reactivity cap")
    ap.add_argument("--pseudo", type=float, default=1.0, help="pseudocount of the ln(count + pseudo) scoring")
    ap.add_argument("--win", type=int, default=1000000000)
    ap.add_argument("--decimals", type=int, default=3)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    os.makedirs(args.outdir, exist_ok=True)
    g = lambda x: f"{x:g}"
    # same attributes, in the same (alphabetical) order, as rf-norm 2.9.6 writes for Ding scoring
    attrs = {"combined": "FALSE", "max": g(args.max), "norm": args.norm, "offset": str(args.win),
             "pseudo": g(args.pseudo), "reactive": args.reactive, "remap": "0", "scoring": args.scoring,
             "tool": "rf-norm", "win": str(args.win)}
    attrs = " ".join(f'{k}="{escape(v)}"' for k, v in attrs.items())
    records = read_react(args.react) if args.react else read_csv(args.csv)
    n_written = n_empty = n_noseq = n_fixed = 0
    used = set()
    for tid, values in records:
        seq = seqs.get(tid)
        if seq is None:
            n_noseq += 1
            print(f"warning: {tid} not in the FASTA, skipped", file=sys.stderr)
            continue
        vals = [fmt(v, args.decimals) for v in values]
        if all(v is None for v in vals):
            n_empty += 1
            continue
        if len(vals) != len(seq):          # pad/trim to the sequence (the last position is always NA)
            n_fixed += 1
            print(f"warning: {tid}: {len(vals)} values for {len(seq)} nt, padded/trimmed with NaN", file=sys.stderr)
            vals = (vals + [None] * len(seq))[:len(seq)]
        seq = seq.upper().replace("U", "T")
        fname = re.sub(r"[^A-Za-z0-9._+-]", "_", tid)
        if fname in used:
            sys.exit(f"error: two transcript IDs give the file name {fname}.xml")
        used.add(fname)
        xid = escape(tid, {'"': "&quot;"})
        # line wrapping as rf-norm: 60 nt, 60 values (each followed by its comma) per line
        seq_txt = "\n\t\t\t".join(seq[k:k + 60] for k in range(0, len(seq), 60))
        r = ["NaN" if v is None else v for v in vals]
        react_txt = ",\n\t\t\t".join(",".join(r[k:k + 60]) for k in range(0, len(r), 60))
        with open(os.path.join(args.outdir, fname + ".xml"), "w") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                    f"<data {attrs}>\n"
                    f'\t<transcript id="{xid}" length="{len(seq)}">\n'
                    f"\t\t<sequence>\n\t\t\t{seq_txt}\n\t\t</sequence>\n"
                    f"\t\t<reactivity>\n\t\t\t{react_txt}\n\t\t</reactivity>\n"
                    "\t</transcript>\n</data>\n")
        n_written += 1
    print(f"{n_written} transcripts written to {args.outdir}; {n_empty} without any value left out; "
          f"{n_noseq} not in the FASTA; {n_fixed} with a length mismatch")


if __name__ == "__main__":
    main()
