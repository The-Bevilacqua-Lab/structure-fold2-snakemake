#!/usr/bin/env python3
"""Concatenate the shards of depth_matched_reactivity.py (--shard/--n-shards), in shard order:
per library reactivity.csv (one header) and matched_plus_draw1.rtsc, and summary.tsv (one header).
The result is identical to an unsharded run."""
import argparse
import os
import shutil

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--shard-dirs", nargs="+", required=True, help="shard output dirs, in shard order; each has <label>/ and summary.tsv")
ap.add_argument("--labels", nargs="+", required=True)
ap.add_argument("--outdirs", nargs="+", required=True, help="final output dir per label")
ap.add_argument("--summary", required=True)
args = ap.parse_args()


def cat(parts, out, header):
    with open(out, "w") as o:
        for k, part in enumerate(parts):
            with open(part) as f:
                if header:
                    first = f.readline()
                    if k == 0:
                        o.write(first)
                shutil.copyfileobj(f, o, 1 << 24)


for label, outdir in zip(args.labels, args.outdirs):
    os.makedirs(outdir, exist_ok=True)
    cat([os.path.join(d, label, "reactivity.csv") for d in args.shard_dirs], os.path.join(outdir, "reactivity.csv"), True)
    cat([os.path.join(d, label, "matched_plus_draw1.rtsc") for d in args.shard_dirs],
        os.path.join(outdir, "matched_plus_draw1.rtsc"), False)
cat([os.path.join(d, "summary.tsv") for d in args.shard_dirs], args.summary, True)
n = sum(1 for _ in open(args.summary)) - 1
print(f"merged {len(args.shard_dirs)} shards: {n:,} transcripts in {args.summary}")
