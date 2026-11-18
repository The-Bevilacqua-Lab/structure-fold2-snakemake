#!/usr/bin/env python3
"""
Fit the 3' coverage-bias model (see bias_model.py) on per-nucleotide read
depth of the training transcripts: log(+DMS depth / -DMS reference depth),
using every position with non-zero depth in both. --model-form picks the form:
spline (default; natural cubic spline in distance from the 3' end + relative
position + natural cubic spline in distance from the 5' end) or linear
(distance from the 3' end, relative position and their interaction).

Outputs
  --model   TSV (name, value): model_form, the spline knots/caps (spline form),
            the coefficients, r2, n_fit_positions, n_training_transcripts
  --binned  TSV per relative-position bin: mean raw log-ratio, fit and
            residual, plus the depth metagenes (+DMS raw / corrected, -DMS)
  --plot    PNG (plot_bias.py): left, raw log-ratio / fit / residual by bin; right, depth
            metagene of +DMS before and after correction vs the -DMS reference
"""

import argparse

import numpy as np
import pandas as pd

from bias_model import MODEL_FORMS, BiasModel, bin_index, empty_model, metagene
from plot_bias import plot_coverage


def read_depth(path, name):
    return pd.read_csv(path, sep="\t", header=None, names=["transcript", "position", name],
                       dtype={"transcript": str, "position": np.int64, name: float})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--depth-plus", required=True, help="transcript, position, depth (+DMS)")
    parser.add_argument("--depth-minus", required=True, help="transcript, position, depth (-DMS reference)")
    parser.add_argument("--training", required=True, help="Training transcripts TSV (transcript, length, ...)")
    parser.add_argument("--model-form", choices=MODEL_FORMS, default="spline")
    parser.add_argument("--bins", type=int, default=50)
    parser.add_argument("--title", default="")
    parser.add_argument("--model", required=True)
    parser.add_argument("--binned", required=True)
    parser.add_argument("--plot", required=True)
    args = parser.parse_args()

    lengths = pd.read_csv(args.training, sep="\t", dtype={"transcript": str}).set_index("transcript")["length"]
    df = read_depth(args.depth_plus, "plus").merge(read_depth(args.depth_minus, "minus"),
                                                   on=["transcript", "position"], how="inner")
    df["length"] = df["transcript"].map(lengths)
    df = df.dropna(subset=["length"])
    df = df[(df["position"] >= 1) & (df["position"] <= df["length"])]

    fit = df[(df["plus"] > 0) & (df["minus"] > 0)]
    X = empty_model(args.model_form).design(fit["position"], fit["length"])
    if len(fit) < X.shape[1]:
        raise SystemExit(f"Only {len(fit)} positions have non-zero +DMS and -DMS depth -- cannot fit the model")
    y = np.log(fit["plus"].to_numpy() / fit["minus"].to_numpy())
    if args.model_form == "linear":
        coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
    else:
        # The spline columns span very different scales; solve on max-abs-scaled columns.
        col_scale = np.abs(X).max(axis=0)
        col_scale[col_scale == 0] = 1.0
        coefs, *_ = np.linalg.lstsq(X / col_scale, y, rcond=None)
        coefs = coefs / col_scale
    fitted_model = BiasModel(args.model_form, coefs)
    y_hat = X @ coefs
    r2 = 1 - np.sum((y - y_hat) ** 2) / np.sum((y - y.mean()) ** 2)

    model = pd.concat([
        fitted_model.to_series(),
        pd.Series({"r2": r2, "n_fit_positions": len(fit), "n_training_transcripts": df["transcript"].nunique()}),
    ]).rename("value")
    model.index.name = "name"
    model.to_csv(args.model, sep="\t")
    print(model.to_string())

    # Log-ratio / fit / residual by bin (fit positions only -- the log-ratio
    # is undefined elsewhere).
    b = bin_index(fit["position"], fit["length"], args.bins)
    by_bin = pd.DataFrame({"bin": b, "log_ratio": y, "fitted": y_hat, "residual": y - y_hat}).groupby("bin").mean()
    by_bin = by_bin.reindex(range(args.bins))

    # Depth metagene before/after, over every training position.
    corrected = df["plus"].to_numpy() / np.exp(fitted_model.predict_log(df["position"], df["length"]))
    meta = {}
    for key, values in [("plus_raw", df["plus"]), ("plus_corrected", corrected), ("minus", df["minus"])]:
        meta[key] = metagene(df["transcript"].to_numpy(), df["position"], df["length"], values, args.bins)

    binned = by_bin.copy()
    for key, (mean, sem) in meta.items():
        binned[f"depth_{key}_mean"] = mean
        binned[f"depth_{key}_sem"] = sem
    binned.index = binned.index + 1
    binned.index.name = "bin"
    binned.to_csv(args.binned, sep="\t")

    plot_coverage(binned, model, args.title, args.plot)


if __name__ == "__main__":
    main()
