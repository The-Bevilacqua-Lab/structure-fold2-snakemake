#!/usr/bin/env python3
"""
Fit the 3' coverage-bias model (see bias_model.py) on per-nucleotide read
depth of the training transcripts: log(+DMS depth / -DMS reference depth)
regressed on distance from the 3' end, relative position, and their
interaction, using every position with non-zero depth in both.

Outputs
  --model   TSV (name, value): the 4 coefficients, r2, n_fit_positions,
            n_training_transcripts
  --binned  TSV per relative-position bin: mean raw log-ratio, fit and
            residual, plus the depth metagenes (+DMS raw / corrected, -DMS)
  --plot    PNG: left, raw log-ratio / fit / residual by bin; right, depth
            metagene of +DMS before and after correction vs the -DMS reference
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bias_model import COEF_NAMES, bin_index, design, metagene


def read_depth(path, name):
    return pd.read_csv(path, sep="\t", header=None, names=["transcript", "position", name],
                       dtype={"transcript": str, "position": np.int64, name: float})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--depth-plus", required=True, help="transcript, position, depth (+DMS)")
    parser.add_argument("--depth-minus", required=True, help="transcript, position, depth (-DMS reference)")
    parser.add_argument("--training", required=True, help="Training transcripts TSV (transcript, length, ...)")
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
    if len(fit) < len(COEF_NAMES):
        raise SystemExit(f"Only {len(fit)} positions have non-zero +DMS and -DMS depth -- cannot fit the model")
    y = np.log(fit["plus"].to_numpy() / fit["minus"].to_numpy())
    X = design(fit["position"], fit["length"])
    coefs, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ coefs
    r2 = 1 - np.sum((y - y_hat) ** 2) / np.sum((y - y.mean()) ** 2)

    model = pd.Series(
        list(coefs) + [r2, len(fit), df["transcript"].nunique()],
        index=COEF_NAMES + ["r2", "n_fit_positions", "n_training_transcripts"],
        name="value",
    )
    model.index.name = "name"
    model.to_csv(args.model, sep="\t")
    print(model.to_string())

    # Log-ratio / fit / residual by bin (fit positions only -- the log-ratio
    # is undefined elsewhere).
    b = bin_index(fit["position"], fit["length"], args.bins)
    by_bin = pd.DataFrame({"bin": b, "log_ratio": y, "fitted": y_hat, "residual": y - y_hat}).groupby("bin").mean()
    by_bin = by_bin.reindex(range(args.bins))

    # Depth metagene before/after, over every training position.
    corrected = df["plus"].to_numpy() / np.exp(design(df["position"], df["length"]) @ coefs)
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

    x = np.arange(1, args.bins + 1)
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))
    ax = axes[0]
    ax.plot(x, binned["log_ratio"], color="tab:orange", label="raw log(+DMS / -DMS depth)")
    ax.plot(x, binned["fitted"], color="tab:blue", linestyle="--", label="regression fit")
    ax.plot(x, binned["residual"], color="tab:green", label="residual")
    ax.axhline(0.0, color="gray", linestyle=":", linewidth=1)
    ax.set_ylabel("log(+DMS depth / -DMS depth)")
    ax.set_title(f"Coverage log-ratio, fit and residual (R$^2$ = {r2:.3f}, n = {len(fit):,})")

    ax = axes[1]
    for key, color, label in [("plus_raw", "tab:gray", "+DMS raw"),
                              ("plus_corrected", "tab:blue", "+DMS corrected"),
                              ("minus", "tab:red", "-DMS reference")]:
        mean, sem = meta[key]
        ax.plot(x, mean, color=color, label=label)
        ax.fill_between(x, mean - sem, mean + sem, color=color, alpha=0.25)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_ylabel("Depth / transcript-mean depth")
    ax.set_title("Coverage metagene before/after correction")

    for ax in axes:
        ax.set_xlabel(f"Bin (5' -> 3', relative position, 1-{args.bins})")
        ax.legend()
    fig.suptitle(f"{args.title}3' bias model on {model['n_training_transcripts']:,.0f} training transcripts")
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)


if __name__ == "__main__":
    main()
