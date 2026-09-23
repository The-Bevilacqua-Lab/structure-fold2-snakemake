"""
Shared pieces of the 3' coverage-bias model (fit_3prime_bias_model.py,
apply_3prime_bias_correction.py).

    log(+DMS / -DMS) = b0 + b1*d + b2*r + b3*d*r

with, for 1-based position p on a transcript of length L,
    d = L - p   (distance from the 3' end)
    r = p / L   (relative position, 5'->3')
The same coefficients are used for every transcript; the fitted multiplicative
bias at a position is exp(prediction).
"""

import numpy as np
import pandas as pd

COEF_NAMES = ["b0_intercept", "b1_distance_from_3prime", "b2_relative_position", "b3_interaction"]


def design(position, length):
    position = np.asarray(position, dtype=float)
    length = np.asarray(length, dtype=float)
    d = length - position
    r = position / length
    return np.column_stack([np.ones_like(d), d, r, d * r])


def predict_log(coefs, position, length):
    return design(position, length) @ np.asarray(coefs)


def read_model(path):
    model = pd.read_csv(path, sep="\t", index_col=0)["value"]
    return model[COEF_NAMES].to_numpy(dtype=float), model


def bin_index(position, length, n_bins):
    return np.minimum((np.asarray(position) / np.asarray(length) * n_bins).astype(int), n_bins - 1)


def metagene(transcript, position, length, values, n_bins):
    """Per transcript, value / that transcript's mean, averaged within n_bins
    relative-position bins, then across transcripts (equal weight each).
    Transcripts with mean 0 are dropped. Returns (mean, sem) per bin."""
    df = pd.DataFrame({"t": transcript, "v": np.asarray(values, dtype=float),
                       "bin": bin_index(position, length, n_bins)})
    df["v"] = df["v"] / df.groupby("t")["v"].transform("mean")
    df = df[np.isfinite(df["v"])]
    per_tx = df.groupby(["t", "bin"])["v"].mean().unstack("bin").reindex(columns=range(n_bins))
    return per_tx.mean(axis=0).to_numpy(), per_tx.sem(axis=0).to_numpy()
