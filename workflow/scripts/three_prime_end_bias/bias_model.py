"""
Shared pieces of the 3' coverage-bias model (fit_3prime_bias_model.py,
apply_3prime_bias_correction.py).

For 1-based position p on a transcript of length L,
    d = L - p   (distance from the 3' end)
    r = p / L   (relative position, 5'->3')
the model is fit to log(+DMS depth / -DMS depth). Two forms:

  linear  log(+/-) = b0 + b1*d + b2*r + b3*d*r
  spline  log(+/-) = b0 + sum_j beta_j N_j(u) + b_r*r + sum_m gamma_m M_m(v)
          u = min(d, D_CAP) / D_CAP, v = min(p, P_CAP) / P_CAP,
          N_j, M_m natural cubic spline bases (Harrell's restricted cubic
          spline parameterisation) with knots SPLINE_KNOTS_3P (on d) and
          SPLINE_KNOTS_5P (on p). The 3' spline carries the bias shape and is
          flat beyond D_CAP nt from the 3' end; the 5' spline absorbs the 5'-end
          depletion and is flat beyond P_CAP nt. See docs/3_prime_end_bias.md.

The same coefficients are used for every transcript; the fitted multiplicative
bias at a position is exp(prediction). model.tsv records the form (and, for
the spline, its knots and caps), so read_model() rebuilds the same basis; a
model.tsv without a model_form row is a linear model.
"""

import numpy as np
import pandas as pd

COEF_NAMES = ["b0_intercept", "b1_distance_from_3prime", "b2_relative_position", "b3_interaction"]
MODEL_FORMS = ("linear", "spline")

D_CAP = 2500.0
P_CAP = 300.0
SPLINE_KNOTS_3P = (0.0,) + tuple(np.geomspace(25, D_CAP, 7))   # nt of d: 0, 25, 53.9, ..., 2500
SPLINE_KNOTS_5P = (0.0, 10.0, 30.0, 100.0, 300.0)               # nt of p


def design(position, length):
    """Linear-form design matrix: 1, d, r, d*r."""
    position = np.asarray(position, dtype=float)
    length = np.asarray(length, dtype=float)
    d = length - position
    r = position / length
    return np.column_stack([np.ones_like(d), d, r, d * r])


def natural_spline_basis(x, knots):
    """Natural cubic spline basis on x for knots t_1 < ... < t_K: the linear term
    x / (t_K - t_1) plus K - 2 cubic terms; linear beyond the boundary knots."""
    x = np.asarray(x, dtype=float)
    t = np.asarray(knots, dtype=float)
    pos3 = lambda v: np.clip(v, 0, None) ** 3
    scale = (t[-1] - t[0]) ** 2
    cols = [x / (t[-1] - t[0])]
    for j in range(len(t) - 2):
        cols.append((pos3(x - t[j])
                     - pos3(x - t[-2]) * (t[-1] - t[j]) / (t[-1] - t[-2])
                     + pos3(x - t[-1]) * (t[-2] - t[j]) / (t[-1] - t[-2])) / scale)
    return np.column_stack(cols)


def spline_design(position, length, knots_3p=SPLINE_KNOTS_3P, knots_5p=SPLINE_KNOTS_5P,
                  d_cap=D_CAP, p_cap=P_CAP):
    """Spline-form design matrix: 1, N(u) (len(knots_3p) - 1 columns), r, M(v) (len(knots_5p) - 1)."""
    position = np.asarray(position, dtype=float)
    length = np.asarray(length, dtype=float)
    u = np.minimum(length - position, d_cap) / d_cap
    v = np.minimum(position, p_cap) / p_cap
    return np.column_stack([
        np.ones_like(u),
        natural_spline_basis(u, np.asarray(knots_3p, dtype=float) / d_cap),
        position / length,
        natural_spline_basis(v, np.asarray(knots_5p, dtype=float) / p_cap),
    ])


def spline_coef_names(knots_3p=SPLINE_KNOTS_3P, knots_5p=SPLINE_KNOTS_5P):
    return (["b0_intercept"]
            + [f"beta_3prime_{j}" for j in range(1, len(knots_3p))]
            + ["b_relative_position"]
            + [f"gamma_5prime_{m}" for m in range(1, len(knots_5p))])


class BiasModel:
    """A fitted model: its form, coefficients and (spline) basis settings."""

    def __init__(self, form, coefs, knots_3p=SPLINE_KNOTS_3P, knots_5p=SPLINE_KNOTS_5P,
                 d_cap=D_CAP, p_cap=P_CAP):
        if form not in MODEL_FORMS:
            raise ValueError(f"model form must be one of {MODEL_FORMS}, got {form!r}")
        self.form = form
        self.coefs = np.asarray(coefs, dtype=float)
        self.knots_3p, self.knots_5p = tuple(knots_3p), tuple(knots_5p)
        self.d_cap, self.p_cap = float(d_cap), float(p_cap)

    def design(self, position, length):
        if self.form == "linear":
            return design(position, length)
        return spline_design(position, length, self.knots_3p, self.knots_5p, self.d_cap, self.p_cap)

    def predict_log(self, position, length):
        return self.design(position, length) @ self.coefs

    def coef_names(self):
        return COEF_NAMES if self.form == "linear" else spline_coef_names(self.knots_3p, self.knots_5p)

    def to_series(self):
        """name -> value rows for model.tsv (coefficients plus the settings needed to rebuild the basis)."""
        rows = {"model_form": self.form}
        if self.form == "spline":
            rows.update({
                "d_cap": self.d_cap,
                "p_cap": self.p_cap,
                "knots_3prime_nt": ",".join(repr(float(k)) for k in self.knots_3p),
                "knots_5prime_nt": ",".join(repr(float(k)) for k in self.knots_5p),
            })
        rows.update(dict(zip(self.coef_names(), self.coefs)))
        return pd.Series(rows, name="value")


def empty_model(form):
    """An unfitted model of the given form (used to build the design before fitting)."""
    n = len(COEF_NAMES) if form == "linear" else len(spline_coef_names())
    return BiasModel(form, np.zeros(n))


def predict_log(model, position, length):
    """Fitted log bias. model is a BiasModel, or a bare linear coefficient vector (legacy)."""
    if isinstance(model, BiasModel):
        return model.predict_log(position, length)
    return design(position, length) @ np.asarray(model)


def _number_or_text(value):
    try:
        return float(value)
    except ValueError:
        return value


def read_model(path):
    """Returns (BiasModel, the raw name -> value Series). A model.tsv without a
    model_form row is a linear model (files written before the spline form)."""
    raw = pd.read_csv(path, sep="\t", index_col=0, dtype={"value": str})["value"]
    raw = raw.map(_number_or_text)
    form = raw.get("model_form", "linear")
    if form == "linear":
        model = BiasModel("linear", raw[COEF_NAMES].astype(float).to_numpy())
    else:
        knots_3p = tuple(float(k) for k in raw["knots_3prime_nt"].split(","))
        knots_5p = tuple(float(k) for k in raw["knots_5prime_nt"].split(","))
        names = spline_coef_names(knots_3p, knots_5p)
        model = BiasModel(form, raw[names].astype(float).to_numpy(), knots_3p, knots_5p,
                          float(raw["d_cap"]), float(raw["p_cap"]))
    return model, raw


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
