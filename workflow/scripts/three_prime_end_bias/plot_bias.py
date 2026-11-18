#!/usr/bin/env python3
"""
3' bias correction figures (plotnine), drawn from the binned TSVs the fit/apply scripts write:

  coverage_residual_metagene.png   left: mean log(+DMS / -DMS depth), model fit and residual by bin;
                                   right: depth metagene of +DMS raw / corrected vs the -DMS reference
  rtsc_metagene_before_after.png   +DMS RT-stop metagene before (left) and after (right) correction
                                   vs the -DMS reference, with sum((mean - ref)^2) per panel

fit_3prime_bias_model.py and apply_3prime_bias_correction.py import plot_coverage / plot_rtsc.
Run as a script to redraw one ID's figures from its saved 3prime_bias_correction/ folder:

  plot_bias.py {output_dir}/{id}/3prime_bias_correction [--title "rep1 (pooled -DMS reference): "]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from plotnine import (aes, element_blank, element_text, facet_wrap, geom_hline, geom_line, geom_ribbon, ggplot,
                      guide_legend, guides, labs, scale_color_manual, scale_fill_manual, scale_linetype_manual,
                      theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


GREY_REF = "#7f7f7f"


def _base_theme():
    return (theme_bw(base_size=10)
            + theme(figure_size=(13, 4.8), panel_grid_minor=element_blank(), legend_position="bottom",
                    legend_title=element_blank(), strip_text=element_text(size=9),
                    strip_background=element_blank(), plot_title=element_text(size=11)))


def _num(model, key, default=np.nan):
    try:
        return float(model[key])
    except (KeyError, TypeError, ValueError):
        return default


def plot_coverage(binned, model, title, output):
    """binned: coverage_residual_metagene.tsv (index = bin); model: model.tsv name -> value."""
    x = binned.index.to_numpy()
    n_bins = int(x.max())
    form = str(model.get("model_form", "linear"))
    left = (f"Coverage log(+DMS / -DMS depth): raw, {form} fit and residual\n"
            f"R² = {_num(model, 'r2'):.3f}, n = {_num(model, 'n_fit_positions'):,.0f} positions")
    right = "Coverage metagene before/after correction\ndepth / transcript mean"
    series = ["raw log-ratio", f"{form} fit", "residual", "+DMS raw", "+DMS corrected", "-DMS reference"]
    rows = []
    for col, label in zip(["log_ratio", "fitted", "residual"], series[:3]):
        rows.append(pd.DataFrame({"bin": x, "value": binned[col].to_numpy(), "sem": np.nan, "series": label,
                                  "panel": left}))
    for key, label in zip(["plus_raw", "plus_corrected", "minus"], series[3:]):
        rows.append(pd.DataFrame({"bin": x, "value": binned[f"depth_{key}_mean"].to_numpy(),
                                  "sem": binned[f"depth_{key}_sem"].to_numpy(), "series": label, "panel": right}))
    df = pd.concat(rows, ignore_index=True)
    df["panel"] = pd.Categorical(df["panel"], categories=[left, right])
    df["series"] = pd.Categorical(df["series"], categories=series)
    colors = dict(zip(series, ["#dd8452", "#222222", "#55a868", "#9a9a9a", "#4c72b0", "#c44e52"]))
    linetypes = dict(zip(series, ["solid", "dashed", "solid", "solid", "solid", "solid"]))
    ref = pd.DataFrame({"panel": pd.Categorical([left, right], categories=[left, right]), "y": [0.0, 1.0]})
    band = df.dropna(subset=["sem"]).assign(lo=lambda d: d["value"] - d["sem"], hi=lambda d: d["value"] + d["sem"])
    p = (ggplot(df, aes("bin", "value", color="series", linetype="series"))
         + geom_hline(ref, aes(yintercept="y"), color=GREY_REF, linetype="dotted")
         + geom_ribbon(band, aes(x="bin", ymin="lo", ymax="hi", fill="series"), alpha=0.25, color="none",
                       inherit_aes=False)
         + geom_line(size=0.8)
         + facet_wrap("panel", scales="free_y", ncol=2)
         + scale_color_manual(values=colors) + scale_fill_manual(values=colors, guide=None)
         + scale_linetype_manual(values=linetypes)
         + guides(color=guide_legend(nrow=1), linetype=guide_legend(nrow=1))
         + labs(x=f"Bin (5′ → 3′, relative position, 1–{n_bins}; ribbons ± SEM)", y="",
                title=f"{title}3′ bias model on {_num(model, 'n_training_transcripts'):,.0f} training transcripts")
         + _base_theme())
    p.save(output, dpi=200, verbose=False)


def plot_rtsc(binned, model, title, output, n_training=None):
    """binned: rtsc_metagene_before_after.tsv (index = bin); model: model.tsv name -> value."""
    x = binned.index.to_numpy()
    n_bins = int(x.max())
    if n_training is None:
        n_training = _num(model, "n_training_transcripts")
    ref = binned["rtsc_minus_mean"].to_numpy()
    score = {k: float(np.nansum((binned[f"rtsc_{k}_mean"].to_numpy() - ref) ** 2)) for k in ("plus_raw", "plus_corrected")}
    panels = {"plus_raw": f"Before correction: raw +DMS RT-stops\nsum((mean − ref)²) = {score['plus_raw']:.3f}",
              "plus_corrected": f"After correction: corrected +DMS RT-stops\nsum((mean − ref)²) = {score['plus_corrected']:.3f}"}
    labels = {"plus": "+DMS RT-stops", "minus": "-DMS RT-stops (reference)"}
    rows = []
    for key, panel in panels.items():
        for k, label in [(key, labels["plus"]), ("minus", labels["minus"])]:
            mean, sem = binned[f"rtsc_{k}_mean"].to_numpy(), binned[f"rtsc_{k}_sem"].to_numpy()
            rows.append(pd.DataFrame({"bin": x, "value": mean, "lo": mean - sem, "hi": mean + sem,
                                      "series": label, "panel": panel}))
    df = pd.concat(rows, ignore_index=True)
    df["panel"] = pd.Categorical(df["panel"], categories=list(panels.values()))
    df["series"] = pd.Categorical(df["series"], categories=list(labels.values()))
    colors = {labels["plus"]: "#4c72b0", labels["minus"]: "#c44e52"}
    p = (ggplot(df, aes("bin", "value", color="series"))
         + geom_hline(yintercept=1.0, color=GREY_REF, linetype="dotted")
         + geom_ribbon(aes(ymin="lo", ymax="hi", fill="series"), alpha=0.3, color="none")
         + geom_line(size=0.8)
         + facet_wrap("panel", ncol=2)
         + scale_color_manual(values=colors) + scale_fill_manual(values=colors, guide=None)
         + labs(x=f"Bin (5′ → 3′, relative position, 1–{n_bins}; ribbons ± SEM)",
                y="RT-stops / transcript mean",
                title=f"{title}RT-stop 3′ bias correction, {n_training:,.0f} training transcripts "
                      f"(coverage model R² = {_num(model, 'r2'):.3f})")
         + _base_theme())
    p.save(output, dpi=200, verbose=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="An ID's 3prime_bias_correction/ folder (model.tsv + the binned TSVs)")
    ap.add_argument("--title", default=None, help="Title prefix (default: '<ID>: ')")
    args = ap.parse_args()
    d = args.folder
    title = args.title if args.title is not None else f"{os.path.basename(os.path.dirname(os.path.abspath(d)))}: "
    model = pd.read_csv(os.path.join(d, "model.tsv"), sep="\t", index_col=0, dtype={"value": str})["value"]
    cov = os.path.join(d, "coverage_residual_metagene.tsv")
    if os.path.exists(cov):
        plot_coverage(pd.read_csv(cov, sep="\t", index_col="bin"), model, title,
                      os.path.join(d, "coverage_residual_metagene.png"))
    rt = os.path.join(d, "rtsc_metagene_before_after.tsv")
    if os.path.exists(rt):
        plot_rtsc(pd.read_csv(rt, sep="\t", index_col="bin"), model, title,
                  os.path.join(d, "rtsc_metagene_before_after.png"))
    print(f"redrew the figures in {d}")


if __name__ == "__main__":
    main()
