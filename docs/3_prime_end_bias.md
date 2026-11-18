# 3′ coverage-bias correction

Optional (`3_prime_bias_correction: true`). This corrects the +DMS RT stops for a coverage
gradient along transcripts that the −DMS library doesn't share, before reactivity is calculated.
Rules and files: [pipeline.md §8](pipeline.md#8-optional-3-coverage-bias-correction). The
reactivity calculation it feeds: [methods.md §8](methods.md#8-reactivity).

## Contents

1. [The problem](#1-the-problem)
2. [Overview](#2-overview)
3. [Training transcripts](#3-training-transcripts)
4. [Read depth](#4-read-depth)
5. [The model](#5-the-model)
6. [Fitting](#6-fitting)
7. [Applying the correction](#7-applying-the-correction)
8. [Outputs and how to read them](#8-outputs-and-how-to-read-them)
9. [Assumptions and limits](#9-assumptions-and-limits)
10. [Settings](#10-settings)

---

## 1. The problem

With poly(A) selection or oligo-dT priming, reads pile up toward the 3′ end of mRNAs. If the
+DMS and −DMS libraries had the same gradient, the reactivity calculation would cancel it, since
each channel is scaled to its own transcript mean and then subtracted. They often don't: the
+DMS library goes through DMS treatment and quenching, and RNA degradation or fragmentation
there can make its gradient steeper or shallower than the −DMS one. Without correction:

- the 3′ end looks more (or less) reactive than the 5′ end, for reasons unrelated to structure;
- the 2–8% normalizer is drawn mostly from one end of the transcript;
- reactivities of the same base are less comparable between transcripts of different length.

The correction removes the **difference** in gradient between +DMS and −DMS, measured on read
depth, from the +DMS RT stops.

## 2. Overview

Per reactivity ID:

```
+DMS .rtsc of every ID ─► pick training transcripts (one set shared by all IDs)
                                     │
+DMS BAMs of the ID  ─┐              ▼
                      ├─► samtools depth on the training transcripts
−DMS reference BAMs  ─┘              │
                                     ▼
          fit  log(+DMS depth / −DMS depth)  ~  f(position, length)    → model.tsv
                                     │
   divide every transcript's +DMS stops by exp(f), rescale to the raw total
                                     │
                                     ▼
            combined_plus_corrected.rtsc  →  every reactivity downstream
```

## 3. Training transcripts

The model is fitted on protein-coding transcripts only. They are the transcripts that went through
poly(A) selection, and they're long enough to show the gradient. Protein-coding transcripts come
from `3_prime_bias_annotation` (default `annotation_gtf`, then `gff`): GTF `transcript_biotype`
or `transcript_type` `protein_coding` (falling back to the gene's), GFF3 `biotype=protein_coding`,
or GFF3 `mRNA` features with no biotype.

Each candidate is scored by its **lowest** +DMS coverage across all IDs (coverage = A/C stops /
A/C bases, as in `coverage.csv`). Then:

- `top_n` (default): the top `3_prime_bias_n_transcripts` (default 2,500) by that score, so
  every training transcript is well covered in every replicate.
- `all_covered`: every transcript with coverage ≥ `3_prime_bias_min_coverage` (default 1) in all
  IDs, with no cap.

One training set is shared by all IDs, so differences between the IDs' models come from the data,
not from different transcript sets. In genotype mode the set is chosen within each genotype.

## 4. Read depth

Bias is measured on **read depth** (every nucleotide each read covers), not on RT stops. Depth is
far less sparse than stops, which by design concentrate at reactive bases, so it estimates the
gradient with much less noise and doesn't depend on structure.

`samtools depth -a` (primary alignments, every position including zeros) over the training
transcripts is summed across:

- **+DMS:** the ID's own +DMS samples (all of them if `pool_replicates` is `plus` or `both`);
- **−DMS reference** (`3_prime_bias_reference`):
  - `pooled` (default): every −DMS sample after `exclude_samples_minus`. This gives the deepest,
    smoothest reference and the same reference for every ID.
  - `independent`: the ID's own −DMS samples. Use it if the −DMS libraries' gradients really
    differ between replicates.

## 5. The model

For 1-based position `p` on a transcript of length `L`:

```
d = L − p                    distance from the 3′ end (nt)
r = p / L                    relative position, 5′ → 3′
u = min(d, 2500) / 2500      capped, scaled distance from the 3′ end
v = min(p, 300) / 300        capped, scaled distance from the 5′ end

log(+DMS depth / −DMS depth) = b0  +  Σ_j β_j N_j(u)  +  b_r · r  +  Σ_m γ_m M_m(v)
```

- `N_j`: a natural cubic spline basis (restricted cubic spline, Harrell's parameterisation) on
  `u`, with knots at 0, 25, 54, 116, 250, 539, 1160 and 2500 nt from the 3′ end. Knots are
  spaced evenly on a log scale, so they are dense where the bias changes fastest (near the 3′ end)
  and sparse further away. It is flat beyond 2,500 nt, so very long transcripts don't extrapolate
  a trend.
- `M_m`: the same kind of basis on `v`, knots at 0, 10, 30, 100 and 300 nt from the 5′ end. It
  absorbs the drop in depth at the very 5′ end (reads can't extend past it, and the first bases
  are under-covered), so that drop isn't attributed to the 3′ term. It is flat beyond 300 nt.
- `b_r · r`: a linear term in relative position for any whole-transcript tilt.
- `b0`: the intercept, mostly the +DMS/−DMS library-size ratio. It doesn't matter for the
  correction because of the rescaling in §7.

The same coefficients apply to every transcript; only `p` and `L` differ. `exp(f)` is the
multiplicative bias at a position. Only ratios matter (§7 rescales the total): if `exp(f)` is
1.3 times higher at the 3′ end than at the 5′ end, +DMS has 30% more depth there relative to
−DMS, and the correction divides that difference out.

The pipeline always uses this spline form. (`bias_model.py` can still read the older linear form,
`b0 + b1·d + b2·r + b3·d·r`, from `model.tsv` files that don't have a `model_form` row.)

## 6. Fitting

Ordinary least squares over every training position where both +DMS and −DMS depth are above 0
(about 3 million positions for 2,500 transcripts). Positions with zero depth in either channel
have no log-ratio and are left out. The design columns are scaled to a maximum of 1 before
solving and the coefficients scaled back, because the spline columns span very different ranges.

`R²` in `model.tsv` is per nucleotide. Depth at single positions varies a lot (by sequence,
mappability and local structure), and the model only captures the smooth gradient, so **R² is
usually low (around 0.04 in the rice test data) even when the correction works.** Judge the fit by
the binned plots (§8), not by R².

## 7. Applying the correction

For every transcript in the transcriptome (mRNA, ncRNA and the positive control alike, not just
the training set), the +DMS stop count at `.rtsc` entry `i` (the stop at base `i`, see
[methods.md §6](methods.md#6-rt-stop-counting)) is divided by `exp(f(p = i, L))`.

All corrected counts are then multiplied by **one global constant** so that their total over
the whole transcriptome equals the raw total. This keeps the +DMS counts on their own library
scale. It matters because reactivity takes `ln(count + 1)` before normalizing, so shrinking every
count by the intercept would change how the log treats them. It doesn't change any position
relative to another. Corrected counts are no longer integers, and `rtsc_to_react` handles that.

Transcripts in the `.rtsc` but not in the FASTA are written unchanged. From here on, every
reactivity, p4p6 plot, ROC and depth-matched reactivity uses `combined_plus_corrected.rtsc`.

## 8. Outputs and how to read them

In `{output_dir}/{id}/3prime_bias_correction/` (also in the HTML report):

| File | What's in it |
|---|---|
| `model.tsv` | `model_form`, knots and caps, the coefficients, `r2`, `n_fit_positions`, `n_training_transcripts` |
| `coverage_residual_metagene.png` / `.tsv` | **Left:** mean log(+DMS / −DMS depth) per relative-position bin (orange), the model's fit (black dashed) and the residual (green). **Right:** the depth metagene (each transcript's depth / its own mean, averaged per bin over the training transcripts, ± SEM) of +DMS raw, +DMS corrected and −DMS. |
| `rtsc_metagene_before_after.png` / `.tsv` | The same metagene for RT stops: +DMS raw (left) and corrected (right) against −DMS, each with `sum((mean − ref)²)` over the bins |
| `{output_dir}/qc/3prime_bias_correction/training_transcripts.tsv` / `.bed` | The training set and each transcript's coverage in every ID |

A good correction looks like this:

- **Left panel:** the fit follows the binned log-ratio and the residual is flat near 0 along the
  whole transcript. A residual that still slopes means the model missed part of the gradient.
- **Right panel:** the corrected +DMS line lies on the −DMS line.
- **RT-stop plot:** `sum((mean − ref)²)` is lower after correction than before. The RT-stop
  metagene is noisier than depth and keeps real differences (+DMS stops sit at A/C), so it won't
  match −DMS exactly.

If the raw log-ratio is already flat, the libraries share their gradient and the correction does
little. That's fine, and you can leave it off.

To redraw an ID's figures from its saved TSVs without rerunning the correction:

```
python3 workflow/scripts/three_prime_end_bias/plot_bias.py {output_dir}/{id}/3prime_bias_correction \
    --title "rep1 (pooled -DMS reference): "
```

## 9. Assumptions and limits

- **The −DMS gradient is the right baseline.** The correction makes +DMS look like −DMS along
  the transcript. If the −DMS library itself is biased, that bias remains, but it cancels in the
  subtraction as long as the two libraries share it.
- **One shape for all transcripts.** The gradient is a function of position and length only. A
  transcript whose bias differs (unusual 3′ UTR, alternative polyadenylation, a different RNA
  class) is corrected with the average shape.
- **ncRNAs and the positive control are corrected with the mRNA model.** It's applied to them so
  every reactivity is built the same way, but their real bias may differ (they weren't
  poly(A)-selected in the same way). Check the p4p6 ROC with and without correction.
- **Depth bias is assumed to equal stop bias.** The model is fitted on depth and applied to stops.
- **Very short transcripts** have few training positions near both ends, and their correction
  rests on the spline's shape at small `d` and `p`.

## 10. Settings

| Key | Default | Meaning |
|---|---|---|
| `3_prime_bias_correction` | `false` | Turn the correction on |
| `3_prime_bias_reference` | `pooled` | −DMS reference depth: `pooled` (every −DMS sample) or `independent` (the ID's own) |
| `3_prime_bias_training` | `top_n` | Training set: `top_n` or `all_covered` |
| `3_prime_bias_n_transcripts` | 2500 | Training transcripts (`top_n`) |
| `3_prime_bias_min_coverage` | 1 | Minimum coverage in every ID (`all_covered`) |
| `3_prime_bias_annotation` | `annotation_gtf`, then `gff` | GTF/GFF3 used to find protein-coding transcripts |
| `3_prime_bias_bins` | 50 | Relative-position bins of the QC metagene plots (doesn't affect the fit) |

The knots and caps (`D_CAP`, `P_CAP`, `SPLINE_KNOTS_3P`, `SPLINE_KNOTS_5P`) are constants in
`workflow/scripts/three_prime_end_bias/bias_model.py`, and every `model.tsv` records the values
it was fitted with.
