# Depth-matched reactivity

Optional (`depth_matching: true`; needs at least 2 reactivity IDs). For each transcript, every
library's +DMS RT stops are thinned at random to the same total, reactivity is recalculated, and
the result is averaged over many draws. This lets reactivities from libraries of different depth
be compared without depth itself showing up as a difference. Rules and files:
[pipeline.md §11](pipeline.md#11-optional-depth-matched-reactivity). The standard reactivity
calculation: [methods.md §8](methods.md#8-reactivity).

The standard `{output_dir}/{id}/reactivity.csv` is still produced. Depth matching adds
`{output_dir}/depth_matched/{id}/reactivity.csv` beside it.

## 1. Why depth changes reactivity

The reactivity of a transcript depends on how many +DMS stops it has, even with the same
structure:

- **The log.** `ln(count + 1)` compresses high counts. At low depth most positions are 0 or 1,
  while at high depth reactive positions separate from the rest, so the shape of the θ profile
  changes with depth.
- **Zeros.** At low depth, reactive bases that happened to get no reads are 0, and the 2–8%
  normalizer is taken from fewer, noisier non-zero values.
- **Noise.** Fewer stops means more sampling noise per base, which lowers correlations with other
  libraries and inflates apparent differences.

So when comparing two conditions, replicates or genotypes, a transcript that's twice as deep in
one library tends to differ in reactivity for that reason alone. Depth matching removes that by
giving every library the same number of +DMS stops on each transcript.

## 2. Method

Inputs, per reactivity ID: the +DMS `.rtsc` that the standard reactivity uses (3′ bias corrected
if that's on) and the ID's −DMS background (the pooled one with `pool_replicates: minus`).

For each transcript present in every ID's +DMS (with at least one stop in each), −DMS and the
FASTA:

1. **Target.** `target = min(+DMS total over the IDs)`, chosen per transcript, not per library.
2. **Round** each library's counts to integers. Bias-corrected counts aren't integers, so each
   position is rounded down and the remaining total is added back one count at a time, at
   positions drawn in proportion to their fractional parts. This is unbiased per position and
   keeps the total.
3. **Thin** each library to `target` with a multivariate hypergeometric draw: `target` stops
   drawn without replacement from that library's stops on that transcript. The library that
   already has `target` stays unchanged.
4. With `depth_matching_match_minus: true` (the default), each distinct −DMS background is thinned
   the same way to the smallest −DMS total. With one shared background (`pool_replicates: minus`)
   there's nothing to match.
5. **Recalculate reactivity** on the thinned counts exactly as in
   [methods.md §8](methods.md#8-reactivity): `ln(c + 1)`, scaling to the mean, −DMS subtraction,
   2–8% normalization over A/C, cap at 7. With `depth_matching_normalization: gu` (default: the
   `reactivity_normalization` value) each channel is scaled to its mean over the transcript's G/U
   positions instead of all positions
   ([methods.md §8, G/U normalization](methods.md#gu-normalization-reactivity_normalization-gu));
   a library gets no scale in a draw when the transcript has fewer than `depth_matching_min_gu`
   G/U positions or no stops at G/U in either channel.
6. **Repeat** `depth_matching_draws` times (default 100) and report each position's **mean**
   reactivity over the draws in which the transcript got a normalizer. A transcript is reported
   for a library only if at least half of the draws had one.

The random seed is fixed (20261014, with the transcript's index as a second seed), so reruns give
identical results.

In genotype mode, matched transcripts of the two genotypes share a key (from
`matched_transcripts.tsv`), so each pair is thinned to the same +DMS total in every library of
both genotypes. Transcripts without a match get no depth-matched reactivity.

## 3. Outputs

| File | What's in it |
|---|---|
| `{output_dir}/depth_matched/{id}/reactivity.csv` | `transcript, position, base, reactivity`, like the standard file: mean over draws, `NA` at G/U and the last base |
| `{output_dir}/depth_matched/{id}/matched_plus_draw1.rtsc` | The first draw's thinned +DMS stops, for RT-stop QC |
| `{output_dir}/depth_matched/summary.tsv` | Per transcript: `target_plus`, each ID's raw +DMS total and the number of draws that gave a reactivity |
| `{output_dir}/qc/depth_matched/rt_stop_replicate_correlation.png` | Replicate RT-stop correlation on the first draw (not in genotype mode) |
| `{output_dir}/qc/depth_matched/p4p6_roc.png`, `p4p6_auc.csv`, `p4p6_all_samples.png` | With the p4p6 control: the same ROC and structure plots on depth-matched reactivity |

With depth matching on, the HTML report's positive-control section shows the depth-matched
results, and genotype comparison uses depth-matched reactivity.

## 4. Things to know

- **Every library comes down to the shallowest.** Depth-matched reactivity of a deep library is
  noisier than its standard reactivity. That's the cost of a fair comparison; use the standard
  files when you look at one library alone.
- **Averaging smooths.** The mean over draws is smoother than any single draw, and its values sit
  below the cap more often. Compare depth-matched files only with other depth-matched files.
- **Per-transcript targets.** Two transcripts in one library are thinned to different totals, so
  depth matching makes libraries comparable, not transcripts.
- **Shallow transcripts drop out.** If one library has only a few stops on a transcript, every
  library is reduced to that, and the transcript may get no normalizer in most draws and go
  unreported. Look at `target_plus` in `summary.tsv` to see how deep the comparison really is.
- **Runtime.** It runs on every transcript shared by all IDs with 8 threads and 32 GB. Time grows
  with `depth_matching_draws`, the number of IDs and the number of transcripts (12 IDs, ~35k
  transcripts, 100 draws: ~40 min on 8 threads). On a cluster, `depth_matching_shards: N` runs it as
  N one-core jobs on contiguous blocks of the transcripts plus a merge; the merged output is
  identical to one job, since each transcript's draws are seeded by its index in the full list.

## 5. Settings

| Key | Default | Meaning |
|---|---|---|
| `depth_matching` | `false` | Turn depth matching on (needs ≥ 2 reactivity IDs) |
| `depth_matching_draws` | 100 | Thinning draws averaged per transcript |
| `depth_matching_match_minus` | `true` | Also thin the −DMS backgrounds to the smallest −DMS total |
| `depth_matching_roc_trim3` | 50 | 3′ nt left out of the depth-matched p4p6 ROC |
| `depth_matching_shards` | 1 | Split the work into this many one-core SLURM jobs (merged afterwards) |
| `depth_matching_normalization` | `reactivity_normalization` (`all`) | `all`: scale each channel to its mean over all positions (StructureFold2); `gu`: to its mean over the G/U positions |
| `depth_matching_min_gu` | `reactivity_min_gu` (20) | With `gu`: fewest G/U positions for a transcript to get a scale |
