###############################################################################
# structure-fold2-snakemake
#
# A generalized Snakemake wrapper around StructureFold2 for DMS-MaPseq /
# Structure-seq-style +DMS/-DMS reactivity calculation on any organism.
#
# Minimal usage:
#   snakemake --sdm conda -j 8 --configfile config/config.yaml
#
# On a SLURM cluster:
#   snakemake --sdm conda -j 100 --executor cluster-generic \
#     --cluster-generic-submit-cmd 'sbatch --time=4:00:00 --ntasks=10 --mem=40gb --partition=standard --account=pcb5_cr_default \
#       --output=slurm_logs/%j.out' \
#     --configfile config/config_in_vitro_nipponbare_23C.yaml --latency-wait 60 --conda-prefix ~/group/kjk6173/local
#
# See config/config.yaml for every available config key, and README.md for
# an overview of what this pipeline does and does not generalize.
###############################################################################

from snakemake.utils import min_version
import pandas as pd
from snakemake.shell import shell

# On HPC systems that use Environment Modules, a loaded module (e.g. an R or
# Python module) can prepend its own lib dir onto LD_LIBRARY_PATH, shadowing
# the newer libstdc++ shipped inside conda envs and breaking compiled
# extensions. `module load gcc` first (a no-op, and harmless, on any system
# without an Environment Modules setup -- the `|| true` swallows the "module:
# command not found" case) gets a modern libstdc++ onto LD_LIBRARY_PATH, then
# the conda env's own lib dir (once activated) takes final precedence.
# Re-declaring "set -euo pipefail;" here is required because shell.prefix()
# overrides Snakemake's default rather than extending it.
shell.prefix(
    "set -euo pipefail; "
    "module load gcc >/dev/null 2>&1 || true; "
    'if [ -n "${{CONDA_PREFIX:-}}" ]; then '
    'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"; '
    "fi; "
)


configfile: "config/config.yaml"


container: "docker://continuumio/miniconda3"


# ---------------------------------------------------------------------------
# Samplesheet
#
# One row per (sample, sequencing run) pair, tab-separated, with columns:
#   sample     – sample identifier (matches fastq naming)
#   condition  – exactly "plus" (+DMS/treated) or "minus" (-DMS/untreated)
#   r1         – path to that (sample, run)'s raw FASTQ
#   ID         – optional; groups samples into biological replicates (e.g.
#                "rep1") -- the reactivity calculation pairs one condition's
#                ID against the other's, so a given ID needs both a "plus"
#                and a "minus" row somewhere in the sheet. If the column is
#                absent, ID defaults to each sample's row position *within
#                its condition* ("rep1", "rep2", ... independently for plus
#                and minus), i.e. the 1st minus row is paired with the 1st
#                plus row as one replicate, the 2nd with the 2nd, and so on
#                -- this assumes equal plus/minus sample counts in that
#                order; give an explicit ID column if that's not your layout.
#   run        – optional; distinguishes multiple sequencing runs of the same
#                sample (resequenced on a second flow cell, etc.), which get
#                concatenated before trimming. Defaults to a single run if
#                the column is absent.
#   temperature – optional; the probing temperature for that row (e.g. "30",
#                "55", or labels like "low"/"high" -- any two distinct
#                values). Enables heat correction (react_heat_correct.py,
#                see workflow/rules/heat_correction.smk) between the two
#                temperature groups, and unlocks pool_replicates:
#                both_by_temperature (see below). Every row of a given ID
#                (its plus row, minus row, and any extra runs) must share
#                the same temperature, since they're the same replicate
#                just probed differently. Exactly two distinct temperature
#                values must appear across the whole samplesheet. Omit the
#                column entirely to skip heat correction.
# ---------------------------------------------------------------------------
samples_runs = pd.read_csv(config["samplesheet"], sep="\t")
if "run" not in samples_runs.columns:
    samples_runs["run"] = "run1"
if "ID" not in samples_runs.columns:
    # Compute once per unique sample (not per run-row, which would
    # over-count multi-run samples), ranked by row order within each
    # condition, then broadcast back across that sample's run-rows.
    sample_condition = samples_runs.drop_duplicates(subset="sample")[
        ["sample", "condition"]
    ].copy()
    sample_condition["ID"] = "rep" + (
        sample_condition.groupby("condition").cumcount() + 1
    ).astype(str)
    samples_runs = samples_runs.merge(
        sample_condition[["sample", "ID"]], on="sample", how="left"
    )
if "temperature" in samples_runs.columns:
    # Normalize to str so comparisons/dict keys/filenames built from this
    # column are consistent regardless of whether pandas inferred it as
    # numeric (e.g. a sheet with only "30"/"55") or left it as text labels
    # (e.g. "ambient"/"heat").
    samples_runs["temperature"] = samples_runs["temperature"].astype(str)

samples = samples_runs.drop_duplicates(subset="sample", keep="first")

SAMPLES = samples["sample"].tolist()
MINUS_SAMPLES = samples.loc[samples["condition"] == "minus", "sample"].tolist()

# ---------------------------------------------------------------------------
# Optional: exclude specific sample(s) from the +DMS or -DMS combining step
# (combine_rtsc_plus/minus and their _pooled/_pooled_by_temperature variants
# in workflow/rules/rtsc.smk) without dropping them from the
# samplesheet entirely. An excluded sample is still trimmed, mapped, and
# rtsc-counted individually, and still shows up in per-sample QC
# (specificity_<sample>.csv, alignment_stats_summary.tsv) -- it's just left
# out of whichever channel's combined_plus.rtsc/combined_minus.rtsc feeds
# the reactivity calculation. Useful for dropping a technically bad
# replicate (e.g. one flagged by
# notebooks/dms_rtsc_normalization_diagnostics.ipynb) from the -DMS
# background, or a bad +DMS sample from the treated signal, without losing
# its QC trail. Each key is independent -- excluding a sample from
# exclude_samples_minus does not affect its use in +DMS, and vice versa.
# Omit either/both keys (or leave as an empty list) to exclude nothing.
#
# If excluding sample(s) leaves a replicate ID with ZERO remaining samples
# for a channel (e.g. that ID only had one -DMS sample and it got
# excluded), that ID is dropped from:
#   - CORRELATION_PLUS_IDS/CORRELATION_MINUS_IDS below, and so skipped by
#     the always-per-replicate correlation QC
#     (rt_stop_replicate_correlation[_minus].png) and the +DMS coverage
#     UpSet plot;
#   - IDS_WITH_PLUS/IDS_WITH_MINUS below, and so skipped by the other
#     per-ID targets that always read the true per-ID channel regardless of
#     pool_replicates -- coverage.csv, counts_minus.csv,
#     specificity_plus/minus.csv, abundance_<mode>.csv;
# since there's nothing left for either to read on that channel for that
# ID. Its reactivity output still gets computed normally if pool_replicates
# provides a pooled background for that channel (e.g.
# "minus"/"both"/"minus_all_plus_by_temperature" for a fully-excluded -DMS
# replicate) -- reactivity.react/.csv is the one per-ID target that CAN
# fall back this way, via get_minus_rtsc_for_reactivity/
# get_plus_rtsc_for_reactivity; otherwise combine_rtsc_plus/minus itself
# raises a clear ValueError naming the empty ID, since there's no signal
# left to combine at all.
# ---------------------------------------------------------------------------
def _validate_exclude_samples(key):
    excluded = config.get(key, [])
    if isinstance(excluded, str):
        excluded = [excluded]
    unknown = set(excluded) - set(SAMPLES)
    if unknown:
        raise ValueError(
            f"config[{key!r}] contains sample(s) not in the samplesheet: {sorted(unknown)}"
        )
    return set(excluded)


EXCLUDE_SAMPLES_PLUS = _validate_exclude_samples("exclude_samples_plus")
EXCLUDE_SAMPLES_MINUS = _validate_exclude_samples("exclude_samples_minus")

# ---------------------------------------------------------------------------
# Optional heat correction / temperature-aware pooling
#
# An optional 'temperature' samplesheet column (see the samplesheet comment
# above) marks each row's probing temperature. All rows sharing a replicate
# ID must agree on temperature, and the sheet must have exactly 2 distinct
# temperature values -- matching react_heat_correct.py's lower/higher
# two-condition design (workflow/rules/heat_correction.smk). The same
# grouping also backs pool_replicates: both_by_temperature below.
# ---------------------------------------------------------------------------
HEAT_CORRECTION = "temperature" in samples_runs.columns
if HEAT_CORRECTION:
    id_temp_counts = samples_runs.groupby("ID")["temperature"].nunique()
    inconsistent_ids = id_temp_counts[id_temp_counts > 1].index.tolist()
    if inconsistent_ids:
        raise ValueError(
            f"Replicate ID(s) {inconsistent_ids} have more than one distinct "
            "'temperature' value across their rows -- every row (plus/minus, "
            "all runs) of a given ID must share the same temperature."
        )

    id_temperature = samples.drop_duplicates(subset="ID").set_index("ID")["temperature"]
    try:
        DISTINCT_TEMPERATURES = sorted(id_temperature.unique().tolist(), key=float)
    except (TypeError, ValueError):
        DISTINCT_TEMPERATURES = sorted(id_temperature.unique().tolist(), key=str)
    if len(DISTINCT_TEMPERATURES) != 2:
        raise ValueError(
            "samplesheet 'temperature' column must have exactly 2 distinct "
            f"values, found {DISTINCT_TEMPERATURES}"
        )

# See workflow/rules/rtsc.smk for the full explanation of what
# each pool_replicates value does to the reactivity calculation.
POOL_REPLICATES = config.get("pool_replicates", "none")

# "both" collapses every replicate's -DMS and +DMS into one pooled pair
# (across the WHOLE samplesheet, temperature included), so there is exactly
# one reactivity output instead of one per replicate ID -- every per-ID rule
# in the pipeline (coverage, specificity, QC plots...) then runs once, on
# that single pooled dataset, under the synthetic ID "pooled". That makes it
# incompatible with heat correction, since pooling across temperature
# destroys the very distinction being corrected for -- use
# "both_by_temperature" instead, which pools +DMS/-DMS within each
# temperature separately, producing one pooled ID per temperature
# ("pooled_<temperature>"), still leaving the two temperatures to
# heat-correct against each other. "minus_all_plus_by_temperature" pools
# -DMS across the whole sheet (temperature included) while still pooling
# +DMS only within each temperature -- same one-ID-per-temperature shape as
# "both_by_temperature", just with a broader -DMS background. Any other
# pool_replicates value leaves IDS as the real replicate-group list.
if POOL_REPLICATES == "both":
    if HEAT_CORRECTION:
        raise ValueError(
            "samplesheet has a 'temperature' column but pool_replicates is "
            "'both' -- that mode pools across every replicate AND every "
            "temperature into one output, leaving nothing to heat-correct. "
            "Use 'both_by_temperature' to pool within each temperature "
            "separately instead."
        )
    IDS = ["pooled"]
elif POOL_REPLICATES in ("both_by_temperature", "minus_all_plus_by_temperature"):
    if not HEAT_CORRECTION:
        raise ValueError(
            f"config['pool_replicates'] is {POOL_REPLICATES!r} but the "
            "samplesheet has no 'temperature' column to pool within."
        )
    IDS = [f"pooled_{t}" for t in DISTINCT_TEMPERATURES]
else:
    IDS = samples["ID"].unique().tolist()

# The true biological-replicate ID list, regardless of pool_replicates --
# unlike IDS, this never collapses to a synthetic "pooled"/"pooled_<temperature>"
# entry. Used by QC that specifically needs to compare replicates against
# each other (e.g. rtsc_stop_correlation/_minus in replicate_correlation.smk),
# since that comparison is meaningless once replicates have already been
# merged together.
REPLICATE_IDS = samples["ID"].unique().tolist()


def _replicate_ids_with_samples(condition):
    """REPLICATE_IDS subset that still has >=1 non-excluded sample for the
    given condition, after exclude_samples_plus/exclude_samples_minus.
    Backs the "always per-replicate, regardless of pool_replicates" QC --
    rtsc_stop_correlation/_minus and upset_covered_transcripts_replicates
    (workflow/rules/replicate_correlation.smk) -- which needs a real, non-empty
    combine_rtsc_plus/minus output for every ID it includes. A replicate
    whose entire +DMS or -DMS channel was excluded has nothing left to
    correlate/cover on that channel, even though it may still get a
    reactivity output via a pooled background (see
    get_minus_rtsc_for_reactivity/get_plus_rtsc_for_reactivity)."""
    excluded = EXCLUDE_SAMPLES_PLUS if condition == "plus" else EXCLUDE_SAMPLES_MINUS
    return [
        id_
        for id_ in REPLICATE_IDS
        if not samples.loc[
            (samples["condition"] == condition) & (samples["ID"] == id_), "sample"
        ]
        .isin(excluded)
        .all()
    ]


CORRELATION_PLUS_IDS = _replicate_ids_with_samples("plus")
CORRELATION_MINUS_IDS = _replicate_ids_with_samples("minus")


def _ids_with_samples(condition):
    """IDS subset that still has a real, non-empty combine_rtsc_plus/minus
    output for the given channel. Unlike reactivity.react/.csv (which can
    fall back to a pooled background via get_minus_rtsc_for_reactivity/
    get_plus_rtsc_for_reactivity when pool_replicates provides one), several
    per-ID QC/summary targets in get_all_targets below -- coverage.csv,
    counts_minus.csv, specificity_plus/minus.csv, abundance_<mode>.csv --
    always read the TRUE per-ID channel regardless of pool_replicates (see
    calculate_stop_coverage/rtsc_total_stops/calculate_specificity/
    calculate_transcript_abundance in workflow/rules/), so
    they can't fall back the same way. A synthetic pooled/pooled_<temperature>
    ID always passes through here unfiltered -- its own
    combine_rtsc_*_pooled[_by_temperature] rule already applies
    exclude_samples_plus/minus and raises its own clear error if that whole
    pooled group ends up empty; only a REAL replicate ID (in REPLICATE_IDS)
    whose own channel was fully excluded needs to be dropped here."""
    usable_real = set(_replicate_ids_with_samples(condition))
    return [id_ for id_ in IDS if id_ not in REPLICATE_IDS or id_ in usable_real]


IDS_WITH_PLUS = _ids_with_samples("plus")
IDS_WITH_MINUS = _ids_with_samples("minus")

# True when pool_replicates has actually collapsed IDS to something other
# than the true replicate list (i.e. "both", "both_by_temperature", or
# "minus_all_plus_by_temperature" -- see IDS above) -- "minus"/"plus"
# channel-only pooling leaves IDS == REPLICATE_IDS, since it only changes
# which .rtsc feeds the reactivity calculation, not the per-ID coverage
# used by upset_covered_transcripts_merged (workflow/rules/coverage.smk).
MERGED_IDS_DIFFER = IDS != REPLICATE_IDS

if HEAT_CORRECTION:
    if POOL_REPLICATES in ("both_by_temperature", "minus_all_plus_by_temperature"):
        # Pooling already collapsed each temperature to a single ID, so
        # there's exactly one "replicate" per temperature left to correct.
        HEAT_CORRECTION_LOWER_IDS = [f"pooled_{DISTINCT_TEMPERATURES[0]}"]
        HEAT_CORRECTION_HIGHER_IDS = [f"pooled_{DISTINCT_TEMPERATURES[1]}"]
    else:
        HEAT_CORRECTION_LOWER_IDS = id_temperature[
            id_temperature == DISTINCT_TEMPERATURES[0]
        ].index.tolist()
        HEAT_CORRECTION_HIGHER_IDS = id_temperature[
            id_temperature == DISTINCT_TEMPERATURES[1]
        ].index.tolist()
    HEAT_CORRECTION_SUFFIX = config.get("heat_correction_suffix", "corrected")
    # Optional comparison run: also heat-correct every transcript present in
    # both temperature conditions' plain reactivity.react, with no RT-stop
    # coverage threshold at all (see heat_correction_all_transcripts_shared_react
    # in workflow/rules/heat_correction.smk), alongside the normal
    # coverage-qualified correction -- to see how much the coverage
    # threshold itself changes the correction.
    HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS = config.get(
        "heat_correction_compare_all_transcripts", False
    )


# ---------------------------------------------------------------------------
# Optional: estimate per-transcript relative abundance (RPKM and/or TPM)
# from the untreated (-DMS) RT-stop counts, via StructureFold2's own
# rtsc_abundances.py (see calculate_transcript_abundance in
# workflow/rules/abundance.smk) -- same -DMS-channel rationale as
# counts_minus.csv (rtsc_total_stops.py's docstring), just converted to a
# length/library-size-normalized abundance metric instead of DESeq2-style
# raw counts. One {id}/abundance_<mode>.csv per requested mode. Omit
# entirely (or leave as an empty list) to skip this output.
# ---------------------------------------------------------------------------
_ALLOWED_ABUNDANCE_MODES = {"RPKM", "TPM"}
ABUNDANCE_MODES = config.get("transcript_abundance", [])
if isinstance(ABUNDANCE_MODES, str):
    ABUNDANCE_MODES = [ABUNDANCE_MODES]
ABUNDANCE_MODES = sorted({m.upper() for m in ABUNDANCE_MODES})
_invalid_abundance_modes = set(ABUNDANCE_MODES) - _ALLOWED_ABUNDANCE_MODES
if _invalid_abundance_modes:
    raise ValueError(
        f"config['transcript_abundance'] must only contain {sorted(_ALLOWED_ABUNDANCE_MODES)}, "
        f"got {sorted(_invalid_abundance_modes)}"
    )


include: "workflow/rules/settings.smk"
include: "workflow/rules/helpers.smk"
include: "workflow/rules/transcriptome.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/read_prep.smk"
include: "workflow/rules/sam_processing.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/rtsc.smk"
include: "workflow/rules/reactivity.smk"
include: "workflow/rules/heat_correction.smk"
include: "workflow/rules/coverage.smk"
include: "workflow/rules/abundance.smk"
include: "workflow/rules/specificity.smk"
include: "workflow/rules/replicate_correlation.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/metagene.smk"


# The p4p6 structure plot rules depend on a hand-annotated coordinate map
# specific to that one construct (see workflow/rules/positive_control.smk) --
# only wired in when the positive control actually is p4p6. Any other
# positive control still gets the generic alignment-% QC from qc.smk.
POSITIVE_CONTROL_IS_P4P6 = config.get("positive_control_name") == "p4p6"
if POSITIVE_CONTROL_IS_P4P6:

    include: "workflow/rules/positive_control.smk"


def get_all_targets(wildcards):
    out = config["output_dir"]
    targets = (
        expand("{out}/{id}/reactivity.react", out=out, id=IDS)
        + expand("{out}/{id}/reactivity.csv", out=out, id=IDS)
        # coverage.csv/specificity_plus.csv always read the true per-ID +DMS
        # channel (no pooled fallback -- see IDS_WITH_PLUS above), so an ID
        # whose entire +DMS was excluded via exclude_samples_plus is skipped
        # here rather than requesting an impossible target.
        + expand("{out}/{id}/coverage.csv", out=out, id=IDS_WITH_PLUS)
        + expand("{out}/{id}/specificity_plus.csv", out=out, id=IDS_WITH_PLUS)
        # counts_minus.csv/specificity_minus.csv: same, for the -DMS channel
        # (IDS_WITH_MINUS) -- this is what exclude_samples_minus actually
        # needed here (see IDS_WITH_MINUS above).
        + expand("{out}/{id}/counts_minus.csv", out=out, id=IDS_WITH_MINUS)
        + expand("{out}/{id}/specificity_minus.csv", out=out, id=IDS_WITH_MINUS)
        + expand("{out}/specificity_{sample}.csv", out=out, sample=SAMPLES)
        + [
            f"{out}/qc/alignment_stats_summary.tsv",
            f"{out}/qc/specificity_plot.png",
            f"{out}/qc/metagene/metagene_depth.png",
            f"{out}/qc/metagene/metagene_rtsc.png",
        ]
        + (
            # Correlation is computed on CORRELATION_PLUS_IDS (true
            # biological replicates that still have a non-excluded +DMS
            # sample), not IDS, so gate on that count too -- otherwise a
            # pooled run (IDS == ["pooled"] or ["pooled_<t>", ...]) would
            # skip this QC even when there are multiple real replicates to
            # compare.
            [f"{out}/qc/rt_stop_replicate_correlation.png"]
            if len(CORRELATION_PLUS_IDS) > 1
            else []
        )
        + (
            # Same, for the -DMS channel -- gated separately since
            # exclude_samples_plus/exclude_samples_minus can leave the two
            # channels with different replicate counts.
            [f"{out}/qc/rt_stop_replicate_correlation_minus.png"]
            if len(CORRELATION_MINUS_IDS) > 1
            else []
        )
        + (
            # Always produced (per-replicate, unpooled -- see
            # upset_covered_transcripts_replicates in coverage.smk,
            # which depends only on each replicate's +DMS channel), same
            # gate as the +DMS correlation plot above: an UpSet plot of a
            # single set has nothing to overlap.
            [f"{out}/qc/covered_transcripts_upset_replicates.png"]
            if len(CORRELATION_PLUS_IDS) > 1
            else []
        )
        + (
            # Only when pool_replicates has actually merged replicates into
            # a different (and still >1-set) grouping -- see
            # MERGED_IDS_DIFFER above and upset_covered_transcripts_merged
            # in coverage.smk. "both" pooling collapses to the single
            # ID "pooled", which can't be UpSet-plotted against anything.
            [f"{out}/qc/covered_transcripts_upset_merged.png"]
            if MERGED_IDS_DIFFER and len(IDS) > 1
            else []
        )
        + (
            [f"{out}/transcript_position_annotations.csv"]
            if config.get("annotation_gtf") and config.get("genome")
            else []
        )
    )

    if ABUNDANCE_MODES:
        # abundance_<mode>.csv is computed from the true per-ID -DMS channel
        # (see calculate_transcript_abundance in abundance.smk), same
        # IDS_WITH_MINUS gating as counts_minus.csv above.
        targets += expand(
            "{out}/{id}/abundance_{mode}.csv", out=out, id=IDS_WITH_MINUS, mode=ABUNDANCE_MODES
        )

    if config.get("positive_control_fasta"):
        targets += [f"{out}/qc/positive_control_alignment_summary.tsv"]

    if POSITIVE_CONTROL_IS_P4P6:
        targets += (
            expand("{out}/{id}/p4p6_react.png", out=out, id=IDS)
            + [
                f"{out}/qc/p4p6_all_samples.png",
            ]
        )

    if config.get("specificity_comparison_runs"):
        targets += [f"{out}/qc/specificity_comparison.png"]

    if TRIMMER == "fastp":
        targets += [f"{out}/qc/multiqc_fastp/multiqc_report.html"]

    if HEAT_CORRECTION:
        targets += expand(
            "{out}/{id}/heat_correction/reactivity_{suffix}.{ext}",
            out=out,
            id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            suffix=HEAT_CORRECTION_SUFFIX,
            ext=["react", "csv"],
        ) + expand(
            "{out}/{id}/heat_correction/reactivity.csv",
            out=out,
            id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
        )

        if HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS:
            targets += expand(
                "{out}/{id}/heat_correction_all_transcripts/reactivity_{suffix}.{ext}",
                out=out,
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                suffix=HEAT_CORRECTION_SUFFIX,
                ext=["react", "csv"],
            )

    return targets


rule all:
    input:
        get_all_targets,
