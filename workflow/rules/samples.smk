############################################################################
# Samplesheet parsing and the sample/replicate-ID lists derived from it
# (exclusions, heat-correction temperature groups, pool_replicates IDs).
############################################################################

import pandas as pd


# Read in the samplesheet
samples_runs = pd.read_csv(config["samplesheet"], sep="\t")

# Make sure temperature is read as a string since temperature could be specified 
# as heat/ambient 
if "temperature" in samples_runs.columns:
    samples_runs["temperature"] = samples_runs["temperature"].astype(str)

samples = samples_runs.drop_duplicates(subset="sample", keep="first")

SAMPLES = samples["sample"].tolist()
MINUS_SAMPLES = samples.loc[samples["condition"] == "minus", "sample"].tolist()


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
# An optional 'temperature' samplesheet column marks each row's probing temperature. All rows sharing a replicate
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
    per-ID QC/summary targets in get_all_targets (targets.smk) -- coverage.csv,
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
