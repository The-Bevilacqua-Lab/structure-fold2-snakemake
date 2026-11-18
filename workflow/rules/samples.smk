############################################################################
# Samplesheet parsing and the sample / replicate-ID lists derived from it.
# Columns: docs/samplesheet.md.
############################################################################

import pandas as pd


samples_runs = pd.read_csv(config["samplesheet"], sep="\t", dtype=str)
samples_runs = samples_runs.apply(lambda col: col.str.strip())   # stray spaces after a path etc.

_missing = [c for c in ("sample", "condition", "r1") if c not in samples_runs.columns]
if _missing:
    raise ValueError(f"samplesheet {config['samplesheet']} is missing column(s) {_missing}")
_bad_condition = sorted(set(samples_runs["condition"]) - {"plus", "minus"})
if _bad_condition:
    raise ValueError(f"samplesheet 'condition' must be 'plus' or 'minus', got {_bad_condition}")

# Optional columns. run: one run per sample. ID: the n-th minus sample and the n-th plus sample
# (in sheet order) form replicate "rep<n>".
if "run" not in samples_runs.columns:
    samples_runs["run"] = "run1"
if "ID" not in samples_runs.columns:
    _order = samples_runs.drop_duplicates("sample")
    _rank = _order.groupby("condition").cumcount() + 1
    samples_runs["ID"] = samples_runs["sample"].map(dict(zip(_order["sample"], "rep" + _rank.astype(str))))

_dup = samples_runs[samples_runs.duplicated(["sample", "run"], keep=False)]
if not _dup.empty:
    raise ValueError(f"samplesheet has more than one row for (sample, run): "
                     f"{sorted(set(zip(_dup['sample'], _dup['run'])))}")
for _col in ("condition", "ID") + (("temperature",) if "temperature" in samples_runs.columns else ()):
    _mixed = samples_runs.groupby("sample")[_col].nunique()
    if (_mixed > 1).any():
        raise ValueError(f"sample(s) {_mixed[_mixed > 1].index.tolist()} have rows with different '{_col}'")

samples = samples_runs.drop_duplicates(subset="sample", keep="first")   # one row per sample
SAMPLES = samples["sample"].tolist()


def _validate_exclude_samples(key):
    excluded = config.get(key) or []
    if isinstance(excluded, str):
        excluded = [excluded]
    unknown = set(excluded) - set(SAMPLES)
    if unknown:
        raise ValueError(f"config[{key!r}] contains sample(s) not in the samplesheet: {sorted(unknown)}")
    return set(excluded)


# Left out of the +DMS / -DMS combining (and everything after it), but still aligned and QC'd.
EXCLUDE_SAMPLES_PLUS = _validate_exclude_samples("exclude_samples_plus")
EXCLUDE_SAMPLES_MINUS = _validate_exclude_samples("exclude_samples_minus")

# ---------------------------------------------------------------------------
# Heat correction: an optional 'temperature' column with exactly two values, one per replicate ID.
# ---------------------------------------------------------------------------
HEAT_CORRECTION = "temperature" in samples_runs.columns
if HEAT_CORRECTION:
    _id_temps = samples_runs.groupby("ID")["temperature"].nunique()
    if (_id_temps > 1).any():
        raise ValueError(
            f"replicate ID(s) {_id_temps[_id_temps > 1].index.tolist()} have more than one temperature; "
            "every row of an ID must share one")
    id_temperature = samples.drop_duplicates(subset="ID").set_index("ID")["temperature"]
    try:
        DISTINCT_TEMPERATURES = sorted(id_temperature.unique().tolist(), key=float)
    except ValueError:   # labels such as ambient/heat: sorted alphabetically
        DISTINCT_TEMPERATURES = sorted(id_temperature.unique().tolist())
    if len(DISTINCT_TEMPERATURES) != 2:
        raise ValueError(
            f"samplesheet 'temperature' column must have exactly 2 distinct values, found {DISTINCT_TEMPERATURES}")

# ---------------------------------------------------------------------------
# Reactivity IDs (IDS). pool_replicates (rtsc.smk) decides whether there is one per replicate or
# pooled ones: "both" -> "pooled"; "both_by_temperature" / "minus_all_plus_by_temperature" ->
# "pooled_<temperature>". "none", "minus" and "plus" keep one ID per replicate.
# ---------------------------------------------------------------------------
POOL_REPLICATES = config.get("pool_replicates", "none")
if POOL_REPLICATES == "both":
    if HEAT_CORRECTION:
        raise ValueError("pool_replicates: both pools across temperatures, leaving nothing to heat-correct; "
                         "use both_by_temperature")
    IDS = ["pooled"]
elif POOL_REPLICATES in ("both_by_temperature", "minus_all_plus_by_temperature"):
    if not HEAT_CORRECTION:
        raise ValueError(f"pool_replicates: {POOL_REPLICATES} needs a 'temperature' samplesheet column")
    IDS = [f"pooled_{t}" for t in DISTINCT_TEMPERATURES]
else:
    IDS = samples["ID"].unique().tolist()

# The real replicate IDs, whatever the pooling (replicate-vs-replicate QC uses these).
REPLICATE_IDS = samples["ID"].unique().tolist()
if set(REPLICATE_IDS) & {"pooled"} or any(i.startswith("pooled_") for i in REPLICATE_IDS):
    raise ValueError("replicate IDs 'pooled' and 'pooled_*' are reserved for pool_replicates; rename them")


def _replicate_ids_with_samples(condition):
    """REPLICATE_IDS that keep at least one sample of this condition after exclude_samples_*."""
    excluded = EXCLUDE_SAMPLES_PLUS if condition == "plus" else EXCLUDE_SAMPLES_MINUS
    keep = samples.loc[(samples["condition"] == condition) & ~samples["sample"].isin(excluded), "ID"]
    return [i for i in REPLICATE_IDS if i in set(keep)]


CORRELATION_PLUS_IDS = _replicate_ids_with_samples("plus")
CORRELATION_MINUS_IDS = _replicate_ids_with_samples("minus")


def _ids_with_samples(condition):
    """IDS whose own combined_<condition>.rtsc exists (per-ID QC reads that file, never a pooled
    stand-in). Pooled IDs always do; a real ID drops out if all its samples were excluded."""
    usable = set(_replicate_ids_with_samples(condition))
    return [i for i in IDS if i not in REPLICATE_IDS or i in usable]


IDS_WITH_PLUS = _ids_with_samples("plus")
IDS_WITH_MINUS = _ids_with_samples("minus")

# True when pooling replaced the replicate IDs (both / *_by_temperature)
MERGED_IDS_DIFFER = IDS != REPLICATE_IDS

if HEAT_CORRECTION:
    if POOL_REPLICATES in ("both_by_temperature", "minus_all_plus_by_temperature"):
        HEAT_CORRECTION_LOWER_IDS = [f"pooled_{DISTINCT_TEMPERATURES[0]}"]
        HEAT_CORRECTION_HIGHER_IDS = [f"pooled_{DISTINCT_TEMPERATURES[1]}"]
    else:
        HEAT_CORRECTION_LOWER_IDS = id_temperature[id_temperature == DISTINCT_TEMPERATURES[0]].index.tolist()
        HEAT_CORRECTION_HIGHER_IDS = id_temperature[id_temperature == DISTINCT_TEMPERATURES[1]].index.tolist()
    HEAT_CORRECTION_SUFFIX = config.get("heat_correction_suffix", "corrected")
    # Also heat-correct every transcript in both temperatures' reactivity.react, without the coverage filter
    HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS = _flag(config.get("heat_correction_compare_all_transcripts"), False)

# ---------------------------------------------------------------------------
# Genotype mode (config genotypes, settings.smk): a 'genotype' column assigns every row to one of
# the two genotypes; all rows of a replicate ID share it.
# ---------------------------------------------------------------------------
if GENOTYPE_MODE:
    if "genotype" not in samples_runs.columns:
        raise ValueError("config['genotypes'] is set, so the samplesheet needs a 'genotype' column")
    _unknown = sorted(set(samples_runs["genotype"]) - set(GENOTYPE_NAMES))
    if _unknown:
        raise ValueError(f"samplesheet genotype(s) {_unknown} are not in config['genotypes'] ({GENOTYPE_NAMES})")
    for _col in ("sample", "ID"):
        _mixed = samples_runs.groupby(_col)["genotype"].nunique()
        if (_mixed > 1).any():
            raise ValueError(f"{_col}(s) {_mixed[_mixed > 1].index.tolist()} have rows from more than one genotype")
    if HEAT_CORRECTION:
        raise ValueError("heat correction (a 'temperature' column) isn't supported together with config['genotypes']")
    if POOL_REPLICATES not in ("none", "minus"):
        raise ValueError(
            f"pool_replicates {POOL_REPLICATES!r} would pool across genotypes; with config['genotypes'] use "
            "'none' or 'minus' (-DMS is pooled within each genotype; add pool_minus_across_genotypes: true "
            "to also pool it across genotypes)")
    samples = samples_runs.drop_duplicates(subset="sample", keep="first")
    SAMPLE_GENOTYPE = samples.set_index("sample")["genotype"].to_dict()
    ID_GENOTYPE = samples.drop_duplicates(subset="ID").set_index("ID")["genotype"].to_dict()
    for _g in GENOTYPE_NAMES:
        for _cond in ("plus", "minus"):
            if not ((samples["genotype"] == _g) & (samples["condition"] == _cond)).any():
                raise ValueError(f"genotype {_g!r} has no {_cond} samples in the samplesheet")
    # replicate QC (RT-stop correlation, UpSet) is per genotype: the transcripts have different names
    GENOTYPE_PLUS_IDS = {g: [i for i in CORRELATION_PLUS_IDS if ID_GENOTYPE[i] == g] for g in GENOTYPE_NAMES}
    GENOTYPE_MINUS_IDS = {g: [i for i in CORRELATION_MINUS_IDS if ID_GENOTYPE[i] == g] for g in GENOTYPE_NAMES}
else:
    SAMPLE_GENOTYPE, ID_GENOTYPE = {}, {}
