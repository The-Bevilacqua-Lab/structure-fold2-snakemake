############################################################################
# Helper functions shared by the rule files (input functions, sample lookups).
# Kept separate from the rules -- see snakemake --lint's "mixed rules and
# functions" check.
############################################################################


# ---------------------------------------------------------------------------
# Sample sheet helpers: which runs / samples belong to what
# ---------------------------------------------------------------------------

def get_runs_for_sample(wildcards):
    """Return every sequencing run identifier for a sample, sorted."""
    rows = samples_runs.loc[samples_runs["sample"] == wildcards.sample]
    if rows.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in samplesheet")
    return sorted(rows["run"].unique().tolist())


def get_fastq_for_sample_run(wildcards):
    """Return the single raw fastq for one (sample, run) pair."""
    rows = samples_runs.loc[
        (samples_runs["sample"] == wildcards.sample) & (samples_runs["run"] == wildcards.run)
    ]
    if rows.empty:
        raise ValueError(f"No row for sample '{wildcards.sample}', run '{wildcards.run}' in samplesheet")
    return rows["r1"].iloc[0]


def _excluded_samples_for_condition(condition):
    """exclude_samples_plus/exclude_samples_minus (see the Snakefile) for one
    condition -- used to drop sample(s) from the +DMS/-DMS combining step
    while leaving them in every per-sample rule untouched."""
    return EXCLUDE_SAMPLES_PLUS if condition == "plus" else EXCLUDE_SAMPLES_MINUS


def get_samples_by_condition(condition):
    excluded = _excluded_samples_for_condition(condition)
    rows = samples.loc[(samples["condition"] == condition) & (~samples["sample"].isin(excluded))]
    if rows.empty:
        raise ValueError(
            f"No samples found for condition '{condition}' after applying "
            f"exclude_samples_{condition} (excluded: {sorted(excluded)})"
        )
    return rows["sample"].tolist()


def get_samples_by_condition_and_temperature(condition, temperature):
    excluded = _excluded_samples_for_condition(condition)
    rows = samples.loc[
        (samples["condition"] == condition)
        & (samples["temperature"] == temperature)
        & (~samples["sample"].isin(excluded))
    ]
    if rows.empty:
        raise ValueError(
            f"No samples found for condition '{condition}' and temperature '{temperature}' "
            f"after applying exclude_samples_{condition} (excluded: {sorted(excluded)})"
        )
    return rows["sample"].tolist()


# ---------------------------------------------------------------------------
# Transcriptome preparation helpers
# ---------------------------------------------------------------------------

def _source_transcriptome(wildcards):
    """Return the raw transcriptome FASTA before extra sequences are appended."""
    if TRANSCRIPT_SOURCE == "gffread":
        return f"{TMP}/resources/transcriptome_gffread_filtered.fa"
    return config["transcriptome"]


def _extra_transcriptome_fastas(wildcards):
    """
    Return extra FASTAs to append to the transcriptome before mapping:
    positive-control sequences (config key positive_control_fasta) and/or
    ncRNA sequences (config key ncrna_fasta). Either, both, or neither may
    be set.
    """
    extras = []
    if config.get("positive_control_fasta"):
        extras.append(config["positive_control_fasta"])
    if config.get("ncrna_fasta"):
        extras.append(config["ncrna_fasta"])
    return extras


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def get_run_sams_for_sample(wildcards):
    runs = get_runs_for_sample(wildcards)
    return [
        f"{TMP}/output/se_runs/{wildcards.sample}_{r}/{wildcards.sample}_{r}_trimmed_mapped.sam"
        for r in runs
    ]


def get_run_mapping_logs_for_sample(wildcards):
    runs = get_runs_for_sample(wildcards)
    return [
        f"{TMP}/output/se_runs/{wildcards.sample}_{r}/{wildcards.sample}_{r}_mapping_log.txt"
        for r in runs
    ]


# ---------------------------------------------------------------------------
# StructureFold2 (RT-stop counting / reactivity) helpers
# ---------------------------------------------------------------------------

def get_rtsc_counts_by_condition(condition):
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in get_samples_by_condition(condition)
    ]


def get_rtsc_counts_by_condition_and_id(condition, id):
    excluded = _excluded_samples_for_condition(condition)
    rows = samples.loc[
        (samples["condition"] == condition)
        & (samples["ID"] == id)
        & (~samples["sample"].isin(excluded))
    ]
    if rows.empty:
        raise ValueError(
            f"No samples found for condition '{condition}' and ID '{id}' after "
            f"applying exclude_samples_{condition} (excluded: {sorted(excluded)})"
        )
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in rows["sample"].tolist()
    ]


def get_rtsc_counts_by_condition_and_temperature(condition, temperature):
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in get_samples_by_condition_and_temperature(condition, temperature)
    ]


def get_minus_rtsc_for_reactivity(wildcards):
    if POOL_REPLICATES in ("minus", "both", "minus_all_plus_by_temperature"):
        return f"{TMP}/output/se/pooled/combined_minus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_minus.rtsc"


def get_raw_plus_rtsc(id):
    """The uncorrected +DMS .rtsc that feeds id's reactivity. A synthetic
    pooled/pooled_<temperature> ID always reads its own pooled file (this is
    also how the heat-correction lower-temperature scale asks for one)."""
    if not id.startswith("pooled") and POOL_REPLICATES in ("plus", "both"):
        return f"{TMP}/output/se/pooled/combined_plus.rtsc"
    return f"{TMP}/output/se/{id}/combined_plus.rtsc"


def get_plus_rtsc_for_id(id):
    """+DMS .rtsc used for id's reactivity: the 3' bias-corrected copy when
    3_prime_bias_correction is on (bias_correction.smk), else the raw one."""
    if BIAS_CORRECTION:
        return f"{BIAS_TMP}/{id}/combined_plus_corrected.rtsc"
    return get_raw_plus_rtsc(id)


def get_plus_rtsc_for_reactivity(wildcards):
    return get_plus_rtsc_for_id(wildcards.id)


# ---------------------------------------------------------------------------
# 3' bias correction helpers (workflow/rules/bias_correction.smk). The
# samples behind each channel mirror the combine_rtsc_* rules exactly
# (exclude_samples_plus/minus applied), so the depth the model is fit on
# comes from the same libraries as the RT-stops it corrects.
# ---------------------------------------------------------------------------

def _samples_for_condition_and_id(condition, id):
    if id == "pooled":
        return get_samples_by_condition(condition)
    if id.startswith("pooled_"):
        return get_samples_by_condition_and_temperature(condition, id[len("pooled_") :])
    excluded = _excluded_samples_for_condition(condition)
    rows = samples.loc[
        (samples["condition"] == condition)
        & (samples["ID"] == id)
        & (~samples["sample"].isin(excluded))
    ]
    if rows.empty:
        raise ValueError(
            f"No samples found for condition '{condition}' and ID '{id}' after "
            f"applying exclude_samples_{condition} (excluded: {sorted(excluded)})"
        )
    return rows["sample"].tolist()


def get_bias_samples(condition, id):
    """Samples whose BAMs give the +DMS (sample) or -DMS (reference) depth
    for id's 3' bias model."""
    if condition == "plus":
        if not id.startswith("pooled") and POOL_REPLICATES in ("plus", "both"):
            return get_samples_by_condition("plus")
        return _samples_for_condition_and_id("plus", id)
    if BIAS_REFERENCE == "pooled":
        return get_samples_by_condition("minus")
    return _samples_for_condition_and_id("minus", id)


def get_bias_bams(wildcards):
    return [
        f"{TMP}/output/se/{s}/{s}_trimmed_mapped_filtered_sorted.bam"
        for s in get_bias_samples(wildcards.condition, wildcards.id)
    ]


def get_bias_reference_minus_rtsc(wildcards):
    """-DMS RT-stops matching the reference depth, for the before/after plot."""
    if BIAS_REFERENCE == "pooled":
        return f"{TMP}/output/se/pooled/combined_minus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_minus.rtsc"


def get_sams_by_condition_and_id(condition, id):
    # "pooled" is the synthetic ID used only when pool_replicates == "both"
    # (see IDS in the main Snakefile) -- it's not a real ID in the
    # samplesheet, so fall back to every sample of that condition, mirroring
    # get_rtsc_counts_by_condition's pooling behavior. "pooled_<temperature>"
    # is the analogous synthetic ID for pool_replicates ==
    # "both_by_temperature", mirroring get_rtsc_counts_by_condition_and_temperature.
    if id == "pooled":
        rows = samples.loc[samples["condition"] == condition]
    elif id.startswith("pooled_"):
        temperature = id[len("pooled_") :]
        rows = samples.loc[
            (samples["condition"] == condition)
            & (samples["temperature"] == temperature)
        ]
    else:
        rows = samples.loc[(samples["condition"] == condition) & (samples["ID"] == id)]
    if rows.empty:
        raise ValueError(f"No samples found for condition '{condition}' and ID '{id}'")
    return [
        f"{TMP}/output/se/{s}/{s}_trimmed_mapped_filtered_sorted.bam"
        for s in rows["sample"].tolist()
    ]


def _specificity_comparison_inputs(wildcards):
    runs = config.get("specificity_comparison_runs", [])
    inputs = []
    for run in runs:
        sheet = pd.read_csv(run["samplesheet"], sep="\t")
        for sample in sheet["sample"].tolist():
            inputs.append(f"{run['output_dir']}/specificity_{sample}.csv")
    return inputs
