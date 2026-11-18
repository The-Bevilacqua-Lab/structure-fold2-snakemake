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
    """exclude_samples_plus / exclude_samples_minus (samples.smk) for one condition."""
    return EXCLUDE_SAMPLES_PLUS if condition == "plus" else EXCLUDE_SAMPLES_MINUS


def samples_to_combine(condition, id):
    """
    The samples whose RT-stops make up {id}/combined_<condition>.rtsc, leaving out exclude_samples_*:
    one replicate ID's samples, every sample (id "pooled"), or every sample of one temperature
    (id "pooled_<temperature>").
    """
    excluded = _excluded_samples_for_condition(condition)
    rows = samples.loc[(samples["condition"] == condition) & ~samples["sample"].isin(excluded)]
    if id == "pooled":
        pass
    elif id.startswith("pooled_"):
        rows = rows.loc[rows["temperature"] == id[len("pooled_"):]]
    else:
        rows = rows.loc[rows["ID"] == id]
    if rows.empty:
        raise ValueError(
            f"No {condition} samples for ID '{id}' after exclude_samples_{condition} (excluded: {sorted(excluded)})")
    return rows["sample"].tolist()


def get_samples_by_condition(condition):
    """Every sample of one condition, leaving out exclude_samples_*."""
    return samples_to_combine(condition, "pooled")


def sample_rtsc(sample):
    return f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"


def sample_bam(sample):
    return f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered_sorted.bam"


# ---------------------------------------------------------------------------
# Transcriptome preparation helpers
# ---------------------------------------------------------------------------

def _source_transcriptome(wildcards):
    """Return the raw transcriptome FASTA before extra sequences are appended."""
    if GENOTYPE_MODE:
        # both genotypes' transcripts, renamed <genotype>__<id> (genotypes.smk)
        return f"{GENOTYPE_DIR}/all_genotypes_transcripts.fa"
    if TRANSCRIPT_SOURCE == "gffread":
        if CDS_ONLY:
            return f"{TMP}/resources/transcriptome_gffread_filtered.fa"
        return f"{TMP}/resources/transcriptome_gffread.fa"
    return config["transcriptome"]


def _extra_transcriptome_fastas(wildcards):
    """
    Return extra FASTAs to append to the transcriptome before mapping:
    positive-control sequences (config key positive_control_fasta), rRNA, ncRNA
    and lncRNA sequences (EXTRA_RNA, settings.smk: rrna_/ncrna_/lncrna_ fasta,
    or gff/gtf extracted with gffread). Any combination may be set.
    """
    extras = []
    if config.get("positive_control_fasta"):
        extras.append(config["positive_control_fasta"])
    for kind in EXTRA_RNA_CATEGORY:
        if kind in EXTRA_RNA:
            extras.append(EXTRA_RNA[kind]["fasta"])
    return extras


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def get_bowtie2_index_for_sample(wildcards):
    """Bowtie2 index (.1.bt2) a sample's reads map to: its genotype's own transcripts (+ shared
    positive control / ncRNAs) in genotype mode, else the one transcriptome."""
    if GENOTYPE_MODE:
        return f"{GENOTYPE_DIR}/{SAMPLE_GENOTYPE[wildcards.sample]}/mapping.fa.1.bt2"
    return f"{TRANSCRIPTOME}.1.bt2"


def get_rtsc_fasta_for_sample(wildcards):
    """Uppercase FASTA sam_to_rtsc writes a sample's RT-stops over (every transcript it lists):
    only the sample's own genotype in genotype mode, so its .rtsc never lists the other genotype."""
    if GENOTYPE_MODE:
        return f"{GENOTYPE_DIR}/{SAMPLE_GENOTYPE[wildcards.sample]}/mapping_upper.fa"
    return TRANSCRIPTOME_UPPER


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

def genotype_pooled_minus_rtsc(genotype):
    """Genotype mode: every -DMS sample of one genotype pooled (combine_rtsc_minus_pooled_genotype)."""
    return f"{TMP}/output/se/pooled/{genotype}/combined_minus.rtsc"


def get_minus_rtsc_for_reactivity(wildcards):
    if POOL_MINUS_ACROSS_GENOTYPES:
        # -DMS pooled within each genotype + the other genotype's at aligned positions (genotypes.smk)
        return CROSS_GENOTYPE_MINUS
    if GENOTYPE_MODE and POOL_REPLICATES == "minus":
        # the ID's own genotype's -DMS only, also on the shared positive control / ncRNAs / lncRNAs
        return genotype_pooled_minus_rtsc(ID_GENOTYPE[wildcards.id])
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

def get_bias_samples(condition, id):
    """Samples whose BAMs give the +DMS (sample) or -DMS (reference) depth
    for id's 3' bias model."""
    if condition == "plus":
        if not id.startswith("pooled") and POOL_REPLICATES in ("plus", "both"):
            return get_samples_by_condition("plus")
        return samples_to_combine("plus", id)
    if BIAS_REFERENCE == "pooled":
        if GENOTYPE_MODE:
            # every -DMS sample of the ID's own genotype (the other genotype's BAMs hold other transcripts)
            return [s for s in get_samples_by_condition("minus") if SAMPLE_GENOTYPE[s] == ID_GENOTYPE[id]]
        return get_samples_by_condition("minus")
    return samples_to_combine("minus", id)


def get_bias_bams(wildcards):
    return [sample_bam(s) for s in get_bias_samples(wildcards.condition, wildcards.id)]


def get_bias_reference_minus_rtsc(wildcards):
    """-DMS RT-stops matching the reference depth, for the before/after plot."""
    if BIAS_REFERENCE == "pooled":
        if GENOTYPE_MODE:
            return genotype_pooled_minus_rtsc(ID_GENOTYPE[wildcards.id])
        return f"{TMP}/output/se/pooled/combined_minus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_minus.rtsc"


def get_bams_by_condition_and_id(condition, id):
    """Sorted BAMs of the samples in {id}/combined_<condition>.rtsc."""
    return [sample_bam(s) for s in samples_to_combine(condition, id)]


def _specificity_comparison_inputs(wildcards):
    runs = config.get("specificity_comparison_runs", [])
    inputs = []
    for run in runs:
        sheet = pd.read_csv(run["samplesheet"], sep="\t", dtype=str)
        for sample in dict.fromkeys(sheet["sample"].str.strip()):
            inputs.append(f"{run['output_dir']}/specificity_{sample}.csv")
    return inputs
