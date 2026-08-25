############################################################################
# StructureFold2-specific rules for RT stop counting and reactivity calculation
############################################################################


def get_rtsc_counts_by_condition(condition):
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in get_samples_by_condition(condition)
    ]


def get_rtsc_counts_by_condition_and_id(condition, id):
    rows = samples.loc[(samples["condition"] == condition) & (samples["ID"] == id)]
    if rows.empty:
        raise ValueError(f"No samples found for condition '{condition}' and ID '{id}'")
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in rows["sample"].tolist()
    ]


def get_rtsc_counts_by_condition_and_temperature(condition, temperature):
    return [
        f"{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc"
        for sample in get_samples_by_condition_and_temperature(condition, temperature)
    ]


# ---------------------------------------------------------------------------
# Replicate pooling before the reactivity calculation
#
# config['pool_replicates'] controls which channel(s), if any, get pooled
# across ALL replicates (all IDs) before rtsc_to_react runs, instead of each
# replicate keeping its own -DMS and +DMS RT-stop counts (today's default,
# "none"):
#
#   "none"  -- each replicate keeps its own +DMS and -DMS (default).
#   "minus" -- pool -DMS across every replicate into one shared background;
#              +DMS stays per-replicate. Still one reactivity output per ID.
#   "plus"  -- pool +DMS across every replicate into one shared treated
#              signal; -DMS stays per-replicate. Still one reactivity
#              output per ID.
#   "both"  -- pool BOTH channels across every replicate. There is only one
#              resulting reactivity output in this case (not one per
#              replicate, since both inputs are already collapsed to a
#              single pair) -- see IDS in the main Snakefile, which becomes
#              the single synthetic ID "pooled" for this mode. Every
#              per-ID rule elsewhere in the pipeline (coverage, specificity,
#              QC plots...) then runs on that single pooled dataset for
#              free, with no rule changes needed. Requires the samplesheet
#              to NOT have a 'temperature' column (see HEAT_CORRECTION in
#              the Snakefile) -- pooling across temperature would destroy
#              the distinction being corrected for.
#   "both_by_temperature" -- like "both", but pools +DMS/-DMS within each
#              temperature group separately instead of across the whole
#              sheet. Requires a 'temperature' samplesheet column. Produces
#              one synthetic ID per temperature ("pooled_<temperature>"),
#              so heat correction still has two (single-replicate) groups
#              to correct against each other.
#   "minus_all_plus_by_temperature" -- pool -DMS across EVERY replicate,
#              temperature included, into one shared background (like
#              "minus"), while +DMS is pooled only WITHIN each temperature
#              group separately (like the +DMS half of "both_by_temperature").
#              Useful when the -DMS (untreated) signal is assumed stable
#              across temperature and pooling it more broadly reduces noise,
#              while +DMS still needs to stay temperature-resolved for heat
#              correction. Requires a 'temperature' samplesheet column;
#              produces one synthetic ID per temperature ("pooled_<temperature>"),
#              same as "both_by_temperature".
#
# QC rules (coverage, specificity, nucleotide coverage, ...) report on
# whatever granularity IDS ends up at, so under "both"/"both_by_temperature"/
# "minus_all_plus_by_temperature" they report on the pooled dataset -- there
# is no finer-grained +DMS data left for them to report on at that point.
#
# Replicate correlation (rtsc_stop_correlation/_minus below) is the one
# exception: it always uses REPLICATE_IDS (the Snakefile's true
# biological-replicate list), never IDS, regardless of pool_replicates --
# comparing replicates against each other is meaningless once they've
# already been merged, so this check has to run on the still-unpooled
# per-replicate combine_rtsc_plus/minus output.
# ---------------------------------------------------------------------------
_ALLOWED_POOL_REPLICATES = {
    "none",
    "minus",
    "plus",
    "both",
    "both_by_temperature",
    "minus_all_plus_by_temperature",
}
if POOL_REPLICATES not in _ALLOWED_POOL_REPLICATES:
    raise ValueError(
        f"config['pool_replicates'] must be one of {sorted(_ALLOWED_POOL_REPLICATES)}, "
        f"got {POOL_REPLICATES!r}"
    )


def get_minus_rtsc_for_reactivity(wildcards):
    if POOL_REPLICATES in ("minus", "both", "minus_all_plus_by_temperature"):
        return f"{TMP}/output/se/pooled/combined_minus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_minus.rtsc"


def get_plus_rtsc_for_reactivity(wildcards):
    if POOL_REPLICATES in ("plus", "both"):
        return f"{TMP}/output/se/pooled/combined_plus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_plus.rtsc"


# ---------------------------------------------------------------------------
# reactivity_region (see the Snakefile comment near REACTIVITY_REGION):
# restricts ONLY the main reactivity.react/.csv calculation to one mRNA
# region. ACTIVE_TRANSCRIPTOME/get_*_rtsc_active are used exclusively by
# rtsc_to_react and convert_react_to_csv below -- every other rule
# (rtsc_to_react_plus_only, raw_reactivity, heat correction, QC) keeps using
# the full-transcript TRANSCRIPTOME_UPPER/get_*_rtsc_for_reactivity, since
# only the flagship reactivity output is in scope for this option.
# ---------------------------------------------------------------------------
ACTIVE_TRANSCRIPTOME = (
    f"{TMP}/resources/transcriptome_region.fa" if REACTIVITY_REGION else TRANSCRIPTOME_UPPER
)


def get_minus_rtsc_active(wildcards):
    path = get_minus_rtsc_for_reactivity(wildcards)
    return path.replace(".rtsc", ".region.rtsc") if REACTIVITY_REGION else path


def get_plus_rtsc_active(wildcards):
    path = get_plus_rtsc_for_reactivity(wildcards)
    return path.replace(".rtsc", ".region.rtsc") if REACTIVITY_REGION else path


if REACTIVITY_REGION:

    rule region_rtsc:
        """
        Slice a combined RT-stop-count file down to just reactivity_region,
        renumbering positions to start at 1 within the region -- consumed by
        rtsc_to_react (via get_*_rtsc_active) instead of the full-transcript
        combined_plus/minus.rtsc. Matches whichever rule produced the
        combined_{condition}.rtsc for this {id} (per-replicate, pooled, or
        pooled_by_temperature), since they all share this same output path
        shape.
        """
        input:
            rtsc=f"{TMP}/output/se/{{id}}/combined_{{condition}}.rtsc",
            coords=f"{config['output_dir']}/region_coordinates.tsv",
        output:
            f"{TMP}/output/se/{{id}}/combined_{{condition}}.region.rtsc",
        wildcard_constraints:
            condition="plus|minus",
        log:
            "logs/region_rtsc/{id}_{condition}.log",
        conda:
            "../envs/biopython.yaml"
        shell:
            """
            python3 workflow/scripts/extract_region_rtsc.py \
                --rtsc {input.rtsc} --coords {input.coords} --output {output} \
                2>{log}
            """


rule sam_to_rtsc:
    """
    Count RT-stops from the SAM file.

    Depends on TRANSCRIPTOME_UPPER (not just via params) so it's the first
    post-mapping rule to require it -- every downstream rule that reads
    TRANSCRIPTOME_UPPER only via params (not input:) relies on this edge to
    guarantee it already exists, since every later per-transcript file
    (.rtsc, .react, ...) traces back to this rule's output.
    """
    input:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
        tx_upper=TRANSCRIPTOME_UPPER,
    output:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    log:
        "logs/sam_to_rtsc/{sample}.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/sam_to_rtsc.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} -single {input.sam} {params.transcriptome} >{log} 2>&1
        """


rule combine_rtsc_plus:
    """
    Combine the RT-stop counts of the sample(s) making up one +DMS
    replicate (ID). wildcard_constraints excludes the literal ID "pooled"
    and anything starting with "pooled_" so this rule can never collide
    with combine_rtsc_plus_pooled or combine_rtsc_plus_pooled_by_temperature
    below.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("plus", wc.id),
    output:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus/{id}.log",
    wildcard_constraints:
        id="(?!pooled(_|$)).+",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus:
    """
    Combine the RT-stop counts of the sample(s) making up one -DMS
    replicate (ID). See combine_rtsc_plus for the wildcard_constraints note.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("minus", wc.id),
    output:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus/{id}.log",
    wildcard_constraints:
        id="(?!pooled(_|$)).+",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_plus_pooled:
    """
    Pool the +DMS RT-stop counts across ALL replicates into a single shared
    file. Always defined (harmless if unused); consumed by
    get_plus_rtsc_for_reactivity when pool_replicates is "plus" or "both",
    and directly by every other per-ID rule when pool_replicates == "both"
    (IDS becomes ["pooled"] in that mode -- see the main Snakefile).
    """
    input:
        get_rtsc_counts_by_condition("plus"),
    output:
        f"{TMP}/output/se/pooled/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus_pooled.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus_pooled:
    """
    Pool the -DMS RT-stop counts across ALL replicates into a single shared
    background file. See combine_rtsc_plus_pooled.
    """
    input:
        get_rtsc_counts_by_condition("minus"),
    output:
        f"{TMP}/output/se/pooled/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus_pooled.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_plus_pooled_by_temperature:
    """
    Pool the +DMS RT-stop counts across every replicate within one
    temperature group (pool_replicates: both_by_temperature). Always
    defined (harmless if unused). Its output lands at the same
    {TMP}/output/se/{id}/combined_plus.rtsc path every other per-ID rule
    already expects, with id == "pooled_<temperature>" -- see IDS in the
    main Snakefile under this mode -- so no other rule needs to change.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_temperature("plus", wc.temperature),
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus_pooled_by_temperature/{temperature}.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus_pooled_by_temperature:
    """
    Pool the -DMS RT-stop counts across every replicate within one
    temperature group. See combine_rtsc_plus_pooled_by_temperature.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_temperature("minus", wc.temperature),
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus_pooled_by_temperature/{temperature}.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule rtsc_to_react:
    """
    Convert the RT-stop counts to reactivities. When reactivity_region is
    set, plus/minus are region-restricted (get_*_rtsc_active) and
    transcriptome is the matching region-sliced fasta (ACTIVE_TRANSCRIPTOME)
    -- see the reactivity_region comment above get_plus_rtsc_active -- so
    both the 2-8% scale and the reported reactivities cover only that
    region.

    rtsc_to_react.py's own -trim3 already excludes the trimmed 3' positions
    from 2-8% scale generation, but still reports a real computed value for
    them -- mask_trim3_react.py re-NAs those same positions in the output so
    they're excluded from the reported reactivity too, and (via that NA)
    from anything downstream that sums/reads non-NA values, e.g. heat
    correction's scaling-factor calculation. trim3=0 is a no-op.
    """
    input:
        plus=get_plus_rtsc_active,
        minus=get_minus_rtsc_active,
        transcriptome=ACTIVE_TRANSCRIPTOME,
    output:
        f"{config['output_dir']}/{{id}}/reactivity.react",
    log:
        "logs/rtsc_to_react/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_to_react.py",
        output_prefix=f"{config['output_dir']}/{{id}}/reactivity",
        trim3=TRIM3,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python {params.workdir}/{params.script} {input.minus} {input.plus} {input.transcriptome} \
                -name {params.output_prefix} -trim3 {params.trim3} >{log} 2>&1 \
            && python3 workflow/scripts/mask_trim3_react.py \
                --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
            && mv {output}.masked {output}
        """


rule rtsc_to_react_plus_only:
    """
    Convert RT-stop counts to reactivities from the +DMS channel alone, with
    no -DMS (background) subtraction. Same formula as rtsc_to_react
    otherwise -- natural log, sum-normalize-by-length, 2-8% percentile
    scale, threshold cap -- just applied to one channel instead of
    subtracting a background channel from it. Same trim3 tail-masking as
    rtsc_to_react -- see that rule's docstring.
    """
    input:
        plus=get_plus_rtsc_for_reactivity,
    output:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.react",
    log:
        "logs/rtsc_to_react_plus_only/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/reactivity_plus_only",
        trim3=TRIM3,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/rtsc_to_react_plus_only.py {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -trim3 {params.trim3} >{log} 2>&1 \
            && python3 workflow/scripts/mask_trim3_react.py \
                --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
            && mv {output}.masked {output}
        """


rule rtsc_to_raw_react_plus_only:
    """
    Convert RT-stop counts to raw (un-2-8%-scaled, uncapped) +DMS-only
    reactivities (no -DMS subtraction) for one replicate/ID.
    """
    input:
        plus=get_plus_rtsc_for_reactivity,
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.react",
    log:
        "logs/rtsc_to_raw_react_plus_only/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only",
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/rtsc_to_react_plus_only.py \
                {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -nrm_off -threshold 1000000 >{log} 2>&1
        """


rule rtsc_to_raw_react:
    """
    Convert RT-stop counts to raw (un-2-8%-scaled, uncapped) reactivities
    for one replicate/ID. -nrm_off skips the per-replicate self-scale (but
    keeps its transcript QC filter); the very high -threshold keeps values
    uncapped, useful for feeding an external pooled-normalization step.
    """
    input:
        plus=get_plus_rtsc_for_reactivity,
        minus=get_minus_rtsc_for_reactivity,
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react",
    log:
        "logs/rtsc_to_raw_react/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity",
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -nrm_off -threshold 1000000 >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Optional: heat correction between two temperature conditions
#
# react_heat_correct.py rescales two sets of reactivity.react files (e.g. the
# same construct/transcriptome probed at a lower and a higher temperature) so
# both conditions share the same overall mean signal -- correcting for heat
# itself increasing overall reagent reactivity, independent of any real
# structural difference. It requires every input .react to share the exact
# same transcript set, AND its own docstring requires those transcripts to
# already be restricted to ones above an accepted coverage threshold -- a
# transcript with little to no +DMS signal in one temperature condition is
# not a meaningful member of either condition's shared signal scale.
#
# So before heat correction, restrict both conditions to transcripts with
# RT-stop coverage >= 1 (StructureFold2's standard AC-specific
# stops-per-base metric) in the POOLED +DMS samples of BOTH temperatures --
# pooled regardless of config['pool_replicates'], so a low-coverage
# replicate can't be masked by pooling, and the same qualifying set is used
# whether or not replicates are pooled for the reactivity calculation
# itself. The shared/overlapping transcript list itself is generated with
# coverage_overlap.py (workflow/scripts/StructureFold2/deprecated/) --
# upstream StructureFold2 has since folded that functionality into
# rtsc_coverage.py's own -ol flag, but this project calls the original
# script directly.
#
# Following Su et al. 2018 PNAS SI (Materials and Methods, "Determination
# of DMS reactivity", steps 3a/3b): the per-transcript 2-8% normalization
# scale is generated ONCE from the pooled LOWER-temperature data only, then
# that SAME scale is applied to the higher-temperature data too, rather
# than each temperature computing its own independent scale. DMS is
# intrinsically more reactive at higher temperature; letting the higher
# temperature renormalize itself would partly cancel that increase back
# out during normalization, masking the very heat-induced reactivity
# change react_heat_correct.py (step 4 in the SI) is meant to reveal. See
# heat_correction_lower_temperature_scale / rtsc_to_react_heat_shared below.
#
# Enabled by an optional 'temperature' column in the samplesheet -- see the
# Snakefile's samplesheet comment for the exact semantics. HEAT_CORRECTION,
# HEAT_CORRECTION_LOWER_IDS/HIGHER_IDS, and HEAT_CORRECTION_SUFFIX are
# computed there (from samples_runs) before this file is included.
# ---------------------------------------------------------------------------
if HEAT_CORRECTION:

    rule heat_correction_plus_coverage:
        """
        Per-transcript RT-stop coverage (rtsc_coverage.py, AC-specific
        stops-per-base) for the pooled +DMS samples of both temperature
        conditions (combine_rtsc_plus_pooled_by_temperature pools +DMS
        across every replicate of one temperature regardless of
        pool_replicates). Feeds heat_correction_shared_transcripts below.

        Each temperature's combined_plus.rtsc is copied to a bare
        <temperature>.rtsc name in a scratch dir first, then rtsc_coverage.py
        is run from there, so the resulting coverage.csv's column headers
        are plain temperature labels rather than full paths -- required by
        coverage_overlap.py below, which derives its own output filename
        from those headers and can't handle a '/' in one.
        """
        input:
            lower=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[0]}/combined_plus.rtsc",
            higher=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[1]}/combined_plus.rtsc",
        output:
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv",
        log:
            "logs/heat_correction_plus_coverage.log",
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/rtsc_coverage.py",
            transcriptome=os.path.abspath(TRANSCRIPTOME_UPPER),
            tmpdir=f"{TMP}/structurefold2/heat_correction_plus_coverage",
            lower_name=f"{DISTINCT_TEMPERATURES[0]}.rtsc",
            higher_name=f"{DISTINCT_TEMPERATURES[1]}.rtsc",
        shell:
            """
            mkdir -p {params.tmpdir} $(dirname {output}) $(dirname {log}) \
                && cp {input.lower} {params.tmpdir}/{params.lower_name} \
                && cp {input.higher} {params.tmpdir}/{params.higher_name} \
                && OUT_ABS=$(readlink -f {output}) \
                && LOG_ABS=$(readlink -f {log}) \
                && cd {params.tmpdir} \
                && python {params.workdir}/{params.script} \
                    {params.transcriptome} -f {params.lower_name} {params.higher_name} \
                    -name "$OUT_ABS" >"$LOG_ABS" 2>&1 \
                && cd - >/dev/null \
                && rm -rf {params.tmpdir}
            """

    rule heat_correction_shared_transcripts:
        """
        List transcripts with RT-stop coverage >= 1 in the pooled +DMS
        samples of BOTH temperature conditions (heat_correction_plus_coverage),
        using StructureFold2's coverage_overlap.py, for restricting the
        reactivity calculation feeding heat correction to a shared,
        coverage-qualified transcript set.

        coverage_overlap.py names its own output from the input CSV's
        column headers and always writes into the current directory, so
        this runs from a scratch dir and the result is moved to the
        declared output.
        """
        input:
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv",
        output:
            f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
        log:
            "logs/heat_correction_shared_transcripts.log",
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/deprecated/coverage_overlap.py",
            tmpdir=f"{TMP}/structurefold2/heat_correction_shared_transcripts",
            threshold=1.0,
        shell:
            """
            mkdir -p {params.tmpdir} $(dirname {output}) $(dirname {log}) \
                && IN_ABS=$(readlink -f {input}) \
                && LOG_ABS=$(readlink -f {log}) \
                && cd {params.tmpdir} \
                && python {params.workdir}/{params.script} -f "$IN_ABS" -n {params.threshold} >"$LOG_ABS" 2>&1 \
                && cd - >/dev/null \
                && mv {params.tmpdir}/*_overlap_{params.threshold}.txt {output} \
                && rm -rf {params.tmpdir}
            """

    rule heat_correction_lower_temperature_scale:
        """
        Generate the per-transcript 2-8% normalization scale (Su et al. 2018
        PNAS SI step 3a) from the FULLY POOLED lower-temperature dataset
        only (combine_rtsc_plus/minus_pooled_by_temperature -- defined
        regardless of pool_replicates), restricted to the coverage-qualified
        shared transcripts. rtsc_to_react_heat_shared below applies this
        SAME scale to every ID on both sides of the correction (step 3b),
        rather than each ID computing its own independent scale, matching
        the paper's method. Using the fully pooled lower-temperature data
        (rather than each lower replicate's own noisier scale) also gives a
        single, well-defined scale to share when pool_replicates isn't
        both_by_temperature and there are multiple lower/higher replicates.
        """
        input:
            plus=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[0]}/combined_plus.rtsc",
            minus=f"{TMP}/output/se/pooled_{DISTINCT_TEMPERATURES[0]}/combined_minus.rtsc",
            restrict=f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
        output:
            react=f"{config['output_dir']}/qc/heat_correction_lower_temperature.react",
            scale=f"{config['output_dir']}/qc/heat_correction_lower_temperature.scale",
        log:
            "logs/heat_correction_lower_temperature_scale.log",
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/qc/heat_correction_lower_temperature",
            trim3=TRIM3,
        shell:
            """
            mkdir -p $(dirname {output.react}) \
                && python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                    -name {params.output_prefix} -restrict {input.restrict} -trim3 {params.trim3} >{log} 2>&1 \
                && python3 workflow/scripts/mask_trim3_react.py \
                    --input {output.react} --output {output.react}.masked --trim3 {params.trim3} >>{log} 2>&1 \
                && mv {output.react}.masked {output.react}
            """

    rule rtsc_to_react_heat_shared:
        """
        Same reactivity calculation as rtsc_to_react, restricted (-restrict)
        to the transcripts in heat_correction_shared_transcripts, and
        normalized with the shared lower-temperature 2-8% scale
        (heat_correction_lower_temperature_scale) instead of each ID
        computing its own -- see the "Optional: heat correction" comment
        block above for why. Coverage alone doesn't guarantee an identical
        transcript set across conditions though -- rtsc_to_react.py also
        independently drops transcripts with zero +DMS signal or an
        unresolvable normalization scale -- so this is an intermediate
        file, further narrowed by heat_correction_shared_react below.

        A -scale is supplied here, so rtsc_to_react.py's own -trim3 flag
        would be a no-op (it only affects scale *generation*, skipped
        entirely when -scale is given) -- the mask_trim3_react.py step
        below is what actually excludes the trim3'd tail from this file's
        reported values, using config's trim3 directly rather than relying
        on the python script's -trim3 flag. Without this, those tail
        positions would carry real values into heat_correct_reactivity's
        scaling-factor calculation (react_heat_correct.py's sum_react)
        even though trim3 already excluded them from the scale that
        produced heat_correction_lower_temperature.scale above.
        """
        input:
            plus=get_plus_rtsc_for_reactivity,
            minus=get_minus_rtsc_for_reactivity,
            restrict=f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
            scale=f"{config['output_dir']}/qc/heat_correction_lower_temperature.scale",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection.react",
        log:
            "logs/rtsc_to_react_heat_shared/{id}.log",
        wildcard_constraints:
            id="[^/]+",
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection",
            trim3=TRIM3,
        shell:
            """
            mkdir -p $(dirname {output}) \
                && python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                    -name {params.output_prefix} -restrict {input.restrict} -scale {input.scale} >{log} 2>&1 \
                && python3 workflow/scripts/mask_trim3_react.py \
                    --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
                && mv {output}.masked {output}
            """

    rule heat_correction_shared_react:
        """
        Narrow every rtsc_to_react_heat_shared output (both temperatures,
        all replicates) down to the transcripts present in ALL of them, so
        react_heat_correct.py's "exact same transcripts" requirement is
        actually satisfied -- the coverage restriction alone isn't enough,
        see rtsc_to_react_heat_shared.
        """
        input:
            expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection.react",
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            ),
        output:
            expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            ),
        log:
            "logs/heat_correction_shared_react.log",
        conda:
            "../envs/structurefold.yaml"
        params:
            script=f"{workflow.basedir}/workflow/scripts/react_intersect_transcripts.py",
        shell:
            """
            python {params.script} -in {input} -out {output} >{log} 2>&1
            """

    rule convert_pre_heat_correction_react_to_csv:
        """
        Convert the pre-correction reactivities (heat_correction_shared_react
        -- coverage- and transcript-set-restricted, but not yet rescaled) to
        CSV form (transcript, position, base, reactivity -- one row per
        position), same as convert_react_to_csv. Lets the corrected and
        uncorrected values be compared directly, since both sides cover the
        exact same transcripts/positions.
        """
        input:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.csv",
        log:
            "logs/convert_pre_heat_correction_react_to_csv/{id}.log",
        wildcard_constraints:
            id="[^/]+",
        conda:
            "../envs/biopython.yaml"
        params:
            transcriptome=TRANSCRIPTOME_UPPER,
        shell:
            """
            python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
            """

    rule heat_correct_reactivity:
        """
        Scale coverage- and transcript-set-restricted reactivity.react
        files (heat_correction_shared_react) from a lower- and a
        higher-temperature condition onto a common overall-signal scale,
        supporting multiple replicates per temperature. Uses
        react_heat_correct.py by default (correction factors computed from
        every coverage-qualified position), or
        react_heat_correct_positive_only.py when
        config['heat_correction_positive_bases_only'] is set -- same
        correction mechanics, but the factors are computed only from
        positions with reactivity > 0 at BOTH temperatures (excluding exact
        0.0s, not just NA). Either way the correction itself is still
        applied to every position.
        """
        input:
            lower=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_LOWER_IDS,
            ),
            higher=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_HIGHER_IDS,
            ),
        output:
            lower=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                id=HEAT_CORRECTION_LOWER_IDS,
            ),
            higher=expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                id=HEAT_CORRECTION_HIGHER_IDS,
            ),
        log:
            f"{config['output_dir']}/qc/heat_correction_scale_factors.log",
        conda:
            "../envs/structurefold.yaml"
        params:
            script=(
                f"{workflow.basedir}/workflow/scripts/react_heat_correct_positive_only.py"
                if HEAT_CORRECTION_POSITIVE_BASES_ONLY
                else f"{workflow.basedir}/workflow/scripts/StructureFold2/react_heat_correct.py"
            ),
            suffix=HEAT_CORRECTION_SUFFIX,
        shell:
            """
            python {params.script} \
                -lower {input.lower} \
                -higher {input.higher} \
                -suffix {params.suffix} >{log} 2>&1
            """

    rule convert_heat_corrected_react_to_csv:
        """
        Convert heat-corrected reactivities to CSV form (transcript,
        position, base, reactivity -- one row per position), same as
        convert_react_to_csv.
        """
        input:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.csv",
        log:
            "logs/convert_heat_corrected_react_to_csv/{id}.log",
        wildcard_constraints:
            id="[^/]+",
        conda:
            "../envs/biopython.yaml"
        params:
            transcriptome=TRANSCRIPTOME_UPPER,
        shell:
            """
            python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
            """

    # -----------------------------------------------------------------------
    # Optional comparison run: heat-correct on EVERY transcript present in
    # both temperature conditions' PLAIN reactivity.react (rtsc_to_react),
    # with no RT-stop coverage threshold at all -- to see how much the
    # coverage-qualification step above (heat_correction_shared_transcripts)
    # actually changes the correction, vs. just intersecting whatever
    # rtsc_to_react.py could compute a reactivity for. Mirrors
    # heat_correction_shared_react / heat_correct_reactivity /
    # convert_heat_corrected_react_to_csv exactly, minus the coverage
    # restriction. Enabled by config['heat_correction_compare_all_transcripts'].
    # -----------------------------------------------------------------------
    if HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS:

        rule heat_correction_all_transcripts_shared_react:
            """
            Narrow the plain, unrestricted reactivity.react (rtsc_to_react)
            for every heat-correction ID down to the transcripts present in
            ALL of them -- no coverage threshold, just whatever
            rtsc_to_react.py itself could compute a reactivity for. The
            no-coverage-filter counterpart to heat_correction_shared_react.
            """
            input:
                expand(
                    f"{config['output_dir']}/{{id}}/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                ),
            output:
                expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                ),
            log:
                "logs/heat_correction_all_transcripts_shared_react.log",
            conda:
                "../envs/structurefold.yaml"
            params:
                script=f"{workflow.basedir}/workflow/scripts/react_intersect_transcripts.py",
            shell:
                """
                mkdir -p $(dirname {output[0]}) \
                    && python {params.script} -in {input} -out {output} >{log} 2>&1
                """

        rule heat_correct_reactivity_all_transcripts:
            """
            The no-coverage-filter counterpart to heat_correct_reactivity --
            same react_heat_correct.py rescaling, run on every transcript
            present in both conditions regardless of coverage
            (heat_correction_all_transcripts_shared_react).
            """
            input:
                lower=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS,
                ),
                higher=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_HIGHER_IDS,
                ),
            output:
                lower=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                    id=HEAT_CORRECTION_LOWER_IDS,
                ),
                higher=expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
                    id=HEAT_CORRECTION_HIGHER_IDS,
                ),
            log:
                f"{config['output_dir']}/qc/heat_correction_all_transcripts_scale_factors.log",
            conda:
                "../envs/structurefold.yaml"
            params:
                workdir=f"{workflow.basedir}/workflow",
                script="scripts/StructureFold2/react_heat_correct.py",
                suffix=HEAT_CORRECTION_SUFFIX,
            shell:
                """
                python {params.workdir}/{params.script} \
                    -lower {input.lower} \
                    -higher {input.higher} \
                    -suffix {params.suffix} >{log} 2>&1
                """

        rule convert_heat_correction_all_transcripts_to_csv:
            """
            Convert the all-transcripts (no coverage filter) heat-corrected
            reactivities to CSV form, same as convert_heat_corrected_react_to_csv.
            """
            input:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react",
            output:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.csv",
            log:
                "logs/convert_heat_correction_all_transcripts_to_csv/{id}.log",
            wildcard_constraints:
                id="[^/]+",
            conda:
                "../envs/biopython.yaml"
            params:
                transcriptome=TRANSCRIPTOME_UPPER,
            shell:
                """
                python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
                """


# ---------------------------------------------------------------------------
# Normalization comparison: all 8 combinations of trim / ln / nrm options
# (experimental, on-demand only -- not wired to pool_replicates; always uses
# each replicate's own per-ID rtsc regardless of the pooling config)
# ---------------------------------------------------------------------------
# Wildcard encoding:
#   {trim}  = "trim50"  -> -trim3 50      | "notrim"  -> (no flag, default 0)
#   {nlog}  = "noln"    -> -ln_off        | "ln"      -> (no flag, ln is ON)
#   {norm}  = "nonrm"   -> -nrm_off       | "nrm"     -> (no flag, nrm is ON)

NORM_COMBOS = expand(
    "{trim}_{nlog}_{norm}",
    trim=["notrim", "trim125"],
    nlog=["noln", "ln"],
    norm=["nonrm", "nrm"],
)


rule rtsc_to_react_comparison:
    input:
        plus=f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
        minus=f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react",
    log:
        "logs/rtsc_to_react_comparison/{id}_{trim}_{nlog}_{norm}.log",
    wildcard_constraints:
        trim="notrim|trim125",
        nlog="noln|ln",
        norm="nonrm|nrm",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=lambda wc: f"{config['output_dir']}/{wc.id}/norm_comparison/{wc.trim}_{wc.nlog}_{wc.norm}/reactivity",
        trim_flag=lambda wc: "-trim3 125" if wc.trim == "trim125" else "",
        trim3=lambda wc: 125 if wc.trim == "trim125" else 0,
        ln_flag=lambda wc: "-ln_off" if wc.nlog == "noln" else "",
        nrm_flag=lambda wc: "-nrm_off" if wc.norm == "nonrm" else "",
    shell:
        """
        python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
            -name {params.output_prefix} \
            {params.trim_flag} {params.ln_flag} {params.nrm_flag} >{log} 2>&1 \
        && python3 workflow/scripts/mask_trim3_react.py \
            --input {output} --output {output}.masked --trim3 {params.trim3} >>{log} 2>&1 \
        && mv {output}.masked {output}
        """


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


rule samtools_depth_by_condition:
    """Run samtools depth (all positions) on SAM files for one condition."""
    input:
        lambda wc: get_sams_by_condition_and_id(wc.condition, wc.id),
    output:
        f"{config['output_dir']}/{{id}}/depth_{{condition}}.txt",
    log:
        "logs/samtools_depth_by_condition/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
        condition="plus|minus",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && samtools depth -a {input} >{output} 2>{log}
        """


rule calculate_stop_coverage:
    input:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
    output:
        coverage=f"{config['output_dir']}/{{id}}/coverage.csv",
    log:
        "logs/calculate_stop_coverage/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        coverage_name=f"{config['output_dir']}/{{id}}/coverage",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_coverage.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        mkdir -p $(dirname {output.coverage}) \
            && python {params.workdir}/{params.script} -f {input} -name {params.coverage_name} {params.transcriptome} >{log} 2>&1
        """


rule rtsc_total_stops:
    """
    Sum RT-stop counts per transcript from the untreated (-DMS) channel to
    use as a raw read-count proxy for differential expression (DESeq2-style
    tools expect raw, un-normalized counts and do their own library-size
    normalization).
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/counts_minus.csv",
    log:
        "logs/rtsc_total_stops/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/upset.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/rtsc_total_stops.py --rtsc {input} --output {output} >{log} 2>&1
        """


rule calculate_transcript_abundance:
    """
    Estimate per-transcript relative abundance (RPKM or TPM) from the
    untreated (-DMS) RT-stop counts, using StructureFold2's own
    rtsc_abundances.py. Same -DMS-channel rationale as rtsc_total_stops
    above (unstructured RT stops as a raw-count proxy for expression), but
    here converted to a length- and library-size-normalized abundance
    metric instead of DESeq2-style raw counts. Enabled per-mode via
    config['transcript_abundance'] (see the main Snakefile).

    rtsc_abundances.py always names its output by replacing the input
    file's own '.rtsc' extension (e.g. combined_minus_RPKM.csv next to
    combined_minus.rtsc) and can't be told an output path directly, so this
    runs from a scratch dir on a bare-named copy of the input, then moves
    the result to the declared output -- same pattern as combine_rtsc_plus.
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/abundance_{{mode}}.csv",
    log:
        "logs/calculate_transcript_abundance/{id}_{mode}.log",
    wildcard_constraints:
        id="[^/]+",
        mode="RPKM|TPM",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_abundances.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_abundance_{{mode}}",
    shell:
        """
        (
            mkdir -p {params.tmpdir} $(dirname {output}) \
                && cp {input} {params.tmpdir}/combined_minus.rtsc \
                && cd {params.tmpdir} \
                && python {params.workdir}/{params.script} {wildcards.mode} -f combined_minus.rtsc \
                && cd - >/dev/null \
                && mv {params.tmpdir}/combined_minus_{wildcards.mode}.csv {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule calculate_nucleotide_coverage:
    input:
        rtsc=f"{TMP}/output/se/{{id}}/combined_{{condition}}.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/nucleotide_coverage_{{condition}}.csv",
    log:
        "logs/calculate_nucleotide_coverage/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/rtsc_to_nucleotide_coverage.py \
            --rtsc {input.rtsc} \
            --fasta {params.transcriptome} \
            --output {output} >{log} 2>&1
        """


rule calculate_specificity:
    """
    Calculate the A/C specificity for the RT-stops
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_{{treatment}}.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/specificity_{{treatment}}.csv",
    log:
        "logs/calculate_specificity/{id}_{treatment}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        output=f"{config['output_dir']}/{{id}}/specificity_{{treatment}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome} >{log} 2>&1
        """


rule convert_react_to_csv:
    """
    Convert the reactivities to CSV form. Uses ACTIVE_TRANSCRIPTOME (the
    same region-sliced fasta rtsc_to_react used, when reactivity_region is
    set) since react_to_csv.py looks up each position's base by indexing
    into this fasta and skips any transcript whose length doesn't match the
    .react file.
    """
    input:
        react=f"{config['output_dir']}/{{id}}/reactivity.react",
        transcriptome=ACTIVE_TRANSCRIPTOME,
    output:
        f"{config['output_dir']}/{{id}}/reactivity.csv",
    log:
        "logs/convert_react_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input.react} --output {output} --fasta {input.transcriptome} >{log} 2>&1
        """


rule convert_react_plus_only_to_csv:
    """
    Convert the +DMS-only (no -DMS subtraction) reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.react",
    output:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.csv",
    log:
        "logs/convert_react_plus_only_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """


rule convert_raw_react_plus_only_to_csv:
    """
    Convert raw (un-2-8%-scaled) +DMS-only reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.react",
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.csv",
    log:
        "logs/convert_raw_react_plus_only_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """


rule convert_raw_react_to_csv:
    """
    Convert raw (un-2-8%-scaled) reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react",
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.csv",
    log:
        "logs/convert_raw_react_to_csv/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """


rule convert_react_comparison_to_csv:
    input:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react",
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.csv",
    log:
        "logs/convert_react_comparison_to_csv/{id}_{trim}_{nlog}_{norm}.log",
    wildcard_constraints:
        trim="notrim|trim125",
        nlog="noln|ln",
        norm="nonrm|nrm",
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome} >{log} 2>&1
        """


rule calculate_specificity_by_sample:
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    output:
        f"{config['output_dir']}/specificity_{{sample}}.csv",
    log:
        "logs/calculate_specificity_by_sample/{sample}.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        output=f"{config['output_dir']}/specificity_{{sample}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome} >{log} 2>&1
        """


rule filter_covered_transcripts:
    """
    Intersect coverage CSVs across all replicates/IDs and write a plain-text
    list of transcript names with mean RT-stop coverage >= 1 in every
    replicate. Used as the -restrict input to rtsc_correlation.py.
    """
    input:
        expand(
            "{out}/{id}/coverage.csv",
            out=config["output_dir"],
            id=IDS,
        ),
    output:
        f"{config['output_dir']}/qc/covered_transcripts.txt",
    log:
        "logs/filter_covered_transcripts.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/filter_covered_transcripts.py \
                --coverages {input} \
                --threshold 1 \
                --output {output} >{log} 2>&1
        """


rule upset_covered_transcripts_replicates:
    """
    UpSet plot of transcripts with RT-stop coverage >= 1 (calculate_stop_coverage's
    AC-specific stops-per-base metric) in the +DMS RT-stops of each TRUE
    biological replicate (REPLICATE_IDS, not IDS) -- always the real,
    unpooled per-replicate coverage regardless of config['pool_replicates'],
    same rationale as rtsc_stop_correlation below (this comparison is
    meaningless once replicates are already merged). Visualizes every
    replicate-overlap combination, including the single "covered
    everywhere" bar that filter_covered_transcripts.py above already
    reduces to a plain-text list.
    """
    input:
        expand(f"{config['output_dir']}/{{id}}/coverage.csv", id=REPLICATE_IDS),
    output:
        f"{config['output_dir']}/qc/covered_transcripts_upset_replicates.png",
    log:
        "logs/upset_covered_transcripts_replicates.log",
    conda:
        "../envs/upset.yaml"
    params:
        labels=REPLICATE_IDS,
        threshold=1.0,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/plot_covered_transcripts_upset.py \
                --coverages {input} \
                --labels {params.labels} \
                --threshold {params.threshold} \
                --output {output} >{log} 2>&1
        """


if MERGED_IDS_DIFFER and len(IDS) > 1:

    rule upset_covered_transcripts_merged:
        """
        Same UpSet plot as upset_covered_transcripts_replicates, but over
        the pool_replicates-merged grouping (IDS) instead of the true
        replicate list -- e.g. one set per pooled_<temperature> under
        both_by_temperature/minus_all_plus_by_temperature. Only defined
        when pooling actually collapses IDS to something other than
        REPLICATE_IDS (MERGED_IDS_DIFFER, set in the main Snakefile) *and*
        that still leaves more than one set to compare -- "both" pooling
        collapses everything to the single ID "pooled", which has nothing
        left to overlap with and can't be UpSet-plotted.
        """
        input:
            expand(f"{config['output_dir']}/{{id}}/coverage.csv", id=IDS),
        output:
            f"{config['output_dir']}/qc/covered_transcripts_upset_merged.png",
        log:
            "logs/upset_covered_transcripts_merged.log",
        conda:
            "../envs/upset.yaml"
        params:
            labels=IDS,
            threshold=1.0,
        shell:
            """
            mkdir -p $(dirname {output}) \
                && python3 workflow/scripts/plot_covered_transcripts_upset.py \
                    --coverages {input} \
                    --labels {params.labels} \
                    --threshold {params.threshold} \
                    --output {output} >{log} 2>&1
            """


rule rtsc_stop_correlation:
    """
    Run StructureFold2 rtsc_correlation.py on the combined-plus RTSC files
    for every TRUE biological replicate (REPLICATE_IDS, not IDS) -- this
    check is meaningless on already-merged/pooled data, so it always uses
    each real replicate's own combine_rtsc_plus output regardless of
    config['pool_replicates'].
    """
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
            id=REPLICATE_IDS,
        ),
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv",
    log:
        "logs/rtsc_stop_correlation.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_correlation.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} \
            {input.rtsc} \
            -spec AC \
            -fasta {params.transcriptome} \
            -name {output} >{log} 2>&1
        """


rule plot_rt_stop_correlation:
    """
    Pairwise hexbin scatter of log10(RT-stop + 1) counts at A and C positions
    across all replicates, with Pearson r annotated on each panel.
    """
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv",
    output:
        report(
            f"{config['output_dir']}/qc/rt_stop_replicate_correlation.png",
            category="QC",
            subcategory="Replicate correlation",
            caption="../report_captions/rt_stop_replicate_correlation.rst",
        ),
    log:
        "logs/plot_rt_stop_correlation.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output} >{log} 2>&1
        """


rule rtsc_stop_correlation_minus:
    """
    Same as rtsc_stop_correlation, but for the combined-minus (untreated,
    -DMS/background) RTSC files -- lets replicate agreement in the
    background channel be checked the same way as the treated channel.
    Also uses REPLICATE_IDS (true biological replicates), not IDS, for the
    same reason as rtsc_stop_correlation.
    """
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
            id=REPLICATE_IDS,
        ),
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv",
    log:
        "logs/rtsc_stop_correlation_minus.log",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_correlation.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} \
            {input.rtsc} \
            -spec AC \
            -fasta {params.transcriptome} \
            -name {output} >{log} 2>&1
        """


rule plot_rt_stop_correlation_minus:
    """
    Same as plot_rt_stop_correlation, but for the combined-minus (untreated,
    -DMS/background) RT-stop counts.
    """
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv",
    output:
        report(
            f"{config['output_dir']}/qc/rt_stop_replicate_correlation_minus.png",
            category="QC",
            subcategory="Replicate correlation",
            caption="../report_captions/rt_stop_replicate_correlation_minus.rst",
        ),
    log:
        "logs/plot_rt_stop_correlation_minus.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output} >{log} 2>&1
        """


rule rtsc_end_coverage_sweep:
    """
    Run rtsc_end_coverage.py in TP mode on the combined minus RTSC for trim
    values 0..200, recording the mean TP coverage at each trim amount.
    Only non-NA transcripts (those long enough for the given trim + tp_l) are
    included in each mean. On-demand only -- not in the default target list.
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.tsv",
    log:
        "logs/rtsc_end_coverage_sweep/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_end_coverage.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_end_coverage_sweep",
        length=50,
        tp_l=300,
    shell:
        r"""
        mkdir -p {params.tmpdir}
        mkdir -p $(dirname {output})
        printf 'trim\tmean_tp_coverage\n' >{output}
        : >{log}
        for trim in $(seq 0 200); do
            tmpout={params.tmpdir}/trim_$trim.csv
            python {params.workdir}/{params.script} {input} TP \
                -length {params.length} -tp_l {params.tp_l} -trim $trim \
                -name $tmpout >>{log} 2>&1
            mean=$(awk -F',' \
                'NR>1 && $2 != "NA" {{sum+=$2; count++}} END {{if(count>0) printf "%.6f", sum/count; else print "NA"}}' \
                $tmpout)
            printf '%d\t%s\n' $trim $mean >>{output}
        done
        rm -rf {params.tmpdir}
        """


rule plot_end_coverage_sweep:
    """
    Plot mean 3' TP coverage vs. trim amount from the sweep TSV produced by
    rtsc_end_coverage_sweep. On-demand only.
    """
    input:
        f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.tsv",
    output:
        report(
            f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.png",
            category="QC",
            subcategory="End coverage sweep",
            labels={"id": "{id}"},
            caption="../report_captions/end_coverage_trim_sweep.rst",
        ),
    log:
        "logs/plot_end_coverage_sweep/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_end_coverage_sweep.py \
            --input {input} \
            --output {output} \
            --id {wildcards.id} >{log} 2>&1
        """


rule plot_specificity:
    input:
        expand(
            "{output_dir}/specificity_{sample}.csv",
            output_dir=config["output_dir"],
            sample=SAMPLES,
        ),
    output:
        report(
            f"{config['output_dir']}/qc/specificity_plot.png",
            category="QC",
            subcategory="Specificity",
            caption="../report_captions/specificity_plot.rst",
        ),
    log:
        "logs/plot_specificity.log",
    conda:
        "../envs/plotting.yaml"
    params:
        samples=SAMPLES,
        conditions=samples["condition"].tolist(),
    shell:
        """
        python3 workflow/scripts/plot_specificity.py \
            --inputs {input} \
            --samples {params.samples} \
            --conditions {params.conditions} \
            --output {output} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Optional: combined specificity figure across several named comparison runs
# (e.g. different organisms, genotypes, or conditions, each its own
# pipeline run with its own output_dir/samplesheet). Not organism-specific.
#
# Requires config key `specificity_comparison_runs`, e.g.:
#
#   specificity_comparison_runs:
#     - name: "Wild-type"
#       output_dir: results_wildtype
#       samplesheet: config/samples_wildtype.tsv
#     - name: "Mutant"
#       output_dir: results_mutant
#       samplesheet: config/samples_mutant.tsv
#
# Columns in the figure: comparison group x condition (+/- DMS)
# Rows in the figure:    replicate (from the ID column of each samplesheet)
# ---------------------------------------------------------------------------


def _specificity_comparison_inputs(wildcards):
    runs = config.get("specificity_comparison_runs", [])
    inputs = []
    for run in runs:
        sheet = pd.read_csv(run["samplesheet"], sep="\t")
        for sample in sheet["sample"].tolist():
            inputs.append(f"{run['output_dir']}/specificity_{sample}.csv")
    return inputs


rule plot_specificity_comparison:
    input:
        _specificity_comparison_inputs,
    output:
        plot=report(
            f"{config['output_dir']}/qc/specificity_comparison.png",
            category="QC",
            subcategory="Specificity",
            caption="../report_captions/specificity_comparison.rst",
        ),
        manifest=f"{config['output_dir']}/qc/specificity_comparison_manifest.tsv",
    run:
        import pandas as pd

        rows = []
        for run in config["specificity_comparison_runs"]:
            sheet = pd.read_csv(run["samplesheet"], sep="\t")
            for _, row in sheet.iterrows():
                rows.append(
                    {
                        "file": f"{run['output_dir']}/specificity_{row['sample']}.csv",
                        "group": run["name"],
                        "condition": row["condition"],
                        "replicate": row.get("ID", row["sample"]),
                    }
                )
        pd.DataFrame(rows).to_csv(output.manifest, sep="\t", index=False)
        shell(
            "Rscript workflow/scripts/plot_specificity_combined.R"
            " {output.manifest} {output.plot}"
        )
