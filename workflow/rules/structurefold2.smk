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
#
# QC rules (coverage, specificity, nucleotide coverage, replicate
# correlation, ...) always report on the true per-replicate data except
# under "both"/"both_by_temperature", where there is no finer-grained data
# left to report on.
# ---------------------------------------------------------------------------
_ALLOWED_POOL_REPLICATES = {"none", "minus", "plus", "both", "both_by_temperature"}
if POOL_REPLICATES not in _ALLOWED_POOL_REPLICATES:
    raise ValueError(
        f"config['pool_replicates'] must be one of {sorted(_ALLOWED_POOL_REPLICATES)}, "
        f"got {POOL_REPLICATES!r}"
    )


def get_minus_rtsc_for_reactivity(wildcards):
    if POOL_REPLICATES in ("minus", "both"):
        return f"{TMP}/output/se/pooled/combined_minus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_minus.rtsc"


def get_plus_rtsc_for_reactivity(wildcards):
    if POOL_REPLICATES in ("plus", "both"):
        return f"{TMP}/output/se/pooled/combined_plus.rtsc"
    return f"{TMP}/output/se/{wildcards.id}/combined_plus.rtsc"


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
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/sam_to_rtsc.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python {params.workdir}/{params.script} -single {input.sam} {params.transcriptome}
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
        lambda wc: get_rtsc_counts_by_condition_and_id("plus", wc.id)
    output:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc"
    wildcard_constraints:
        id="(?!pooled(_|$)).+"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_plus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} &&
        mv {params.tmpdir}/combined_plus.rtsc {output} && \
        rm -rf {params.tmpdir}
        """


rule combine_rtsc_minus:
    """
    Combine the RT-stop counts of the sample(s) making up one -DMS
    replicate (ID). See combine_rtsc_plus for the wildcard_constraints note.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("minus", wc.id)
    output:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc"
    wildcard_constraints:
        id="(?!pooled(_|$)).+"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_minus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} &&
        mv {params.tmpdir}/combined_minus.rtsc {output} && \
        rm -rf {params.tmpdir}
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
        get_rtsc_counts_by_condition("plus")
    output:
        f"{TMP}/output/se/pooled/combined_plus.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_plus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} &&
        mv {params.tmpdir}/combined_plus.rtsc {output} && \
        rm -rf {params.tmpdir}
        """


rule combine_rtsc_minus_pooled:
    """
    Pool the -DMS RT-stop counts across ALL replicates into a single shared
    background file. See combine_rtsc_plus_pooled.
    """
    input:
        get_rtsc_counts_by_condition("minus")
    output:
        f"{TMP}/output/se/pooled/combined_minus.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_minus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} &&
        mv {params.tmpdir}/combined_minus.rtsc {output} && \
        rm -rf {params.tmpdir}
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
        lambda wc: get_rtsc_counts_by_condition_and_temperature("plus", wc.temperature)
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_plus.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_plus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} &&
        mv {params.tmpdir}/combined_plus.rtsc {output} && \
        rm -rf {params.tmpdir}
        """


rule combine_rtsc_minus_pooled_by_temperature:
    """
    Pool the -DMS RT-stop counts across every replicate within one
    temperature group. See combine_rtsc_plus_pooled_by_temperature.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_temperature("minus", wc.temperature)
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_minus.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_minus"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        cp {input} {params.tmpdir}/ && \
        python {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} &&
        mv {params.tmpdir}/combined_minus.rtsc {output} && \
        rm -rf {params.tmpdir}
        """


rule rtsc_to_react:
    """
    Convert the RT-stop counts to reactivities
    """
    input:
        plus=get_plus_rtsc_for_reactivity,
        minus=get_minus_rtsc_for_reactivity
    output:
        f"{config['output_dir']}/{{id}}/reactivity.react"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/reactivity",
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} -name {params.output_prefix}
        """


rule rtsc_to_react_plus_only:
    """
    Convert RT-stop counts to reactivities from the +DMS channel alone, with
    no -DMS (background) subtraction. Same formula as rtsc_to_react
    otherwise -- natural log, sum-normalize-by-length, 2-8% percentile
    scale, threshold cap -- just applied to one channel instead of
    subtracting a background channel from it.
    """
    input:
        plus=get_plus_rtsc_for_reactivity
    output:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.react"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/reactivity_plus_only",
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python3 workflow/scripts/rtsc_to_react_plus_only.py {input.plus} {params.transcriptome} -name {params.output_prefix}
        """


rule rtsc_to_raw_react_plus_only:
    """
    Convert RT-stop counts to raw (un-2-8%-scaled, uncapped) +DMS-only
    reactivities (no -DMS subtraction) for one replicate/ID.
    """
    input:
        plus=get_plus_rtsc_for_reactivity
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.react"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only",
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python3 workflow/scripts/rtsc_to_react_plus_only.py \
            {input.plus} {params.transcriptome} \
            -name {params.output_prefix} -nrm_off -threshold 1000000
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
        minus=get_minus_rtsc_for_reactivity
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_to_react.py",
        transcriptome=TRANSCRIPTOME_UPPER,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity",
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
            -name {params.output_prefix} -nrm_off -threshold 1000000
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
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv"
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
            mkdir -p {params.tmpdir} $(dirname {output}) && \
            cp {input.lower} {params.tmpdir}/{params.lower_name} && \
            cp {input.higher} {params.tmpdir}/{params.higher_name} && \
            OUT_ABS=$(readlink -f {output}) && \
            cd {params.tmpdir} && \
            python {params.workdir}/{params.script} \
                {params.transcriptome} -f {params.lower_name} {params.higher_name} \
                -name "$OUT_ABS" && \
            cd - > /dev/null && \
            rm -rf {params.tmpdir}
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
            f"{config['output_dir']}/qc/heat_correction_plus_coverage.csv"
        output:
            f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt"
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/deprecated/coverage_overlap.py",
            tmpdir=f"{TMP}/structurefold2/heat_correction_shared_transcripts",
            threshold=1.0,
        shell:
            """
            mkdir -p {params.tmpdir} $(dirname {output}) && \
            IN_ABS=$(readlink -f {input}) && \
            cd {params.tmpdir} && \
            python {params.workdir}/{params.script} -f "$IN_ABS" -n {params.threshold} && \
            cd - > /dev/null && \
            mv {params.tmpdir}/*_overlap_{params.threshold}.txt {output} && \
            rm -rf {params.tmpdir}
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
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/qc/heat_correction_lower_temperature",
        shell:
            """
            mkdir -p $(dirname {output.react}) && \
            python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -restrict {input.restrict}
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
        """
        input:
            plus=get_plus_rtsc_for_reactivity,
            minus=get_minus_rtsc_for_reactivity,
            restrict=f"{config['output_dir']}/qc/heat_correction_shared_transcripts.txt",
            scale=f"{config['output_dir']}/qc/heat_correction_lower_temperature.scale",
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection.react"
        wildcard_constraints:
            id="[^/]+"
        conda:
            "../envs/structurefold.yaml"
        params:
            workdir=f"{workflow.basedir}/workflow",
            script="scripts/StructureFold2/rtsc_to_react.py",
            transcriptome=TRANSCRIPTOME_UPPER,
            output_prefix=f"{config['output_dir']}/{{id}}/heat_correction/reactivity_precorrection",
        shell:
            """
            mkdir -p $(dirname {output}) && \
            python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
                -name {params.output_prefix} -restrict {input.restrict} -scale {input.scale}
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
            )
        output:
            expand(
                f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react",
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            )
        conda:
            "../envs/structurefold.yaml"
        params:
            script=f"{workflow.basedir}/workflow/scripts/react_intersect_transcripts.py",
        shell:
            """
            python {params.script} -in {input} -out {output}
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
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react"
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity.csv"
        wildcard_constraints:
            id="[^/]+"
        conda:
            "../envs/biopython.yaml"
        params:
            transcriptome=TRANSCRIPTOME_UPPER,
        shell:
            """
            python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
            """


    rule heat_correct_reactivity:
        """
        Scale coverage- and transcript-set-restricted reactivity.react
        files (heat_correction_shared_react) from a lower- and a
        higher-temperature condition onto a common overall-signal scale
        with react_heat_correct.py, supporting multiple replicates per
        temperature.
        """
        input:
            lower=expand(f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react", id=HEAT_CORRECTION_LOWER_IDS),
            higher=expand(f"{config['output_dir']}/{{id}}/heat_correction/reactivity.react", id=HEAT_CORRECTION_HIGHER_IDS),
        output:
            lower=expand(f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react", id=HEAT_CORRECTION_LOWER_IDS),
            higher=expand(f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react", id=HEAT_CORRECTION_HIGHER_IDS),
        log:
            f"{config['output_dir']}/qc/heat_correction_scale_factors.log"
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
                -suffix {params.suffix} > {log} 2>&1
            """


    rule convert_heat_corrected_react_to_csv:
        """
        Convert heat-corrected reactivities to CSV form (transcript,
        position, base, reactivity -- one row per position), same as
        convert_react_to_csv.
        """
        input:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.react"
        output:
            f"{config['output_dir']}/{{id}}/heat_correction/reactivity_{HEAT_CORRECTION_SUFFIX}.csv"
        wildcard_constraints:
            id="[^/]+"
        conda:
            "../envs/biopython.yaml"
        params:
            transcriptome=TRANSCRIPTOME_UPPER,
        shell:
            """
            python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
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
                )
            output:
                expand(
                    f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react",
                    id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                )
            conda:
                "../envs/structurefold.yaml"
            params:
                script=f"{workflow.basedir}/workflow/scripts/react_intersect_transcripts.py",
            shell:
                """
                mkdir -p $(dirname {output[0]}) && \
                python {params.script} -in {input} -out {output}
                """


        rule heat_correct_reactivity_all_transcripts:
            """
            The no-coverage-filter counterpart to heat_correct_reactivity --
            same react_heat_correct.py rescaling, run on every transcript
            present in both conditions regardless of coverage
            (heat_correction_all_transcripts_shared_react).
            """
            input:
                lower=expand(f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react", id=HEAT_CORRECTION_LOWER_IDS),
                higher=expand(f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity.react", id=HEAT_CORRECTION_HIGHER_IDS),
            output:
                lower=expand(f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react", id=HEAT_CORRECTION_LOWER_IDS),
                higher=expand(f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react", id=HEAT_CORRECTION_HIGHER_IDS),
            log:
                f"{config['output_dir']}/qc/heat_correction_all_transcripts_scale_factors.log"
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
                    -suffix {params.suffix} > {log} 2>&1
                """


        rule convert_heat_correction_all_transcripts_to_csv:
            """
            Convert the all-transcripts (no coverage filter) heat-corrected
            reactivities to CSV form, same as convert_heat_corrected_react_to_csv.
            """
            input:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.react"
            output:
                f"{config['output_dir']}/{{id}}/heat_correction_all_transcripts/reactivity_{HEAT_CORRECTION_SUFFIX}.csv"
            wildcard_constraints:
                id="[^/]+"
            conda:
                "../envs/biopython.yaml"
            params:
                transcriptome=TRANSCRIPTOME_UPPER,
            shell:
                """
                python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
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
        minus=f"{TMP}/output/se/{{id}}/combined_minus.rtsc"
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react"
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
        ln_flag=lambda wc: "-ln_off" if wc.nlog == "noln" else "",
        nrm_flag=lambda wc: "-nrm_off" if wc.norm == "nonrm" else "",
    shell:
        """
        python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
            -name {params.output_prefix} \
            {params.trim_flag} {params.ln_flag} {params.nrm_flag}
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
        temperature = id[len("pooled_"):]
        rows = samples.loc[(samples["condition"] == condition) & (samples["temperature"] == temperature)]
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
        lambda wc: get_sams_by_condition_and_id(wc.condition, wc.id)
    output:
        f"{config['output_dir']}/{{id}}/depth_{{condition}}.txt"
    wildcard_constraints:
        id="[^/]+",
        condition="plus|minus",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        samtools depth -a {input} > {output}
        """


rule calculate_stop_coverage:
    input:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc"
    output:
        coverage=f"{config['output_dir']}/{{id}}/coverage.csv",
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/structurefold.yaml"
    params:
        coverage_name=f"{config['output_dir']}/{{id}}/coverage",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_coverage.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        mkdir -p $(dirname {output.coverage}) && \
        python {params.workdir}/{params.script} -f {input} -name {params.coverage_name} {params.transcriptome}
        """


rule rtsc_total_stops:
    """
    Sum RT-stop counts per transcript from the untreated (-DMS) channel to
    use as a raw read-count proxy for differential expression (DESeq2-style
    tools expect raw, un-normalized counts and do their own library-size
    normalization).
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc"
    output:
        f"{config['output_dir']}/{{id}}/counts_minus.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/upset.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python3 workflow/scripts/rtsc_total_stops.py --rtsc {input} --output {output}
        """


rule calculate_nucleotide_coverage:
    input:
        rtsc=f"{TMP}/output/se/{{id}}/combined_{{condition}}.rtsc"
    output:
        f"{config['output_dir']}/{{id}}/nucleotide_coverage_{{condition}}.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/structurefold.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/rtsc_to_nucleotide_coverage.py \
            --rtsc {input.rtsc} \
            --fasta {params.transcriptome} \
            --output {output}
        """


rule calculate_specificity:
    """
    Calculate the A/C specificity for the RT-stops
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_{{treatment}}.rtsc"
    output:
        f"{config['output_dir']}/{{id}}/specificity_{{treatment}}.csv"
    wildcard_constraints:
        id="[^/]+"
    params:
        output=f"{config['output_dir']}/{{id}}/specificity_{{treatment}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    conda:
        "../envs/structurefold.yaml"
    shell:
        """
        python {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome}
        """

rule convert_react_to_csv:
    """
    Convert the reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/reactivity.react"
    output:
        f"{config['output_dir']}/{{id}}/reactivity.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
        """


rule convert_react_plus_only_to_csv:
    """
    Convert the +DMS-only (no -DMS subtraction) reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.react"
    output:
        f"{config['output_dir']}/{{id}}/reactivity_plus_only.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
        """


rule convert_raw_react_plus_only_to_csv:
    """
    Convert raw (un-2-8%-scaled) +DMS-only reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.react"
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity_plus_only.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
        """


rule convert_raw_react_to_csv:
    """
    Convert raw (un-2-8%-scaled) reactivities to CSV form
    """
    input:
        f"{config['output_dir']}/{{id}}/raw_reactivity.react"
    output:
        f"{config['output_dir']}/{{id}}/raw_reactivity.csv"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/biopython.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
        """


rule convert_react_comparison_to_csv:
    input:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.react"
    output:
        f"{config['output_dir']}/{{id}}/norm_comparison/{{trim}}_{{nlog}}_{{norm}}/reactivity.csv"
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
        python3 workflow/scripts/react_to_csv.py --react {input} --output {output} --fasta {params.transcriptome}
        """

rule calculate_specificity_by_sample:
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc"
    output:
        f"{config['output_dir']}/specificity_{{sample}}.csv"
    params:
        output=f"{config['output_dir']}/specificity_{{sample}}",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/rtsc_specificity.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    conda:
        "../envs/structurefold.yaml"
    shell:
        """
        python {params.workdir}/{params.script} -rtsc {input} -name {params.output} -index {params.transcriptome}
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
        )
    output:
        f"{config['output_dir']}/qc/covered_transcripts.txt"
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python3 workflow/scripts/filter_covered_transcripts.py \
            --coverages {input} \
            --threshold 1 \
            --output {output}
        """


rule rtsc_stop_correlation:
    """
    Run StructureFold2 rtsc_correlation.py on the combined-plus RTSC files
    for all replicates/IDs, filtered to A/C positions and covered transcripts.
    """
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
            id=IDS,
        )
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv"
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
            -name {output}
        """


rule plot_rt_stop_correlation:
    """
    Pairwise hexbin scatter of log10(RT-stop + 1) counts at A and C positions
    across all replicates, with Pearson r annotated on each panel.
    """
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv"
    output:
        f"{config['output_dir']}/qc/rt_stop_replicate_correlation.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output}
        """


rule rtsc_stop_correlation_minus:
    """
    Same as rtsc_stop_correlation, but for the combined-minus (untreated,
    -DMS/background) RTSC files -- lets replicate agreement in the
    background channel be checked the same way as the treated channel.
    """
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
            id=IDS,
        )
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv"
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
            -name {output}
        """


rule plot_rt_stop_correlation_minus:
    """
    Same as plot_rt_stop_correlation, but for the combined-minus (untreated,
    -DMS/background) RT-stop counts.
    """
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv"
    output:
        f"{config['output_dir']}/qc/rt_stop_replicate_correlation_minus.png"
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output}
        """


rule rtsc_end_coverage_sweep:
    """
    Run rtsc_end_coverage.py in TP mode on the combined minus RTSC for trim
    values 0..200, recording the mean TP coverage at each trim amount.
    Only non-NA transcripts (those long enough for the given trim + tp_l) are
    included in each mean. On-demand only -- not in the default target list.
    """
    input:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc"
    output:
        f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.tsv"
    wildcard_constraints:
        id="[^/]+"
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
        printf 'trim\tmean_tp_coverage\n' > {output}
        for trim in $(seq 0 200); do
            tmpout={params.tmpdir}/trim_$trim.csv
            python {params.workdir}/{params.script} {input} TP \
                -length {params.length} -tp_l {params.tp_l} -trim $trim \
                -name $tmpout
            mean=$(awk -F',' \
                'NR>1 && $2 != "NA" {{sum+=$2; count++}} END {{if(count>0) printf "%.6f", sum/count; else print "NA"}}' \
                $tmpout)
            printf '%d\t%s\n' $trim $mean >> {output}
        done
        rm -rf {params.tmpdir}
        """


rule plot_end_coverage_sweep:
    """
    Plot mean 3' TP coverage vs. trim amount from the sweep TSV produced by
    rtsc_end_coverage_sweep. On-demand only.
    """
    input:
        f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.tsv"
    output:
        f"{config['output_dir']}/{{id}}/end_coverage_trim_sweep.png"
    wildcard_constraints:
        id="[^/]+"
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_end_coverage_sweep.py \
            --input {input} \
            --output {output} \
            --id {wildcards.id}
        """


rule plot_specificity:
    input:
        expand("{output_dir}/specificity_{sample}.csv", output_dir=config["output_dir"], sample=SAMPLES)
    output:
        f"{config['output_dir']}/qc/specificity_plot.png"
    params:
        samples=SAMPLES,
        conditions=samples["condition"].tolist(),
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_specificity.py \
            --inputs {input} \
            --samples {params.samples} \
            --conditions {params.conditions} \
            --output {output}
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
        _specificity_comparison_inputs
    output:
        plot=f"{config['output_dir']}/qc/specificity_comparison.png",
        manifest=f"{config['output_dir']}/qc/specificity_comparison_manifest.tsv",
    run:
        import pandas as pd
        rows = []
        for run in config["specificity_comparison_runs"]:
            sheet = pd.read_csv(run["samplesheet"], sep="\t")
            for _, row in sheet.iterrows():
                rows.append({
                    "file":      f"{run['output_dir']}/specificity_{row['sample']}.csv",
                    "group":     run["name"],
                    "condition": row["condition"],
                    "replicate": row.get("ID", row["sample"]),
                })
        pd.DataFrame(rows).to_csv(output.manifest, sep="\t", index=False)
        shell(
            "Rscript workflow/scripts/plot_specificity_combined.R"
            " {output.manifest} {output.plot}"
        )
