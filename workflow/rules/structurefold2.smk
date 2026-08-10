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
#              free, with no rule changes needed.
#
# QC rules (coverage, specificity, nucleotide coverage, replicate
# correlation, ...) always report on the true per-replicate data except
# under "both", where there is no per-replicate data left to report on.
# ---------------------------------------------------------------------------
_ALLOWED_POOL_REPLICATES = {"none", "minus", "plus", "both"}
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
    Count RT-stops from the SAM file
    """
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam"
    output:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc"
    conda:
        "../envs/structurefold.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/sam_to_rtsc.py",
        transcriptome=TRANSCRIPTOME,
    shell:
        """
        python {params.workdir}/{params.script} -single {input} {params.transcriptome}
        """


rule combine_rtsc_plus:
    """
    Combine the RT-stop counts of the sample(s) making up one +DMS
    replicate (ID). wildcard_constraints excludes the literal ID "pooled"
    so this rule can never collide with combine_rtsc_plus_pooled below.
    """
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("plus", wc.id)
    output:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc"
    wildcard_constraints:
        id="(?!pooled$).+"
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
        id="(?!pooled$).+"
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
        output_prefix=f"{config['output_dir']}/{{id}}/raw_reactivity",
    shell:
        """
        mkdir -p $(dirname {output}) && \
        python {params.workdir}/{params.script} {input.minus} {input.plus} {params.transcriptome} \
            -name {params.output_prefix} -nrm_off -threshold 1000000
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
        transcriptome=TRANSCRIPTOME,
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
    # get_rtsc_counts_by_condition's pooling behavior.
    if id == "pooled":
        rows = samples.loc[samples["condition"] == condition]
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
        transcriptome=TRANSCRIPTOME,
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
