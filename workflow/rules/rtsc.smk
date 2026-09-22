############################################################################
# RT-stop counting and combining replicates
############################################################################

# config['pool_replicates'] pools RT-stop counts across replicates before the
# reactivity calculation: none (default) | minus | plus | both |
# both_by_temperature | minus_all_plus_by_temperature. "both" gives the single
# ID "pooled"; the by-temperature modes give "pooled_<temperature>" IDs (they
# need a 'temperature' samplesheet column). Per-ID QC rules run on whatever
# IDS ends up as; replicate correlation always uses the unpooled replicates.
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



rule sam_to_rtsc:
    """Count RT-stops per transcript position from a filtered SAM."""
    input:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
        tx_upper=TRANSCRIPTOME_UPPER,
    output:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    log:
        "logs/sam_to_rtsc/{sample}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/sam_to_rtsc.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 {params.workdir}/{params.script} -single {input.sam} {params.transcriptome} >{log} 2>&1
        """


rule combine_rtsc_plus:
    """Combine the per-sample RT-stop counts of one +DMS replicate."""
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("plus", wc.id),
    output:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus/{id}.log",
    wildcard_constraints:
        id="(?!pooled(_|$)).+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus:
    """Combine the per-sample RT-stop counts of one -DMS replicate."""
    input:
        lambda wc: get_rtsc_counts_by_condition_and_id("minus", wc.id),
    output:
        f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus/{id}.log",
    wildcard_constraints:
        id="(?!pooled(_|$)).+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/{{id}}_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_plus_pooled:
    """Pool +DMS RT-stop counts across all replicates."""
    input:
        get_rtsc_counts_by_condition("plus"),
    output:
        f"{TMP}/output/se/pooled/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus_pooled.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus_pooled:
    """Pool -DMS RT-stop counts across all replicates."""
    input:
        get_rtsc_counts_by_condition("minus"),
    output:
        f"{TMP}/output/se/pooled/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus_pooled.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_plus_pooled_by_temperature:
    """Pool +DMS RT-stop counts within one temperature group."""
    input:
        lambda wc: get_rtsc_counts_by_condition_and_temperature("plus", wc.temperature),
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_plus.rtsc",
    log:
        "logs/combine_rtsc_plus_pooled_by_temperature/{temperature}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_plus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_plus {input} \
                && mv {params.tmpdir}/combined_plus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """


rule combine_rtsc_minus_pooled_by_temperature:
    """Pool -DMS RT-stop counts within one temperature group."""
    input:
        lambda wc: get_rtsc_counts_by_condition_and_temperature("minus", wc.temperature),
    output:
        f"{TMP}/output/se/pooled_{{temperature}}/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus_pooled_by_temperature/{temperature}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{temperature}}_combined_minus",
    shell:
        """
        (
            mkdir -p {params.tmpdir} \
                && cp {input} {params.tmpdir}/ \
                && python3 {params.workdir}/{params.script} -name {params.tmpdir}/combined_minus {input} \
                && mv {params.tmpdir}/combined_minus.rtsc {output} \
                && rm -rf {params.tmpdir}
        ) >{log} 2>&1
        """
