############################################################################
# RT-stop counting and combining samples into replicates (or pools)
############################################################################

# pool_replicates decides which combined .rtsc files feed the reactivity (samples.smk, helpers.smk):
#   none                          each replicate ID's own +DMS and -DMS (default)
#   minus / plus                  that channel pooled over every replicate; the other per replicate
#   both                          both channels pooled; one reactivity ID, "pooled"
#   both_by_temperature           both pooled within each temperature; IDs "pooled_<temperature>"
#   minus_all_plus_by_temperature -DMS pooled over everything, +DMS within each temperature
_ALLOWED_POOL_REPLICATES = {"none", "minus", "plus", "both", "both_by_temperature", "minus_all_plus_by_temperature"}
if POOL_REPLICATES not in _ALLOWED_POOL_REPLICATES:
    raise ValueError(
        f"config['pool_replicates'] must be one of {sorted(_ALLOWED_POOL_REPLICATES)}, got {POOL_REPLICATES!r}")


rule sam_to_rtsc:
    """Count RT stops (read 5' ends) per transcript position (methods.md §6)."""
    input:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
        tx_upper=get_rtsc_fasta_for_sample,
    output:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    log:
        "logs/sam_to_rtsc/{sample}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    shell:
        "python3 workflow/scripts/StructureFold3/sam_to_rtsc.py -single {input.sam} {input.tx_upper} > {log} 2>&1"


rule combine_rtsc:
    """
    Sum the RT stops of the samples in one replicate ID, or in a pool ("pooled", or
    "pooled_<temperature>"), for one condition (helpers.smk: samples_to_combine).
    """
    input:
        lambda wc: [sample_rtsc(s) for s in samples_to_combine(wc.condition, wc.id)],
    output:
        f"{TMP}/output/se/{{id}}/combined_{{condition}}.rtsc",
    log:
        "logs/combine_rtsc/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
        condition="plus|minus",
    conda:
        "../envs/rtsc_tools.yaml"
    shell:
        "python3 workflow/scripts/StructureFold3/rtsc_combine.py -name {output} {input} > {log} 2>&1"
