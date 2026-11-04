############################################################################
# Depth-matched reactivity (config key depth_matching)
#
############################################################################

from types import SimpleNamespace

DEPTH_MATCH_DIR = f"{config['output_dir']}/depth_matched"
DEPTH_MATCH_QC = f"{config['output_dir']}/qc/depth_matched"

# Matching to the smallest depth needs at least two reactivity IDs; with one
# (a single replicate, or pool_replicates: both) there is nothing to match.
if len(IDS) < 2:
    raise ValueError(
        f"depth_matching needs at least 2 reactivity IDs, got {IDS} -- "
        "turn it off or use a pool_replicates mode that keeps several IDs."
    )


def _depth_match_minus(id):
    return get_minus_rtsc_for_reactivity(SimpleNamespace(id=id))


rule depth_matched_reactivity:
    """Thin every ID's +DMS to the same RT-stop total per transcript; mean reactivity over the draws."""
    input:
        plus=[get_plus_rtsc_for_id(id) for id in IDS],
        minus=[_depth_match_minus(id) for id in IDS],
        fasta=TRANSCRIPTOME_UPPER,
    output:
        react=expand(f"{DEPTH_MATCH_DIR}/{{id}}/reactivity.csv", id=IDS),
        rtsc=expand(f"{DEPTH_MATCH_DIR}/{{id}}/matched_plus_draw1.rtsc", id=IDS),
        summary=f"{DEPTH_MATCH_DIR}/summary.tsv",
    log:
        "logs/depth_matched_reactivity.log",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=360,
    conda:
        "../envs/metagene.yaml"
    params:
        ids=IDS,
        outdirs=[f"{DEPTH_MATCH_DIR}/{id}" for id in IDS],
        draws=DEPTH_MATCHING_DRAWS,
        match_minus="--match-minus" if DEPTH_MATCHING_MATCH_MINUS else "",
    shell:
        """
        python3 workflow/scripts/depth_matched_reactivity.py \
            --labels {params.ids} --plus {input.plus} --minus {input.minus} --fasta {input.fasta} \
            --outdirs {params.outdirs} --summary {output.summary} \
            --draws {params.draws} {params.match_minus} --threads {threads} >{log} 2>&1
        """


rule depth_matched_rtsc_correlation:
    """RT-stop correlation of the matched +DMS (first draw) across IDs."""
    input:
        expand(f"{DEPTH_MATCH_DIR}/{{id}}/matched_plus_draw1.rtsc", id=IDS),
    output:
        f"{DEPTH_MATCH_QC}/rt_stop_correlation.csv",
    log:
        "logs/depth_matched_rtsc_correlation.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_correlation.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        mkdir -p $(dirname {output})
        python3 {params.workdir}/{params.script} {input} -spec AC -fasta {params.transcriptome} \
            -name {output} >{log} 2>&1
        """


rule plot_depth_matched_rt_stop_correlation:
    """Pairwise replicate correlation plot of the matched +DMS RT-stops."""
    input:
        f"{DEPTH_MATCH_QC}/rt_stop_correlation.csv",
    output:
        f"{DEPTH_MATCH_QC}/rt_stop_replicate_correlation.png",
    log:
        "logs/plot_depth_matched_rt_stop_correlation.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py --input {input} --output {output} >{log} 2>&1
        """
