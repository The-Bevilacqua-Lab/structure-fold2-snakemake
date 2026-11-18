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


if DEPTH_MATCHING_SHARDS == 1:

    rule depth_matched_reactivity:
        """Thin every ID's +DMS to the same RT-stop total per transcript; mean reactivity over the draws."""
        input:
            plus=[get_plus_rtsc_for_id(id) for id in IDS],
            minus=[_depth_match_minus(id) for id in IDS],
            fasta=TRANSCRIPTOME_UPPER,
            # genotype mode: matched transcripts of the two genotypes share a key (genotypes.smk)
            key_map=[f"{GENOTYPE_DIR}/depth_matching_key_map.tsv"] if GENOTYPE_MODE else [],
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
            normalization=f"--normalization {DEPTH_MATCHING_NORMALIZATION} --min-gu {DEPTH_MATCHING_MIN_GU}",
            key_map=f"--key-map {GENOTYPE_DIR}/depth_matching_key_map.tsv" if GENOTYPE_MODE else "",
        shell:
            """
            python3 workflow/scripts/depth_matching/depth_matched_reactivity.py \
                --labels {params.ids} --plus {input.plus} --minus {input.minus} --fasta {input.fasta} \
                --outdirs {params.outdirs} --summary {output.summary} \
                --draws {params.draws} {params.match_minus} {params.normalization} {params.key_map} --threads {threads} >{log} 2>&1
            """

else:
    # depth_matching_shards: N one-thread jobs, each the depth matching of one contiguous block of
    # the sorted transcript list (seeds per transcript index, so the merge equals one unsharded run).
    DEPTH_MATCH_SHARD_DIR = f"{TMP}/depth_matched_shards"
    DEPTH_MATCH_SHARD_IDS = [f"{k:03d}" for k in range(DEPTH_MATCHING_SHARDS)]

    rule depth_matched_reactivity_shard:
        """Depth matching (as depth_matched_reactivity) for one shard of the transcripts."""
        input:
            plus=[get_plus_rtsc_for_id(id) for id in IDS],
            minus=[_depth_match_minus(id) for id in IDS],
            fasta=TRANSCRIPTOME_UPPER,
            key_map=[f"{GENOTYPE_DIR}/depth_matching_key_map.tsv"] if GENOTYPE_MODE else [],
        output:
            react=temp(expand(f"{DEPTH_MATCH_SHARD_DIR}/{{{{shard}}}}/{{id}}/reactivity.csv", id=IDS)),
            rtsc=temp(expand(f"{DEPTH_MATCH_SHARD_DIR}/{{{{shard}}}}/{{id}}/matched_plus_draw1.rtsc", id=IDS)),
            summary=temp(f"{DEPTH_MATCH_SHARD_DIR}/{{shard}}/summary.tsv"),
        log:
            "logs/depth_matched_reactivity_shard/{shard}.log",
        wildcard_constraints:
            shard="[0-9]+",
        threads: 1
        resources:
            mem_mb=24000,   # every shard loads all the libraries (~6-13 GB for 12 rice libraries)
            runtime=240,
        conda:
            "../envs/metagene.yaml"
        params:
            ids=IDS,
            outdirs=lambda wc: [f"{DEPTH_MATCH_SHARD_DIR}/{wc.shard}/{id}" for id in IDS],
            draws=DEPTH_MATCHING_DRAWS,
            n_shards=DEPTH_MATCHING_SHARDS,
            match_minus="--match-minus" if DEPTH_MATCHING_MATCH_MINUS else "",
            normalization=f"--normalization {DEPTH_MATCHING_NORMALIZATION} --min-gu {DEPTH_MATCHING_MIN_GU}",
            key_map=f"--key-map {GENOTYPE_DIR}/depth_matching_key_map.tsv" if GENOTYPE_MODE else "",
        shell:
            """
            python3 workflow/scripts/depth_matching/depth_matched_reactivity.py \
                --labels {params.ids} --plus {input.plus} --minus {input.minus} --fasta {input.fasta} \
                --outdirs {params.outdirs} --summary {output.summary} \
                --draws {params.draws} {params.match_minus} {params.normalization} {params.key_map} --threads 1 \
                --shard $((10#{wildcards.shard})) --n-shards {params.n_shards} >{log} 2>&1
            """

    rule depth_matched_reactivity:
        """Concatenate the depth-matching shards, in shard order."""
        input:
            react=expand(f"{DEPTH_MATCH_SHARD_DIR}/{{shard}}/{{id}}/reactivity.csv", shard=DEPTH_MATCH_SHARD_IDS, id=IDS),
            rtsc=expand(f"{DEPTH_MATCH_SHARD_DIR}/{{shard}}/{{id}}/matched_plus_draw1.rtsc", shard=DEPTH_MATCH_SHARD_IDS, id=IDS),
            summary=expand(f"{DEPTH_MATCH_SHARD_DIR}/{{shard}}/summary.tsv", shard=DEPTH_MATCH_SHARD_IDS),
        output:
            react=expand(f"{DEPTH_MATCH_DIR}/{{id}}/reactivity.csv", id=IDS),
            rtsc=expand(f"{DEPTH_MATCH_DIR}/{{id}}/matched_plus_draw1.rtsc", id=IDS),
            summary=f"{DEPTH_MATCH_DIR}/summary.tsv",
        log:
            "logs/depth_matched_reactivity.log",
        conda:
            "../envs/metagene.yaml"
        params:
            shard_dirs=[f"{DEPTH_MATCH_SHARD_DIR}/{k}" for k in DEPTH_MATCH_SHARD_IDS],
            ids=IDS,
            outdirs=[f"{DEPTH_MATCH_DIR}/{id}" for id in IDS],
        shell:
            """
            python3 workflow/scripts/depth_matching/merge_depth_matched_shards.py \
                --shard-dirs {params.shard_dirs} --labels {params.ids} --outdirs {params.outdirs} \
                --summary {output.summary} >{log} 2>&1
            """


rule depth_matched_to_rnaframework_xml:
    """Depth-matched reactivity.csv -> rnaframework_xml/<transcript>.xml (rf-norm format; config rnaframework_xml)."""
    input:
        csv=f"{DEPTH_MATCH_DIR}/{{id}}/reactivity.csv",
        transcriptome=TRANSCRIPTOME_UPPER,
    output:
        directory(f"{DEPTH_MATCH_DIR}/{{id}}/rnaframework_xml"),
    log:
        "logs/depth_matched_to_rnaframework_xml/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 workflow/scripts/reactivity/react_to_rnaframework_xml.py --csv {input.csv} \
            --fasta {input.transcriptome} --outdir {output} >{log} 2>&1
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
        "../envs/plotnine.yaml"
    shell:
        """
        python3 workflow/scripts/qc_plots/plot_rt_stop_correlation.py --input {input} --output {output} >{log} 2>&1
        """


# p4p6 structure plot and ROC on the depth-matched reactivities (targets requested in targets.smk
# when the positive control is p4p6); the ROC leaves out the last depth_matching_roc_trim3 nt.
if POSITIVE_CONTROL_IS_P4P6:

    rule plot_depth_matched_p4p6_all_samples:
        """Depth-matched p4p6 reactivity on the P4-P6 structure, every ID side by side."""
        input:
            react=expand(f"{DEPTH_MATCH_DIR}/{{id}}/reactivity.csv", id=IDS),
            xrna=P4P6_XRNA,
        output:
            f"{DEPTH_MATCH_QC}/p4p6_all_samples.png",
        log:
            "logs/plot_depth_matched_p4p6_all_samples.log",
        conda:
            "../envs/plotnine.yaml"
        params:
            labels=IDS,
        shell:
            """
            python3 workflow/scripts/qc_plots/plot_p4p6_structure.py \
                --inputs {input.react} --labels {params.labels} --xrna {input.xrna} \
                --output {output} >{log} 2>&1
            """

    rule plot_depth_matched_p4p6_roc:
        """ROC of depth-matched p4p6 reactivity: unpaired A/C (positive) vs WCF-paired A/C (negative)."""
        input:
            react=expand(f"{DEPTH_MATCH_DIR}/{{id}}/reactivity.csv", id=IDS),
            xrna=P4P6_XRNA,
        output:
            plot=f"{DEPTH_MATCH_QC}/p4p6_roc.png",
            auc=f"{DEPTH_MATCH_QC}/p4p6_auc.csv",
        log:
            "logs/plot_depth_matched_p4p6_roc.log",
        conda:
            "../envs/plotnine.yaml"
        params:
            labels=IDS,
            trim3=DEPTH_MATCHING_ROC_TRIM3,
        shell:
            """
            python3 workflow/scripts/qc_plots/plot_p4p6_roc.py \
                --inputs {input.react} --labels {params.labels} --xrna {input.xrna} --trim3 {params.trim3} \
                --title "p4p6 (depth-matched): reactivity vs known structure" \
                --plot {output.plot} --auc-csv {output.auc} >{log} 2>&1
            """
