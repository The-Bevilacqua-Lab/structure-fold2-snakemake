############################################################################
# Final outputs requested by rule all, assembled from the samplesheet-derived
# ID lists (samples.smk) and the optional features enabled in config.yaml.
############################################################################


def get_all_targets(wildcards):
    out = config["output_dir"]
    targets = (
        expand("{out}/{id}/reactivity.{ext}", out=out, id=IDS, ext=["react", "csv"])
        # Per-ID QC reads the ID's own combined .rtsc, so an ID whose +DMS (or -DMS) samples were
        # all excluded gets none for that channel (IDS_WITH_PLUS / IDS_WITH_MINUS, samples.smk).
        + expand("{out}/{id}/{f}.csv", out=out, id=IDS_WITH_PLUS, f=["coverage", "specificity_plus"])
        + expand("{out}/{id}/{f}.csv", out=out, id=IDS_WITH_MINUS, f=["counts_minus", "specificity_minus"])
        + expand("{out}/specificity_{sample}.csv", out=out, sample=SAMPLES)
        + [
            f"{out}/qc/alignment_stats_summary.tsv",
            f"{out}/qc/specificity_plot.png",
            f"{out}/qc/metagene/metagene_depth.png",
            f"{out}/qc/metagene/metagene_rtsc.png",
            f"{out}/qc/read_composition_summary.tsv",
            f"{out}/qc/read_composition_pie.png",
            f"{out}/qc/read_composition_pie_mapped.png",
            f"{out}/qc/report.html",
        ]
    )
    # Replicate QC compares the real replicates (unpooled), so it needs at least two of them.
    if len(CORRELATION_PLUS_IDS) > 1:
        targets += [f"{out}/qc/rt_stop_replicate_correlation.png",
                    f"{out}/qc/covered_transcripts_upset_replicates.png"]
    if len(CORRELATION_MINUS_IDS) > 1:
        targets += [f"{out}/qc/rt_stop_replicate_correlation_minus.png"]
    # Coverage overlap of the pooled IDs, when pooling made more than one
    if MERGED_IDS_DIFFER and len(IDS) > 1:
        targets += [f"{out}/qc/covered_transcripts_upset_merged.png"]
    if config.get("annotation_gtf") and config.get("genome"):
        targets += [f"{out}/transcript_position_annotations.csv"]

    if ABUNDANCE_MODES:
        targets += expand(
            "{out}/{id}/abundance_{mode}.csv", out=out, id=IDS_WITH_MINUS, mode=ABUNDANCE_MODES
        )

    if config.get("positive_control_fasta"):
        targets += [f"{out}/qc/positive_control_alignment_summary.tsv"]

    if POSITIVE_CONTROL_IS_P4P6:
        targets += (
            expand("{out}/{id}/p4p6_react.png", out=out, id=IDS)
            + [
                f"{out}/qc/p4p6_all_samples.png",
                f"{out}/qc/p4p6_roc.png",
            ]
        )

    if config.get("specificity_comparison_runs"):
        targets += [f"{out}/qc/specificity_comparison.png"]

    if TRIMMER == "fastp":
        targets += [f"{out}/qc/multiqc_fastp/multiqc_report.html"]

    if BIAS_CORRECTION:
        targets += expand(
            "{out}/{id}/3prime_bias_correction/{f}",
            out=out,
            id=IDS,
            f=[
                "model.tsv",
                "coverage_residual_metagene.png",
                "rtsc_metagene_before_after.png",
            ],
        )

    if DEPTH_MATCHING:
        targets += expand("{out}/depth_matched/{id}/reactivity.csv", out=out, id=IDS)
        if not GENOTYPE_MODE:   # the two genotypes' matched RT-stops sit under different names
            targets += [f"{out}/qc/depth_matched/rt_stop_replicate_correlation.png"]
        if POSITIVE_CONTROL_IS_P4P6:
            targets += [f"{out}/qc/depth_matched/p4p6_roc.png", f"{out}/qc/depth_matched/p4p6_all_samples.png"]

    if HEAT_CORRECTION:
        targets += expand(
            "{out}/{id}/heat_correction/reactivity_{suffix}.{ext}",
            out=out,
            id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
            suffix=HEAT_CORRECTION_SUFFIX,
            ext=["react", "csv"],
        ) + expand(
            "{out}/{id}/heat_correction/reactivity.csv",
            out=out,
            id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
        )

        if HEAT_CORRECTION_COMPARE_ALL_TRANSCRIPTS:
            targets += expand(
                "{out}/{id}/heat_correction_all_transcripts/reactivity_{suffix}.{ext}",
                out=out,
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                suffix=HEAT_CORRECTION_SUFFIX,
                ext=["react", "csv"],
            )

    if RNAFRAMEWORK_XML:
        targets += expand("{out}/{id}/rnaframework_xml", out=out, id=IDS)
        if DEPTH_MATCHING:
            targets += expand("{out}/depth_matched/{id}/rnaframework_xml", out=out, id=IDS)
        if HEAT_CORRECTION:
            targets += expand(
                "{out}/{id}/heat_correction/rnaframework_xml_{suffix}",
                out=out,
                id=HEAT_CORRECTION_LOWER_IDS + HEAT_CORRECTION_HIGHER_IDS,
                suffix=HEAT_CORRECTION_SUFFIX,
            )

    if GENOME_CHECK:
        targets += [f"{out}/genome_check/{f}" for f in ("summary.tsv", "top_loci.tsv", "categories.png")]

    if GENOTYPE_MODE:
        targets = [t for t in targets if t not in GENOTYPE_REPLACED_TARGETS]
        targets += genotype_targets()

    return targets


# Run-wide replicate QC that genotype mode replaces with per-genotype versions (genotypes.smk).
GENOTYPE_REPLACED_TARGETS = {
    f"{config['output_dir']}/qc/{f}"
    for f in ("rt_stop_replicate_correlation.png", "rt_stop_replicate_correlation_minus.png",
              "covered_transcripts_upset_replicates.png")
}


def genotype_targets():
    """Genotype mode: per-genotype replicate QC, the matched transcripts, the reactivity comparison,
    and the per-position annotation (each genotype's own GFF)."""
    out = config["output_dir"]
    t = [GENOTYPE_LIFTOVER, f"{out}/genotype_comparison/reactivity_comparison.png",
         f"{out}/transcript_position_annotations.csv"]
    for g in GENOTYPE_NAMES:
        if len(GENOTYPE_PLUS_IDS[g]) > 1:
            t += [f"{out}/qc/rt_stop_replicate_correlation_plus_{g}.png",
                  f"{out}/qc/covered_transcripts_upset_replicates_{g}.png"]
        if len(GENOTYPE_MINUS_IDS[g]) > 1:
            t += [f"{out}/qc/rt_stop_replicate_correlation_minus_{g}.png"]
    return t
