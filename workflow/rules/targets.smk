############################################################################
# Final outputs requested by rule all, assembled from the samplesheet-derived
# ID lists (samples.smk) and the optional features enabled in config.yaml.
############################################################################


def get_all_targets(wildcards):
    out = config["output_dir"]
    targets = (
        expand("{out}/{id}/reactivity.react", out=out, id=IDS)
        + expand("{out}/{id}/reactivity.csv", out=out, id=IDS)
        # coverage.csv/specificity_plus.csv always read the true per-ID +DMS
        # channel (no pooled fallback -- see IDS_WITH_PLUS in samples.smk), so an ID
        # whose entire +DMS was excluded via exclude_samples_plus is skipped
        # here rather than requesting an impossible target.
        + expand("{out}/{id}/coverage.csv", out=out, id=IDS_WITH_PLUS)
        + expand("{out}/{id}/specificity_plus.csv", out=out, id=IDS_WITH_PLUS)
        # counts_minus.csv/specificity_minus.csv: same, for the -DMS channel
        # (IDS_WITH_MINUS) -- this is what exclude_samples_minus actually
        # needed here (see IDS_WITH_MINUS in samples.smk).
        + expand("{out}/{id}/counts_minus.csv", out=out, id=IDS_WITH_MINUS)
        + expand("{out}/{id}/specificity_minus.csv", out=out, id=IDS_WITH_MINUS)
        + expand("{out}/specificity_{sample}.csv", out=out, sample=SAMPLES)
        + [
            f"{out}/qc/alignment_stats_summary.tsv",
            f"{out}/qc/specificity_plot.png",
            f"{out}/qc/metagene/metagene_depth.png",
            f"{out}/qc/metagene/metagene_rtsc.png",
        ]
        + (
            # Correlation is computed on CORRELATION_PLUS_IDS (true
            # biological replicates that still have a non-excluded +DMS
            # sample), not IDS, so gate on that count too -- otherwise a
            # pooled run (IDS == ["pooled"] or ["pooled_<t>", ...]) would
            # skip this QC even when there are multiple real replicates to
            # compare.
            [f"{out}/qc/rt_stop_replicate_correlation.png"]
            if len(CORRELATION_PLUS_IDS) > 1
            else []
        )
        + (
            # Same, for the -DMS channel -- gated separately since
            # exclude_samples_plus/exclude_samples_minus can leave the two
            # channels with different replicate counts.
            [f"{out}/qc/rt_stop_replicate_correlation_minus.png"]
            if len(CORRELATION_MINUS_IDS) > 1
            else []
        )
        + (
            # Always produced (per-replicate, unpooled -- see
            # upset_covered_transcripts_replicates in coverage.smk,
            # which depends only on each replicate's +DMS channel), same
            # gate as the +DMS correlation plot above: an UpSet plot of a
            # single set has nothing to overlap.
            [f"{out}/qc/covered_transcripts_upset_replicates.png"]
            if len(CORRELATION_PLUS_IDS) > 1
            else []
        )
        + (
            # Only when pool_replicates has actually merged replicates into
            # a different (and still >1-set) grouping -- see
            # MERGED_IDS_DIFFER in samples.smk and upset_covered_transcripts_merged
            # in coverage.smk. "both" pooling collapses to the single
            # ID "pooled", which can't be UpSet-plotted against anything.
            [f"{out}/qc/covered_transcripts_upset_merged.png"]
            if MERGED_IDS_DIFFER and len(IDS) > 1
            else []
        )
        + (
            [f"{out}/transcript_position_annotations.csv"]
            if config.get("annotation_gtf") and config.get("genome")
            else []
        )
    )

    if ABUNDANCE_MODES:
        # abundance_<mode>.csv is computed from the true per-ID -DMS channel
        # (see calculate_transcript_abundance in abundance.smk), same
        # IDS_WITH_MINUS gating as counts_minus.csv above.
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
        targets += [f"{out}/qc/depth_matched/rt_stop_replicate_correlation.png"]
        if POSITIVE_CONTROL_IS_P4P6:
            targets += [f"{out}/qc/depth_matched/p4p6_roc.png"]

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

    return targets
