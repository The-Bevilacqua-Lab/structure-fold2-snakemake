############################################################################
# Read composition by target type, and the self-contained HTML QC report
############################################################################

REPORT_HTML = f"{config['output_dir']}/qc/report.html"
COMPOSITION_DIR = f"{config['output_dir']}/qc/read_composition"
COMPOSITION_SUMMARY = f"{config['output_dir']}/qc/read_composition_summary.tsv"
COMPOSITION_PIE = f"{config['output_dir']}/qc/read_composition_pie.png"
COMPOSITION_PIE_MAPPED = f"{config['output_dir']}/qc/read_composition_pie_mapped.png"
COMPOSITION_PIE_TSV = f"{config['output_dir']}/qc/read_composition_pie.tsv"
# lncRNA is its own type only when an lncRNA set is configured (lncrna_fasta / lncrna_gff / lncrna_gtf)
COMPOSITION_TYPES = (["mRNA", "rRNA", "ncRNA"] + (["lncRNA"] if "lncrna" in EXTRA_RNA else [])
                     + ["positive_control", "other"])
# Pie slices: a read goes to the first of these types its forward-strand alignments hit, else
# "other" (other ncRNA, positive control, unannotated, antisense-only), or "unmapped"
if not config.get("read_composition_priority"):     # default list: keep the types this run has
    READ_COMPOSITION_PRIORITY = [t for t in READ_COMPOSITION_PRIORITY if t in COMPOSITION_TYPES]
_bad_priority = [t for t in READ_COMPOSITION_PRIORITY if t not in COMPOSITION_TYPES or t == "other"]
if _bad_priority:
    raise ValueError(f"read_composition_priority: unknown type(s) {_bad_priority}; "
                     f"choose from {[t for t in COMPOSITION_TYPES if t != 'other']}")
PIE_SLICES = list(READ_COMPOSITION_PRIORITY) + ["other", "unmapped"]


rule transcript_categories:
    """
    Transcript ID -> target type, from the FASTAs the transcriptome was built from:
    base transcriptome = mRNA (records whose header carries an Ensembl-style transcript_biotype
    other than protein_coding: rRNA -> rRNA, else ncRNA), the rRNA set = rRNA, the ncRNA set =
    ncRNA (rRNA where its gff/gtf marks the transcript as rRNA), the lncRNA set = lncRNA,
    positive_control_fasta = positive_control. Later sources override earlier ones.
    """
    input:
        base=_source_transcriptome,
        extras=_extra_transcriptome_fastas,
    output:
        f"{TMP}/resources/transcript_categories.tsv",
    params:
        # gffread headers carry no biotype: label the base records from the annotation(s) they came from
        base_ann=" ".join(
            [f"--base-annotation {g}{GENOTYPE_SEP} {GENOTYPES[g]['gff']}" for g in GENOTYPE_NAMES] if GENOTYPE_MODE
            else ([f"--base-annotation - {config['gff']}"] if TRANSCRIPT_SOURCE == "gffread" else [])),
        sets=" ".join(
            [f"--set ncRNA {EXTRA_RNA['ncrna']['fasta']} {EXTRA_RNA['ncrna'].get('annotation', '')}".rstrip()
             if "ncrna" in EXTRA_RNA else ""]
            + [f"--set lncRNA {EXTRA_RNA['lncrna']['fasta']}" if "lncrna" in EXTRA_RNA else ""]
            + [f"--set rRNA {EXTRA_RNA['rrna']['fasta']}" if "rrna" in EXTRA_RNA else ""]
            + [f"--set positive_control {config['positive_control_fasta']}"
               if config.get("positive_control_fasta") else ""]),
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/transcript_categories.log",
    shell:
        "python3 workflow/scripts/qc_plots/transcript_categories.py --base {input.base} {params.base_ann} {params.sets} "
        "--out {output} > {log} 2>&1"


rule read_composition:
    """
    Reads of one sample by target type (mRNA / rRNA / ncRNA / lncRNA / positive control), from
    every read in the merged, unfiltered SAM (unmapped reads included). Only forward-strand
    alignments give a type, as in the reactivity (sam_filter drops reverse-strand ones):
      unique_<type>  reads with exactly one alignment
      typed_<type>   reads whose forward alignments all hit that one type (reads shared between
                     isoforms or rRNA copies are kept); type_ambiguous = reads hitting >1 type;
                     antisense_only = mapped reads with no forward alignment
      pie_<slice>    every read in one slice: the first read_composition_priority type it hits,
                     else other; unmapped reads -> unmapped
    bowtie2 (-a) writes all alignments of a read together, so reads are grouped in one pass.
    """
    input:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam",
        categories=f"{TMP}/resources/transcript_categories.tsv",
    output:
        f"{COMPOSITION_DIR}/{{sample}}.tsv",
    log:
        "logs/read_composition/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    params:
        awk="workflow/scripts/qc_plots/read_composition.awk",
        types=" ".join(COMPOSITION_TYPES),
        pie=" ".join(READ_COMPOSITION_PRIORITY),
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view {input.sam} 2> {log} \
            | awk -F'\\t' -v OFS='\\t' -v s="{wildcards.sample}" -v types="{params.types}" -v pie="{params.pie}" \
                  -f {params.awk} {input.categories} - \
            > {output} 2>> {log}
        """


rule aggregate_read_composition:
    """One row per sample: total / unmapped / mapped / unique / multimapped reads, reads by target
    type, and the pie slices (pie_<slice> counts, pie_<slice>_pct of all reads)."""
    input:
        expand(f"{COMPOSITION_DIR}/{{sample}}.tsv", sample=SAMPLES),
    output:
        COMPOSITION_SUMMARY,
    log:
        "logs/aggregate_read_composition.log",
    run:
        long = pd.concat([pd.read_csv(f, sep="\t") for f in input], ignore_index=True)
        wide = long.pivot(index="sample", columns="metric", values="reads").reindex(SAMPLES)
        out = pd.DataFrame(index=wide.index)
        out["total_reads"] = wide["total_reads"]
        out["unmapped_reads"] = wide["unmapped_reads"]
        out["mapped_pct"] = 100 * wide["mapped_reads"] / wide["total_reads"]
        out["mapped_reads"] = wide["mapped_reads"]
        out["unique_reads"] = wide["unique_reads"]
        out["multimapped_reads"] = wide["multimapped_reads"]
        out["unique_pct_of_mapped"] = 100 * wide["unique_reads"] / wide["mapped_reads"]
        for t in COMPOSITION_TYPES:
            out[f"unique_{t}_pct"] = 100 * wide[f"unique_{t}"] / wide["unique_reads"]
        for t in COMPOSITION_TYPES:
            out[f"typed_{t}_pct"] = 100 * wide[f"typed_{t}"] / wide["mapped_reads"]
        out["type_ambiguous_pct"] = 100 * wide["type_ambiguous"] / wide["mapped_reads"]
        out["antisense_only_pct"] = 100 * wide["antisense_only"] / wide["mapped_reads"]
        for s in PIE_SLICES:
            out[f"pie_{s}_pct"] = 100 * wide[f"pie_{s}"] / wide["total_reads"]
        for t in COMPOSITION_TYPES:
            out[f"unique_{t}"] = wide[f"unique_{t}"]
            out[f"typed_{t}"] = wide[f"typed_{t}"]
        out["type_ambiguous"] = wide["type_ambiguous"]
        out["antisense_only"] = wide["antisense_only"]
        for s in PIE_SLICES:
            out[f"pie_{s}"] = wide[f"pie_{s}"]
        out.round(3).to_csv(output[0], sep="\t")


rule plot_read_composition_pie:
    """One pie per sample (+ all samples pooled): mRNA / rRNA / lncRNA / other / unmapped; and the same
    for mapped reads only (no unmapped slice)."""
    input:
        COMPOSITION_SUMMARY,
    output:
        png=COMPOSITION_PIE,
        png_mapped=COMPOSITION_PIE_MAPPED,
        tsv=COMPOSITION_PIE_TSV,
    params:
        slices=PIE_SLICES,
    log:
        "logs/plot_read_composition_pie.log",
    conda:
        "../envs/plotnine.yaml"
    script:
        "../scripts/qc_plots/plot_read_composition_pie.py"


def report_inputs(wildcards):
    """Everything the HTML report embeds; mirrors the conditions in get_all_targets (targets.smk)."""
    out = config["output_dir"]
    d = {
        "alignment": f"{out}/qc/alignment_stats_summary.tsv",
        "composition": COMPOSITION_SUMMARY,
        "composition_pie": COMPOSITION_PIE,
        "composition_pie_mapped": COMPOSITION_PIE_MAPPED,
        "specificity_plot": f"{out}/qc/specificity_plot.png",
        "specificity_tables": expand(f"{out}/specificity_{{sample}}.csv", sample=SAMPLES),
        "metagene_depth": f"{out}/qc/metagene/metagene_depth.png",
        "metagene_rtsc": f"{out}/qc/metagene/metagene_rtsc.png",
        "coverage": expand(f"{out}/{{id}}/coverage.csv", id=IDS_WITH_PLUS),
        "logo": "resources/sf2_snakemake_logo.png",
    }
    if len(CORRELATION_PLUS_IDS) > 1:
        d["rt_corr"] = f"{out}/qc/rt_stop_replicate_correlation.png"
        d["upset_replicates"] = f"{out}/qc/covered_transcripts_upset_replicates.png"
    if len(CORRELATION_MINUS_IDS) > 1:
        d["rt_corr_minus"] = f"{out}/qc/rt_stop_replicate_correlation_minus.png"
    if MERGED_IDS_DIFFER and len(IDS) > 1:
        d["upset_merged"] = f"{out}/qc/covered_transcripts_upset_merged.png"
    if config.get("positive_control_fasta"):
        d["positive_control"] = f"{out}/qc/positive_control_alignment_summary.tsv"
    if POSITIVE_CONTROL_IS_P4P6:
        # only the final reactivities: depth-matched when depth matching is on (built from the
        # 3' bias-corrected stops when that is on too), else the standard ones
        p4p6_qc = f"{out}/qc/depth_matched" if DEPTH_MATCHING else f"{out}/qc"
        d["p4p6_all"] = f"{p4p6_qc}/p4p6_all_samples.png"
        d["p4p6_roc"] = f"{p4p6_qc}/p4p6_roc.png"
        d["p4p6_auc"] = f"{p4p6_qc}/p4p6_auc.csv"
    if BIAS_CORRECTION:
        d["bias_models"] = expand(f"{out}/{{id}}/3prime_bias_correction/model.tsv", id=IDS)
        d["bias_coverage_plots"] = expand(f"{out}/{{id}}/3prime_bias_correction/coverage_residual_metagene.png", id=IDS)
        d["bias_rtsc_plots"] = expand(f"{out}/{{id}}/3prime_bias_correction/rtsc_metagene_before_after.png", id=IDS)
    if DEPTH_MATCHING:
        d["depth_matched_corr"] = f"{out}/qc/depth_matched/rt_stop_replicate_correlation.png"
        d["depth_matched_summary"] = f"{out}/depth_matched/summary.tsv"
    if GENOME_CHECK:
        d["genome_check_summary"] = f"{out}/genome_check/summary.tsv"
        d["genome_check_loci"] = f"{out}/genome_check/top_loci.tsv"
        d["genome_check_plot"] = f"{out}/genome_check/categories.png"
    if HEAT_CORRECTION:
        d["heat_shared"] = f"{out}/qc/heat_correction_shared_transcripts.txt"
    if GENOTYPE_MODE:
        # per-genotype replicate QC replaces the run-wide versions; no matched RT-stop correlation
        for key in ("rt_corr", "rt_corr_minus", "upset_replicates", "depth_matched_corr"):
            d.pop(key, None)
        plus_g = [g for g in GENOTYPE_NAMES if len(GENOTYPE_PLUS_IDS[g]) > 1]
        minus_g = [g for g in GENOTYPE_NAMES if len(GENOTYPE_MINUS_IDS[g]) > 1]
        if plus_g:
            d["rt_corr"] = [f"{out}/qc/rt_stop_replicate_correlation_plus_{g}.png" for g in plus_g]
            d["upset_replicates"] = [f"{out}/qc/covered_transcripts_upset_replicates_{g}.png" for g in plus_g]
        if minus_g:
            d["rt_corr_minus"] = [f"{out}/qc/rt_stop_replicate_correlation_minus_{g}.png" for g in minus_g]
        d["genotype_liftover"] = GENOTYPE_LIFTOVER
        d["genotype_summary"] = f"{out}/genotype_comparison/reactivity_comparison_summary.tsv"
        d["genotype_plot"] = f"{out}/genotype_comparison/reactivity_comparison.png"
    return d


def _report_settings():
    s = {
        "output_dir": config["output_dir"],
        "samplesheet": config["samplesheet"],
        "transcript_source": TRANSCRIPT_SOURCE,
        "trimmer": TRIMMER,
        "max_mismatch (sam_filter)": MAX_MISMATCH,
        "pool_replicates": POOL_REPLICATES,
        "3_prime_bias_correction": BIAS_CORRECTION,
        "depth_matching": DEPTH_MATCHING,
        "heat correction (temperature column)": HEAT_CORRECTION,
        "positive_control_name": config.get("positive_control_name", "—"),
        "ncRNA set": EXTRA_RNA.get("ncrna", {}).get("source", "—"),
        "lncRNA set": EXTRA_RNA.get("lncrna", {}).get("source", "—"),
        "rRNA set": EXTRA_RNA.get("rrna", {}).get("source", "—"),
        "read_composition_priority": " > ".join(READ_COMPOSITION_PRIORITY),
    }
    if TRANSCRIPT_SOURCE == "gffread":
        s["cds_only"] = CDS_ONLY
        s["keep_longest_isoform"] = KEEP_LONGEST_ISOFORM
    if TRIMMER == "fastp":
        s["fastp_adapter"] = FASTP_ADAPTER
        s["fastp_min_length"] = FASTP_MIN_LENGTH
    if BIAS_CORRECTION:
        s["3_prime_bias_reference"] = BIAS_REFERENCE
        s["3_prime_bias_training"] = BIAS_TRAINING
    if GENOTYPE_MODE:
        for g in GENOTYPE_NAMES:
            s[f"genotype {g}"] = f"{GENOTYPES[g]['genome']} + {GENOTYPES[g]['gff']}"
        s["transcript_matching"] = TRANSCRIPT_MATCHING or "identical IDs"
        s["pool_minus_across_genotypes"] = POOL_MINUS_ACROSS_GENOTYPES
    if GENOME_CHECK:
        for r, spec in GENOME_CHECK_REFS.items():
            s[f"genome_check ({r})"] = (f"{spec['index'] or spec['genome']} + "
                                       f"{spec['annotation'] or 'no annotation'}")
    s["reactivity_normalization"] = REACTIVITY_NORMALIZATION + (
        f" (>= {REACTIVITY_MIN_GU} G/U positions)" if REACTIVITY_NORMALIZATION == "gu" else "")
    if DEPTH_MATCHING:
        s["depth_matching_draws"] = DEPTH_MATCHING_DRAWS
        s["depth_matching_match_minus"] = DEPTH_MATCHING_MATCH_MINUS
        s["depth_matching_normalization"] = DEPTH_MATCHING_NORMALIZATION + (
            f" (>= {DEPTH_MATCHING_MIN_GU} G/U positions)" if DEPTH_MATCHING_NORMALIZATION == "gu" else "")
    return s


rule qc_report:
    """Self-contained HTML report (tables + embedded figures) of every QC output of the run."""
    input:
        unpack(report_inputs),
    output:
        REPORT_HTML,
    params:
        title=config.get("report_title", "structure-fold2-snakemake QC report"),
        settings=_report_settings(),
        samples=samples_runs.astype(str).to_dict(orient="list"),
        ids=IDS,
        multiqc_link="multiqc_fastp/multiqc_report.html" if TRIMMER == "fastp" else "",
        composition_types=COMPOSITION_TYPES,
        pie_slices=PIE_SLICES,
        genotype_names=GENOTYPE_NAMES,
    log:
        "logs/qc_report.log",
    conda:
        "../envs/plotnine.yaml"
    script:
        "../scripts/qc_plots/make_report.py"
