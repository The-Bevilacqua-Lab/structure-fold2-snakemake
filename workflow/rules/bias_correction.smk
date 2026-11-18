############################################################################
# Optional: Correct for differences in 3' coverage-bias correction of
# the +DMS RT-stops because of poly-(A) purification
#
############################################################################

BIAS_QC = f"{config['output_dir']}/{{id}}/3prime_bias_correction"
BIAS_TRAINING_TSV = f"{config['output_dir']}/qc/3prime_bias_correction/training_transcripts.tsv"


# genotype mode: per-genotype lists, genotypes.smk
if not GENOTYPE_MODE:

    rule bias_protein_coding_transcripts:
        """Protein-coding transcripts (from the annotation) present in the transcriptome."""
        input:
            annotation=BIAS_ANNOTATION,
            fasta=TRANSCRIPTOME_UPPER,
        output:
            f"{BIAS_TMP}/protein_coding_transcripts.txt",
        log:
            "logs/bias_protein_coding_transcripts.log",
        conda:
            "../envs/metagene.yaml"
        shell:
            """
            python3 workflow/scripts/three_prime_end_bias/list_protein_coding_transcripts.py \
                --annotation {input.annotation} --fasta {input.fasta} --output {output} >{log} 2>&1
            """


rule bias_training_transcripts:
    """
    One training set shared by every ID, from each protein-coding transcript's
    lowest +DMS RT-stop coverage across all IDs: the top-N
    (3_prime_bias_training: top_n) or every transcript at or above
    3_prime_bias_min_coverage in all IDs (all_covered).
    """
    input:
        rtsc=[get_raw_plus_rtsc(id_) for id_ in IDS],
        fasta=TRANSCRIPTOME_UPPER,
        protein_coding=f"{BIAS_TMP}/protein_coding_transcripts.txt",
    output:
        tsv=BIAS_TRAINING_TSV,
        bed=f"{BIAS_TMP}/training_transcripts.bed",
    log:
        "logs/bias_training_transcripts.log",
    conda:
        "../envs/metagene.yaml"
    params:
        selection=(
            f"--mode top_n --n {BIAS_N_TRANSCRIPTS}"
            if BIAS_TRAINING == "top_n"
            else f"--mode all_covered --min-coverage {BIAS_MIN_COVERAGE}"
        ),
        labels=IDS,
        # genotype mode: score each transcript over its own genotype's IDs only (top_n per genotype)
        groups=(f"--groups {' '.join(ID_GENOTYPE[i] for i in IDS)} --group-sep {GENOTYPE_SEP}" if GENOTYPE_MODE else ""),
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {output.bed}) \
            && python3 workflow/scripts/three_prime_end_bias/select_bias_training_transcripts.py \
                --rtsc {input.rtsc} --labels {params.labels} --fasta {input.fasta} \
                --protein-coding {input.protein_coding} {params.selection} {params.groups} \
                --output {output.tsv} --bed {output.bed} >{log} 2>&1
        """


rule bias_training_depth:
    """
    Per-nucleotide read depth over the training transcripts, summed across
    the channel's samples (+DMS, or the -DMS reference).
    """
    input:
        bams=get_bias_bams,
        bed=f"{BIAS_TMP}/training_transcripts.bed",
    output:
        temp(f"{BIAS_TMP}/{{id}}/depth_{{condition}}.tsv.gz"),
    log:
        "logs/bias_training_depth/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
        condition="plus|minus",
    conda:
        "../envs/metagene.yaml"
    params:
        # genotype mode: only the ID's own genotype's training transcripts are in its BAMs
        prefix=lambda wc: f"{ID_GENOTYPE[wc.id]}{GENOTYPE_SEP}" if GENOTYPE_MODE else "",
    shell:
        """
        awk -v p="{params.prefix}" 'index($1, p) == 1' {input.bed} > {output}.bed
        samtools depth -a -b {output}.bed {input.bams} 2>{log} \
            | awk -F'\\t' -v OFS='\\t' '{{s = 0; for (i = 3; i <= NF; i++) s += $i; print $1, $2, s}}' \
            | gzip >{output}
        rm -f {output}.bed
        """


rule fit_3prime_bias_model:
    """Fit the global 3' coverage-bias model on the training transcripts."""
    input:
        plus=f"{BIAS_TMP}/{{id}}/depth_plus.tsv.gz",
        minus=f"{BIAS_TMP}/{{id}}/depth_minus.tsv.gz",
        training=BIAS_TRAINING_TSV,
    output:
        model=f"{BIAS_QC}/model.tsv",
        binned=f"{BIAS_QC}/coverage_residual_metagene.tsv",
        plot=report(
            f"{BIAS_QC}/coverage_residual_metagene.png",
            category="3' bias correction",
            subcategory="{id}",
            caption="../report_captions/bias_coverage_residual_metagene.rst",
        ),
    log:
        "logs/fit_3prime_bias_model/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/plotnine.yaml"
    params:
        bins=BIAS_BINS,
        title=lambda wc: f"{wc.id} ({BIAS_REFERENCE} -DMS reference): ",
    shell:
        """
        python3 workflow/scripts/three_prime_end_bias/fit_3prime_bias_model.py \
            --depth-plus {input.plus} --depth-minus {input.minus} --training {input.training} \
            --bins {params.bins} --title "{params.title}" \
            --model {output.model} --binned {output.binned} --plot {output.plot} >{log} 2>&1
        """


rule apply_3prime_bias_correction:
    """Divide every transcript's +DMS RT-stops (mRNA, ncRNA, positive control) by the fitted bias."""
    input:
        model=f"{BIAS_QC}/model.tsv",
        plus=lambda wc: get_raw_plus_rtsc(wc.id),
        minus=get_bias_reference_minus_rtsc,
        fasta=TRANSCRIPTOME_UPPER,
        training=BIAS_TRAINING_TSV,
    output:
        rtsc=f"{BIAS_TMP}/{{id}}/combined_plus_corrected.rtsc",
        binned=f"{BIAS_QC}/rtsc_metagene_before_after.tsv",
        plot=report(
            f"{BIAS_QC}/rtsc_metagene_before_after.png",
            category="3' bias correction",
            subcategory="{id}",
            caption="../report_captions/bias_rtsc_metagene_before_after.rst",
        ),
    log:
        "logs/apply_3prime_bias_correction/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/plotnine.yaml"
    params:
        bins=BIAS_BINS,
        title=lambda wc: f"{wc.id} ({BIAS_REFERENCE} -DMS reference): ",
    shell:
        """
        python3 workflow/scripts/three_prime_end_bias/apply_3prime_bias_correction.py \
            --model {input.model} --plus-rtsc {input.plus} --minus-rtsc {input.minus} \
            --fasta {input.fasta} --training {input.training} \
            --bins {params.bins} --title "{params.title}" \
            --output {output.rtsc} --binned {output.binned} --plot {output.plot} >{log} 2>&1
        """
