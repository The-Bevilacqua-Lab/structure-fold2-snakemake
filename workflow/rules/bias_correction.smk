############################################################################
# Optional 3' coverage-bias correction of the +DMS RT-stops
#
# Enabled by config['3_prime_bias_correction']. Per reactivity ID:
#   1. training set (one for the whole run, shared by every ID): the top
#      3_prime_bias_n_transcripts protein-coding transcripts ranked by their
#      lowest +DMS RT-stop coverage across all IDs;
#   2. fit log(+DMS depth / -DMS depth) ~ distance_from_3' + relative_position
#      + interaction on the training set's per-nucleotide read depth, with the
#      -DMS reference either pooled across every -DMS sample or the ID's own
#      -DMS samples (3_prime_bias_reference);
#   3. divide every protein-coding transcript's +DMS RT-stops by the fitted
#      bias. get_plus_rtsc_for_reactivity (helpers.smk) then hands the
#      corrected .rtsc to the reactivity rules.
# QC lands in {output_dir}/{id}/3prime_bias_correction/; the shared training
# set in {output_dir}/qc/3prime_bias_correction/training_transcripts.tsv.
############################################################################

BIAS_QC = f"{config['output_dir']}/{{id}}/3prime_bias_correction"
BIAS_TRAINING_TSV = f"{config['output_dir']}/qc/3prime_bias_correction/training_transcripts.tsv"


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
        python3 workflow/scripts/list_protein_coding_transcripts.py \
            --annotation {input.annotation} --fasta {input.fasta} --output {output} >{log} 2>&1
        """


rule bias_training_transcripts:
    """
    One training set shared by every ID: the top-N protein-coding transcripts
    ranked by their lowest +DMS RT-stop coverage across all IDs.
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
        n=BIAS_N_TRANSCRIPTS,
        labels=IDS,
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {output.bed}) \
            && python3 workflow/scripts/select_bias_training_transcripts.py \
                --rtsc {input.rtsc} --labels {params.labels} --fasta {input.fasta} \
                --protein-coding {input.protein_coding} --n {params.n} \
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
    shell:
        """
        samtools depth -a -b {input.bed} {input.bams} 2>{log} \
            | awk -F'\\t' -v OFS='\\t' '{{s = 0; for (i = 3; i <= NF; i++) s += $i; print $1, $2, s}}' \
            | gzip >{output}
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
        "../envs/metagene.yaml"
    params:
        bins=BIAS_BINS,
        title=lambda wc: f"{wc.id} ({BIAS_REFERENCE} -DMS reference): ",
    shell:
        """
        python3 workflow/scripts/fit_3prime_bias_model.py \
            --depth-plus {input.plus} --depth-minus {input.minus} --training {input.training} \
            --bins {params.bins} --title "{params.title}" \
            --model {output.model} --binned {output.binned} --plot {output.plot} >{log} 2>&1
        """


rule apply_3prime_bias_correction:
    """Divide every protein-coding transcript's +DMS RT-stops by the fitted bias."""
    input:
        model=f"{BIAS_QC}/model.tsv",
        plus=lambda wc: get_raw_plus_rtsc(wc.id),
        minus=get_bias_reference_minus_rtsc,
        fasta=TRANSCRIPTOME_UPPER,
        protein_coding=f"{BIAS_TMP}/protein_coding_transcripts.txt",
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
        "../envs/metagene.yaml"
    params:
        bins=BIAS_BINS,
        title=lambda wc: f"{wc.id} ({BIAS_REFERENCE} -DMS reference): ",
    shell:
        """
        python3 workflow/scripts/apply_3prime_bias_correction.py \
            --model {input.model} --plus-rtsc {input.plus} --minus-rtsc {input.minus} \
            --fasta {input.fasta} --protein-coding {input.protein_coding} --training {input.training} \
            --bins {params.bins} --title "{params.title}" \
            --output {output.rtsc} --binned {output.binned} --plot {output.plot} >{log} 2>&1
        """
