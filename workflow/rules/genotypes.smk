############################################################################
# Two genotypes (config key genotypes): species, or individuals of one species
#
# Each genotype has its own genome + GFF/GTF. Its transcripts are extracted with gffread (same
# keep_longest_isoform / cds_only options as the single-genome route) and renamed
# "<genotype>__<transcript_id>". Reads map only to their own genotype (a per-genotype Bowtie2
# index of its transcripts + the shared positive control / ncRNAs), and every sample's .rtsc lists
# only its genotype's transcripts. Downstream rules use the combined transcriptome (both genotypes
# + the shared extras), so per-ID reactivity, coverage, specificity, bias correction etc. run
# unchanged.
#
# Matching: transcripts with the same ID in both genotypes, or the rows of transcript_matching.
# match_genotype_transcripts maps positions between each pair (1:1 for equal-length, SNP-only pairs;
# a pairwise alignment otherwise). That map feeds the optional cross-genotype -DMS pooling
# (pool_minus_across_genotypes) and the reactivity comparison ({output_dir}/genotype_comparison/).
############################################################################

GENOTYPE_OUT = f"{config['output_dir']}/genotype_comparison"
GENOTYPE_A, GENOTYPE_B = GENOTYPE_NAMES
CROSS_GENOTYPE_MIN_IDENTITY = float(config.get("cross_genotype_min_identity", 0.9))
GENOTYPE_COMPARE_MIN_POSITIONS = int(config.get("genotype_comparison_min_positions", 20))


wildcard_constraints:
    genotype="|".join(GENOTYPE_NAMES),


def _genotype_gff(wildcards):
    if KEEP_LONGEST_ISOFORM:
        return f"{GENOTYPE_DIR}/{wildcards.genotype}/longest_isoform.gff"
    return GENOTYPES[wildcards.genotype]["gff"]


def _genotype_source(wildcards):
    """The genotype's gffread transcripts (CDS-filtered when cds_only), original IDs."""
    suffix = "_filtered" if CDS_ONLY else ""
    return f"{GENOTYPE_DIR}/{wildcards.genotype}/transcriptome_gffread{suffix}.fa"


if KEEP_LONGEST_ISOFORM:

    rule agat_longest_isoform_genotype:
        """Longest isoform per gene for one genotype's GFF (keep_longest_isoform)."""
        input:
            gff=lambda wc: GENOTYPES[wc.genotype]["gff"],
        output:
            f"{GENOTYPE_DIR}/{{genotype}}/longest_isoform.gff",
        conda:
            "../envs/agat.yaml"
        log:
            "logs/agat_longest_isoform/{genotype}.log",
        shell:
            "agat_sp_keep_longest_isoform.pl --gff {input.gff} --out {output} > {log} 2>&1"


rule gffread_extract_transcripts_genotype:
    """Extract one genotype's transcript sequences from its own genome."""
    input:
        genome=lambda wc: GENOTYPES[wc.genotype]["genome"],
        gff=_genotype_gff,
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/transcriptome_gffread.fa",
    conda:
        "../envs/gffread.yaml"
    log:
        "logs/gffread/{genotype}.log",
    shell:
        "gffread {input.gff} -g {input.genome} -w {output} > {log} 2>&1"


if CDS_ONLY:

    rule filter_extracted_seqs_genotype:
        """Keep only one genotype's transcripts with a CDS (cds_only)."""
        input:
            f"{GENOTYPE_DIR}/{{genotype}}/transcriptome_gffread.fa",
        output:
            f"{GENOTYPE_DIR}/{{genotype}}/transcriptome_gffread_filtered.fa",
        conda:
            "../envs/biopython.yaml"
        log:
            "logs/filter_extracted_seqs/{genotype}.log",
        shell:
            "python3 workflow/scripts/transcriptome/filter_extracted_seqs.py --input {input} --output {output} > {log} 2>&1"


rule prefix_genotype_transcripts:
    """Rename one genotype's transcripts <genotype>__<id> (rest of the header kept)."""
    input:
        _genotype_source,
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/transcripts.fa",
    params:
        prefix=lambda wc: f"{wc.genotype}{GENOTYPE_SEP}",
    log:
        "logs/prefix_genotype_transcripts/{genotype}.log",
    shell:
        """awk -v p="{params.prefix}" '/^>/ {{ sub(/^>/, ">" p) }} {{ print }}' {input} > {output} 2> {log}"""


rule combine_genotype_transcripts:
    """Both genotypes' renamed transcripts: the base of the combined transcriptome (TRANSCRIPTOME)."""
    input:
        expand(f"{GENOTYPE_DIR}/{{genotype}}/transcripts.fa", genotype=GENOTYPE_NAMES),
    output:
        f"{GENOTYPE_DIR}/all_genotypes_transcripts.fa",
    log:
        "logs/combine_genotype_transcripts.log",
    shell:
        "awk 1 {input} > {output} 2> {log}"


rule genotype_mapping_fasta:
    """What one genotype's reads map to: its own transcripts + the shared positive control / ncRNAs."""
    input:
        transcripts=f"{GENOTYPE_DIR}/{{genotype}}/transcripts.fa",
        extras=_extra_transcriptome_fastas,
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/mapping.fa",
    params:
        awk="workflow/scripts/transcriptome/concat_fasta_unique.awk",
    log:
        "logs/genotype_mapping_fasta/{genotype}.log",
    shell:
        "awk -f {params.awk} {input.transcripts} {input.extras} > {output} 2> {log}"


rule uppercase_genotype_mapping_fasta:
    """Uppercase copy of one genotype's mapping FASTA (sam_to_rtsc writes RT-stops over it)."""
    input:
        f"{GENOTYPE_DIR}/{{genotype}}/mapping.fa",
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/mapping_upper.fa",
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/uppercase_genotype_mapping_fasta/{genotype}.log",
    shell:
        "python3 workflow/scripts/transcriptome/uppercase_fasta.py --input {input} --output {output} > {log} 2>&1"


rule build_bowtie2_index_genotype:
    """Bowtie2 index of one genotype's mapping FASTA."""
    input:
        f"{GENOTYPE_DIR}/{{genotype}}/mapping.fa",
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/mapping.fa.1.bt2",
    conda:
        "../envs/bowtie2.yaml"
    threads: 4
    log:
        "logs/bowtie2_index/{genotype}.log",
    shell:
        "bowtie2-build --threads {threads} {input} {input} > {log} 2>&1"


rule match_genotype_transcripts:
    """Matched transcript pairs between the genotypes and the position map within each pair."""
    input:
        a=f"{GENOTYPE_DIR}/{GENOTYPE_A}/transcripts.fa",
        b=f"{GENOTYPE_DIR}/{GENOTYPE_B}/transcripts.fa",
        matching=[TRANSCRIPT_MATCHING] if TRANSCRIPT_MATCHING else [],
    output:
        GENOTYPE_LIFTOVER,
    log:
        "logs/match_genotype_transcripts.log",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    conda:
        "../envs/genotypes.yaml"
    params:
        matching=f"--matching {TRANSCRIPT_MATCHING}" if TRANSCRIPT_MATCHING else "",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 workflow/scripts/genotypes/match_transcripts.py \
            --fasta-a {input.a} --fasta-b {input.b} --genotype-a {GENOTYPE_A} --genotype-b {GENOTYPE_B} \
            --sep {GENOTYPE_SEP} {params.matching} --threads {threads} --output {output} >{log} 2>&1
        """


rule depth_matching_key_map:
    """name -> key for depth matching: both members of a matched pair get the first genotype's name,
    so each matched transcript is thinned to the same +DMS total in every library of both genotypes."""
    input:
        GENOTYPE_LIFTOVER,
    output:
        f"{GENOTYPE_DIR}/depth_matching_key_map.tsv",
    log:
        "logs/depth_matching_key_map.log",
    shell:
        """
        awk -F'\t' -v OFS='\t' 'NR == 1 {{ print "name", "key"; next }} {{ print $1, $1; print $2, $1 }}' \
            {input} > {output} 2> {log}
        """


rule combine_rtsc_minus_pooled_genotype:
    """Pool the -DMS RT-stops of every -DMS sample of one genotype (all temperatures and replicates).
    Without pool_minus_across_genotypes, each ID's reactivity, depth-matched reactivity and 3' bias
    before/after plot use its own genotype's pool, so the shared positive control / ncRNAs / lncRNAs
    also get only that genotype's -DMS."""
    input:
        lambda wc: [f"{TMP}/output/se/{s}/{s}_trimmed_mapped_filtered.rtsc"
                    for s in get_samples_by_condition("minus") if SAMPLE_GENOTYPE[s] == wc.genotype],
    output:
        f"{TMP}/output/se/pooled/{{genotype}}/combined_minus.rtsc",
    log:
        "logs/combine_rtsc_minus_pooled_genotype/{genotype}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_combine.py",
        tmpdir=f"{TMP}/structurefold2/pooled_{{genotype}}_combined_minus",
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


if POOL_MINUS_ACROSS_GENOTYPES:

    rule pool_minus_across_genotypes:
        """
        Every reactivity's -DMS: the -DMS pooled within each genotype, plus, for each matched pair,
        the other genotype's stops at matching positions. Positions come from
        cross_genotype_coordinate_map when set (variant-derived; every pair with a shared ID, no
        identity cutoff), else from the alignment blocks (pairs with identity >=
        cross_genotype_min_identity).
        """
        input:
            minus=f"{TMP}/output/se/pooled/combined_minus.rtsc",
            liftover=GENOTYPE_LIFTOVER,
            coordinate_map=[CROSS_GENOTYPE_COORDINATE_MAP] if CROSS_GENOTYPE_COORDINATE_MAP else [],
        output:
            CROSS_GENOTYPE_MINUS,
        log:
            "logs/pool_minus_across_genotypes.log",
        resources:
            mem_mb=16000,
        conda:
            "../envs/plotnine.yaml"
        params:
            min_identity=CROSS_GENOTYPE_MIN_IDENTITY,
            map_args=(f"--coordinate-map {CROSS_GENOTYPE_COORDINATE_MAP} --map-ref-genotype {CROSS_GENOTYPE_MAP_REF}"
                      if CROSS_GENOTYPE_COORDINATE_MAP else ""),
        shell:
            """
            python3 workflow/scripts/genotypes/pool_minus_across_genotypes.py --minus {input.minus} \
                --liftover {input.liftover} --min-identity {params.min_identity} {params.map_args} \
                --output {output} >{log} 2>&1
            """


def final_reactivity_csv(id):
    """The reactivity a run ends with: depth-matched when depth_matching is on (from the 3' bias-
    corrected stops when that is on too), else the standard reactivity.csv."""
    if DEPTH_MATCHING:
        return f"{config['output_dir']}/depth_matched/{id}/reactivity.csv"
    return f"{config['output_dir']}/{id}/reactivity.csv"


rule compare_genotypes:
    """Per-transcript reactivity agreement between the genotypes over matched transcripts."""
    input:
        react=[final_reactivity_csv(id) for id in IDS],
        liftover=GENOTYPE_LIFTOVER,
    output:
        per_transcript=f"{GENOTYPE_OUT}/reactivity_comparison.tsv",
        per_position=f"{GENOTYPE_OUT}/reactivity_comparison_positions.tsv.gz",
        summary=f"{GENOTYPE_OUT}/reactivity_comparison_summary.tsv",
        plot=f"{GENOTYPE_OUT}/reactivity_comparison.png",
    log:
        "logs/compare_genotypes.log",
    resources:
        mem_mb=48000,
        runtime=240,
    conda:
        "../envs/plotnine.yaml"
    params:
        labels=IDS,
        genotypes=[ID_GENOTYPE[id] for id in IDS],
        min_positions=GENOTYPE_COMPARE_MIN_POSITIONS,
    shell:
        """
        python3 workflow/scripts/genotypes/compare_genotypes.py \
            --inputs {input.react} --labels {params.labels} --genotypes {params.genotypes} \
            --genotype-a {GENOTYPE_A} --genotype-b {GENOTYPE_B} --liftover {input.liftover} \
            --min-positions {params.min_positions} \
            --per-transcript {output.per_transcript} --per-position {output.per_position} \
            --summary {output.summary} --plot {output.plot} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Per-genotype replicate QC: the two genotypes' transcripts have different names, so RT-stop
# correlation and the covered-transcript UpSet compare replicates within each genotype.
# ---------------------------------------------------------------------------
rule rtsc_stop_correlation_genotype:
    """RT-stop correlation across one genotype's replicates (+DMS or -DMS)."""
    input:
        lambda wc: expand(
            f"{TMP}/output/se/{{id}}/combined_{wc.condition}.rtsc",
            id=(GENOTYPE_PLUS_IDS if wc.condition == "plus" else GENOTYPE_MINUS_IDS)[wc.genotype],
        ),
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation_{{condition}}_{{genotype}}.csv",
    wildcard_constraints:
        condition="plus|minus",
    log:
        "logs/rtsc_stop_correlation/{condition}_{genotype}.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        fasta=lambda wc: f"{GENOTYPE_DIR}/{wc.genotype}/mapping_upper.fa",
    shell:
        """
        python3 {params.workdir}/scripts/StructureFold3/rtsc_correlation.py {input} -spec AC \
            -fasta {params.fasta} -name {output} >{log} 2>&1
        """


rule plot_rt_stop_correlation_genotype:
    """Pairwise replicate RT-stop correlation plot for one genotype."""
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation_{{condition}}_{{genotype}}.csv",
    output:
        f"{config['output_dir']}/qc/rt_stop_replicate_correlation_{{condition}}_{{genotype}}.png",
    wildcard_constraints:
        condition="plus|minus",
    log:
        "logs/plot_rt_stop_correlation/{condition}_{genotype}.log",
    conda:
        "../envs/plotnine.yaml"
    shell:
        "python3 workflow/scripts/qc_plots/plot_rt_stop_correlation.py --input {input} --output {output} >{log} 2>&1"


rule upset_covered_transcripts_genotype:
    """UpSet plot of covered transcripts across one genotype's replicates."""
    input:
        lambda wc: expand(f"{config['output_dir']}/{{id}}/coverage.csv", id=GENOTYPE_PLUS_IDS[wc.genotype]),
    output:
        f"{config['output_dir']}/qc/covered_transcripts_upset_replicates_{{genotype}}.png",
    log:
        "logs/upset_covered_transcripts/{genotype}.log",
    conda:
        "../envs/upset.yaml"
    params:
        labels=lambda wc: GENOTYPE_PLUS_IDS[wc.genotype],
    shell:
        """
        python3 workflow/scripts/qc_plots/plot_covered_transcripts_upset.py --coverages {input} \
            --labels {params.labels} --threshold 1.0 --output {output} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Per-genotype annotation: protein-coding list (3' bias training) and per-position UTR/CDS labels,
# each made on the genotype's own IDs and then renamed <genotype>__<id>.
# ---------------------------------------------------------------------------
if BIAS_CORRECTION:

    rule bias_protein_coding_transcripts_genotype:
        """One genotype's protein-coding transcripts, renamed <genotype>__<id>."""
        input:
            annotation=lambda wc: genotype_annotation(wc.genotype),
            fasta=_genotype_source,
        output:
            f"{BIAS_TMP}/protein_coding_transcripts_{{genotype}}.txt",
        log:
            "logs/bias_protein_coding_transcripts/{genotype}.log",
        conda:
            "../envs/metagene.yaml"
        params:
            prefix=lambda wc: f"{wc.genotype}{GENOTYPE_SEP}",
        shell:
            """
            python3 workflow/scripts/three_prime_end_bias/list_protein_coding_transcripts.py \
                --annotation {input.annotation} --fasta {input.fasta} --output {output}.tmp >{log} 2>&1
            awk -v p="{params.prefix}" 'NF {{ print p $0 }}' {output}.tmp > {output} && rm {output}.tmp
            """

    rule bias_protein_coding_transcripts_combined:
        """Both genotypes' protein-coding transcripts (3' bias training candidates)."""
        input:
            expand(f"{BIAS_TMP}/protein_coding_transcripts_{{genotype}}.txt", genotype=GENOTYPE_NAMES),
        output:
            f"{BIAS_TMP}/protein_coding_transcripts.txt",
        log:
            "logs/bias_protein_coding_transcripts/combined.log",
        shell:
            "cat {input} > {output} 2> {log}"


rule extract_cds_and_ref_transcriptome_genotype:
    """gffread CDS (-x) and reference transcripts (-w) of one genotype."""
    input:
        genome=lambda wc: GENOTYPES[wc.genotype]["genome"],
        annotation=lambda wc: genotype_annotation(wc.genotype),
    output:
        cds=f"{GENOTYPE_DIR}/{{genotype}}/cds.fa",
        ref_tx=f"{GENOTYPE_DIR}/{{genotype}}/ref_transcriptome.fa",
    log:
        "logs/extract_cds_and_ref_transcriptome/{genotype}.log",
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -g {input.genome} {input.annotation} -x {output.cds} -w {output.ref_tx} >{log} 2>&1"


rule annotate_transcript_positions_genotype:
    """Per-position 5UTR/CDS/3UTR labels for one genotype's transcripts, renamed <genotype>__<id>."""
    input:
        ref_tx=f"{GENOTYPE_DIR}/{{genotype}}/ref_transcriptome.fa",
        cds=f"{GENOTYPE_DIR}/{{genotype}}/cds.fa",
        tx=_genotype_source,
    output:
        f"{GENOTYPE_DIR}/{{genotype}}/transcript_position_annotations.csv",
    log:
        "logs/annotate_transcript_positions/{genotype}.log",
    conda:
        "../envs/biopython.yaml"
    params:
        prefix=lambda wc: f"{wc.genotype}{GENOTYPE_SEP}",
    shell:
        """
        python3 workflow/scripts/annotation/annotate_transcript_positions.py --ref_tx_fasta {input.ref_tx} \
            --cds_fasta {input.cds} --tx_fasta {input.tx} --output {output}.tmp 2>{log}
        awk -v p="{params.prefix}" 'NR == 1 {{ print; next }} {{ print p $0 }}' {output}.tmp > {output} && rm {output}.tmp
        """


rule combine_transcript_position_annotations:
    """Both genotypes' per-position annotations in one table."""
    input:
        expand(f"{GENOTYPE_DIR}/{{genotype}}/transcript_position_annotations.csv", genotype=GENOTYPE_NAMES),
    output:
        f"{config['output_dir']}/transcript_position_annotations.csv",
    log:
        "logs/combine_transcript_position_annotations.log",
    shell:
        "awk 'FNR == 1 && NR != 1 {{ next }} {{ print }}' {input} > {output} 2> {log}"
