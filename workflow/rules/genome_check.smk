############################################################################
# Genome check (config key genome_check): where do the reads that did not align to the
# transcriptome come from? Bowtie2-unaligned reads of every sample are aligned to the genome with
# HISAT2 (splice-aware, end-to-end like the transcriptome alignment) and placed on the genome
# annotation by their 5' end. Outputs in {output_dir}/genome_check/:
#   <sample>.bam(.bai)  genome alignments of the sample's transcriptome-unaligned reads (IGV etc.)
#   summary.tsv         per sample: unaligned reads, % on the genome, % per category
#   categories.png      the categories as stacked bars
#   top_loci.tsv        the loci with the most reads, per-sample counts and the gene they hit
# See settings.smk (GENOME_CHECK_*) for the config keys.
############################################################################


def genome_check_ref(sample):
    return SAMPLE_GENOTYPE[sample] if GENOTYPE_MODE else "genome"


def genome_check_index_prefix(ref):
    return GENOME_CHECK_REFS[ref]["index"] or f"{GENOME_CHECK_TMP}/index/{ref}/genome"


def genome_check_index_flag(wildcards):
    """Built index (flag file), or nothing when an existing index was given."""
    ref = genome_check_ref(wildcards.sample)
    return [] if GENOME_CHECK_REFS[ref]["index"] else f"{GENOME_CHECK_TMP}/index/{ref}/index.done"


rule genome_check_hisat2_index:
    """HISAT2 index of one genome (skipped when genome_check_hisat2_index / hisat2_index is set)."""
    input:
        lambda wc: GENOME_CHECK_REFS[wc.ref]["genome"],
    output:
        touch(f"{GENOME_CHECK_TMP}/index/{{ref}}/index.done"),
    params:
        prefix=lambda wc: f"{GENOME_CHECK_TMP}/index/{wc.ref}/genome",
    threads: 8
    resources:
        mem_mb=16000,
    log:
        "logs/genome_check/hisat2_index_{ref}.log",
    conda:
        "../envs/hisat2.yaml"
    message:
        "Building the HISAT2 genome index for {wildcards.ref}"
    shell:
        "hisat2-build -p {threads} {input} {params.prefix} > {log} 2>&1"


rule genome_check_unaligned_reads:
    """Reads Bowtie2 left unaligned (flag 4) in one sample's merged transcriptome SAM."""
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam",
    output:
        temp(f"{GENOME_CHECK_TMP}/{{sample}}/unaligned.fastq.gz"),
    log:
        "logs/genome_check/unaligned_reads_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq -f 4 {input} 2> {log} | gzip -c > {output}"


rule genome_check_align:
    """
    Align one sample's transcriptome-unaligned reads to its genome: HISAT2 end-to-end
    (--no-softclip, as Bowtie2 aligned to the transcriptome), up to 10 alignments per read.
    Keeps a sorted BAM, HISAT2's summary, and the primary alignments as a small table.
    """
    input:
        fastq=f"{GENOME_CHECK_TMP}/{{sample}}/unaligned.fastq.gz",
        index=genome_check_index_flag,
    output:
        bam=f"{GENOME_CHECK_DIR}/{{sample}}.bam",
        bai=f"{GENOME_CHECK_DIR}/{{sample}}.bam.bai",
        summary=f"{GENOME_CHECK_TMP}/{{sample}}/hisat2_summary.txt",
        alignments=temp(f"{GENOME_CHECK_TMP}/{{sample}}/primary_alignments.tsv.gz"),
    params:
        index=lambda wc: genome_check_index_prefix(genome_check_ref(wc.sample)),
    threads: 4
    log:
        "logs/genome_check/align_{sample}.log",
    conda:
        "../envs/hisat2.yaml"
    message:
        "Aligning transcriptome-unaligned reads of {wildcards.sample} to the genome with HISAT2"
    shell:
        r"""
        exec 2> {log}
        hisat2 -p {threads} -x {params.index} -U {input.fastq} --no-softclip -k 10 --no-unal \
            --summary-file {output.summary} \
            | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        # primary alignments: read, flag, chrom, pos, cigar, NH
        samtools view -F 260 {output.bam} \
            | awk -F'\t' -v OFS='\t' '{{ nh = 1; for (i = 12; i <= NF; i++) if ($i ~ /^NH:i:/) nh = substr($i, 6);
                                        print $1, $2, $3, $4, $6, nh }}' \
            | gzip -c > {output.alignments}
        """


rule genome_check_annotation:
    """Exon / gene intervals of one genome's annotation, flagged by whether the transcript is in
    the transcriptome the reads were aligned to."""
    input:
        annotation=lambda wc: GENOME_CHECK_REFS[wc.ref]["annotation"] or [],
        transcriptome=TRANSCRIPTOME,
    output:
        f"{GENOME_CHECK_TMP}/annotation/{{ref}}_intervals.tsv.gz",
    params:
        annotation=lambda wc: GENOME_CHECK_REFS[wc.ref]["annotation"],
        prefix=lambda wc: f"{wc.ref}{GENOTYPE_SEP}" if GENOTYPE_MODE else "",
    log:
        "logs/genome_check/annotation_{ref}.log",
    conda:
        "../envs/plotnine.yaml"
    shell:
        "python3 workflow/scripts/genome_check/annotation_intervals.py --annotation '{params.annotation}' "
        "--transcriptome {input.transcriptome} --prefix '{params.prefix}' -o {output} > {log} 2>&1"


rule genome_check_classify:
    """Place one sample's genome alignments on the annotation (classify_genome_hits.py)."""
    input:
        alignments=f"{GENOME_CHECK_TMP}/{{sample}}/primary_alignments.tsv.gz",
        intervals=lambda wc: f"{GENOME_CHECK_TMP}/annotation/{genome_check_ref(wc.sample)}_intervals.tsv.gz",
    output:
        positions=f"{GENOME_CHECK_TMP}/{{sample}}/positions.tsv.gz",
        counts=f"{GENOME_CHECK_TMP}/{{sample}}/category_counts.tsv",
    params:
        organelles=GENOME_CHECK_ORGANELLES,
    log:
        "logs/genome_check/classify_{sample}.log",
    conda:
        "../envs/plotnine.yaml"
    shell:
        "python3 workflow/scripts/genome_check/classify_genome_hits.py --alignments {input.alignments} "
        "--intervals {input.intervals} --organelles '{params.organelles}' --sample {wildcards.sample} "
        "--positions {output.positions} --counts {output.counts} > {log} 2>&1"


rule genome_check_summary:
    """Per-sample summary, category plot and top loci over all samples."""
    input:
        counts=expand(f"{GENOME_CHECK_TMP}/{{sample}}/category_counts.tsv", sample=SAMPLES),
        positions=expand(f"{GENOME_CHECK_TMP}/{{sample}}/positions.tsv.gz", sample=SAMPLES),
        hisat2=expand(f"{GENOME_CHECK_TMP}/{{sample}}/hisat2_summary.txt", sample=SAMPLES),
        bowtie2=expand(f"{TMP}/output/se/{{sample}}/{{sample}}_mapping_log.txt", sample=SAMPLES),
    output:
        summary=f"{GENOME_CHECK_DIR}/summary.tsv",
        loci=f"{GENOME_CHECK_DIR}/top_loci.tsv",
        plot=f"{GENOME_CHECK_DIR}/categories.png",
    params:
        samples=SAMPLES,
        gap=GENOME_CHECK_LOCUS_GAP,
        top=GENOME_CHECK_TOP_LOCI,
    log:
        "logs/genome_check/summary.log",
    conda:
        "../envs/plotnine.yaml"
    shell:
        "python3 workflow/scripts/genome_check/summarize_genome_check.py --samples {params.samples} "
        "--counts {input.counts} --positions {input.positions} --hisat2-summaries {input.hisat2} "
        "--bowtie2-logs {input.bowtie2} --gap {params.gap} --top {params.top} "
        "--summary {output.summary} --loci {output.loci} --plot {output.plot} > {log} 2>&1"
