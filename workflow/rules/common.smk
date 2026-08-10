############################################################################
# Alignment and QC rules shared by all downstream tools
############################################################################

import os

# ---------------------------------------------------------------------------
# Transcriptome preparation
# ---------------------------------------------------------------------------
TRANSCRIPT_SOURCE = config.get("transcript_source", "file")
TMP = config.get("tmp_dir", "tmp")

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule below.
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"


def _source_transcriptome(wildcards):
    """Return the raw transcriptome FASTA before extra sequences are appended."""
    if TRANSCRIPT_SOURCE == "gffread":
        return f"{TMP}/resources/transcriptome_gffread_filtered.fa"
    return config["transcriptome"]


def _extra_transcriptome_fastas(wildcards):
    """
    Return extra FASTAs to append to the transcriptome before mapping:
    positive-control sequences (config key positive_control_fasta) and/or
    ncRNA sequences (config key ncrna_fasta). Either, both, or neither may
    be set.
    """
    extras = []
    if config.get("positive_control_fasta"):
        extras.append(config["positive_control_fasta"])
    if config.get("ncrna_fasta"):
        extras.append(config["ncrna_fasta"])
    return extras


if TRANSCRIPT_SOURCE == "gffread":
    rule agat_longest_isoform:
        """
        Use AGAT to retain only the longest isoform per gene so that
        gffread produces one transcript per locus.
        """
        input:
            gff=config["gff"]
        output:
            f"{TMP}/resources/longest_isoform.gff"
        conda:
            "../envs/agat.yaml"
        log:
            "logs/agat_longest_isoform/agat.log"
        message:
            "Filtering GFF to longest isoform per gene with AGAT"
        shell:
            "agat_sp_keep_longest_isoform.pl --gff {input.gff} --out {output} > {log} 2>&1"

    rule gffread_extract_transcripts:
        """
        Extract transcript sequences from the genome FASTA using gffread
        and the longest-isoform GFF produced by AGAT.
        """
        input:
            genome=config["genome"],
            gff=f"{TMP}/resources/longest_isoform.gff"
        output:
            f"{TMP}/resources/transcriptome_gffread.fa"
        conda:
            "../envs/gffread.yaml"
        log:
            "logs/gffread/gffread.log"
        message:
            "Extracting transcript sequences with gffread"
        shell:
            "gffread {input.gff} -g {input.genome} -w {output} > {log} 2>&1"

    rule filter_extracted_seqs:
        """
        Keep only those transcripts with a CDS
        """
        input:
            f"{TMP}/resources/transcriptome_gffread.fa"
        output:
            f"{TMP}/resources/transcriptome_gffread_filtered.fa"
        conda:
            "../envs/biopython.yaml"
        shell:
            "python3 workflow/scripts/filter_extracted_seqs.py --input {input} --output {output}"


rule prepare_transcriptome:
    """
    Assemble the canonical transcriptome FASTA at resources/transcriptome.fa.

    Base sequences come from _source_transcriptome (either the gffread
    pipeline output or the user-supplied FASTA, depending on
    transcript_source). Any extras from _extra_transcriptome_fastas
    (positive controls, ncRNAs) are appended after the base sequences.
    """
    input:
        transcriptome=_source_transcriptome,
        extras=_extra_transcriptome_fastas
    output:
        TRANSCRIPTOME
    log:
        "logs/prepare_transcriptome/prepare.log"
    message:
        "Preparing final transcriptome FASTA"
    shell:
        # Plain `cat` silently glues a record onto the next file's header if
        # the preceding file is missing its trailing newline. `awk 1`
        # reprints every line with a newline appended regardless of the
        # source file's ending, so concatenation is always newline-safe.
        "awk 1 {input.transcriptome} {input.extras} > {output} 2> {log}"


# ---------------------------------------------------------------------------
# Alignment (bowtie2 only -- the only aligner this pipeline actually runs
# end-to-end; see README for notes on extending it)
# ---------------------------------------------------------------------------

rule build_bowtie2_index:
    """
    Build a Bowtie2 index from the canonical transcriptome FASTA.
    """
    input:
        TRANSCRIPTOME
    output:
        TRANSCRIPTOME + ".1.bt2"
    conda:
        "../envs/structurefold.yaml"
    log:
        "logs/bowtie2_index/build.log"
    message:
        "Building Bowtie2 index for the transcriptome"
    shell:
        "bowtie2-build {input} {input} > {log} 2>&1"


rule align_with_bowtie2:
    """
    Align trimmed reads to the reference transcriptome using Bowtie2
    via StructureFold2 mapper
    """
    input:
        transcriptome_index=TRANSCRIPTOME + ".1.bt2",
        fastq1=f"{TMP}/reads/trimmed/{{sample}}_trimmed.fastq"
    output:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam",
        mapping_log=f"{TMP}/output/se/{{sample}}/{{sample}}_mapping_log.txt"
    params:
        directory=lambda wc: os.path.abspath(f"{TMP}/output/se/{wc.sample}"),
        abs_index=os.path.abspath(TRANSCRIPTOME),
        abs_log=lambda wc, output: os.path.abspath(output.mapping_log),
        script="scripts/StructureFold2/fastq_mapper.py",
        workdir=f"{workflow.basedir}/workflow"
    conda:
        "../envs/structurefold.yaml"
    message:
        "Aligning trimmed reads for sample {wildcards.sample} with Bowtie2"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.directory}

        cp {input.fastq1} {params.directory}/

        cd {params.directory}
        python {params.workdir}/{params.script} {params.abs_index} 2> {params.abs_log}
        cd {params.workdir}
        """


# ---------------------------------------------------------------------------
# Read preparation: rename -> trim
# ---------------------------------------------------------------------------

def get_fastqs(wildcards):
    """
    Return every raw fastq (one per sequencing run) for a sample, sorted by
    run, so rename_fastq can concatenate them into one canonical file.
    """
    rows = samples_runs.loc[samples_runs["sample"] == wildcards.sample]
    if rows.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in samplesheet")
    return rows.sort_values("run")["r1"].tolist()


def get_samples_by_condition(condition):
    rows = samples.loc[samples["condition"] == condition]
    if rows.empty:
        raise ValueError(f"No samples found for condition '{condition}'")
    return rows["sample"].tolist()


rule rename_fastq:
    """
    Concatenate every sequencing run's FASTQ for a sample (usually just one)
    into a single canonical renamed file for downstream rules.
    """
    input:
        get_fastqs
    output:
        f"{TMP}/reads/renamed/{{sample}}.fastq"
    log:
        "logs/rename_fastq/{sample}.log"
    message:
        "Renaming FASTQ for sample {wildcards.sample}"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


rule trim_reads:
    """
    Trim reads using the trimming script from StructureFold2
    """
    input:
        reads=f"{TMP}/reads/renamed/{{sample}}.fastq"
    output:
        fastq=f"{TMP}/reads/trimmed/{{sample}}_trimmed.fastq"
    threads: 4
    conda:
        "../envs/structurefold.yaml"
    params:
        tmpdir=f"{TMP}/structurefold2/{{sample}}",
        workdir=f"{workflow.basedir}/workflow"
    log:
        "logs/trim_reads/{sample}.log"
    message:
        "Trimming reads for sample {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.tmpdir}
        cp {input.reads} {params.tmpdir}/

        cd {params.tmpdir}
        python {params.workdir}/scripts/StructureFold2/fastq_trimmer.py
        cd {params.workdir}
        cd ..
        mkdir -p $(dirname {output.fastq})
        mv {params.tmpdir}/{wildcards.sample}_trimmed.fastq {output.fastq}
        """


# ---------------------------------------------------------------------------
# SAM/BAM post-processing
# ---------------------------------------------------------------------------

rule filter_sam:
    """
    Filter mapped reads by mismatch count using the StructureFold2 sam_filter
    script.
    """
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"
    output:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
    conda:
        "../envs/structurefold.yaml"
    params:
        sample_dir=lambda wc: os.path.abspath(f"{TMP}/output/se/{wc.sample}"),
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold2/sam_filter.py",
        # Must be absolute: the shell below cd's into sample_dir first, so a
        # relative log path here would resolve against that directory
        # instead of the caller's cwd (only visible when tmp_dir is itself
        # a relative path, e.g. this test suite's tests/tmp).
        log=lambda wc: os.path.abspath(f"{TMP}/output/se/{wc.sample}/{wc.sample}_trimmed_mapped_filtered.log"),
    shell:
        """
        cd {params.sample_dir} &&
        python {params.workdir}/{params.script} -max_mismatch 4 -logname {params.log}
        """


rule add_header_to_filtered_sam:
    """
    sam_filter.py strips SAM headers from the filtered output.
    This rule copies the filtered SAM and prepends the header lines
    from the original mapped SAM so that downstream tools (e.g.
    samtools depth) can parse it correctly.
    """
    input:
        filtered=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
        original=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"
    output:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_with_header.sam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input.original} > {output}
        grep -v '^@' {input.filtered} >> {output}
        """


rule filtered_sam_to_sorted_bam:
    """
    Convert quality-filtered SAM (with header restored) to a sorted, indexed BAM.
    """
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_with_header.sam"
    output:
        bam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_sorted.bam",
        bai=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_sorted.bam.bai",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -bS {input} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """


# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

rule fastqc:
    """
    Run FastQC on trimmed reads for each sample
    """
    input:
        f"{TMP}/reads/trimmed/{{sample}}_trimmed.fastq"
    output:
        html=f"{config['output_dir']}/qc/fastqc/{{sample}}_trimmed_fastqc.html",
        zip=f"{config['output_dir']}/qc/fastqc/{{sample}}_trimmed_fastqc.zip"
    params:
        outdir=f"{config['output_dir']}/qc/fastqc"
    log:
        "logs/fastqc/{sample}.log"
    message:
        "Running FastQC on trimmed reads for sample {wildcards.sample}"
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {params.outdir} && fastqc {input} --outdir {params.outdir} > {log} 2>&1"


rule multiqc:
    """
    Aggregate FastQC reports across all samples with MultiQC
    """
    input:
        expand(f"{config['output_dir']}/qc/fastqc/{{sample}}_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        f"{config['output_dir']}/qc/multiqc/multiqc_report.html"
    params:
        indir=f"{config['output_dir']}/qc/fastqc",
        outdir=f"{config['output_dir']}/qc/multiqc"
    log:
        "logs/multiqc.log"
    message:
        "Running MultiQC"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc {params.indir} --outdir {params.outdir} > {log} 2>&1"


rule alignment_stats:
    """
    Extract per-sample alignment statistics (Bowtie2 log format): total
    reads, unique mapped %, multi-mapped %, and overall alignment %.
    """
    input:
        mapping_log=f"{TMP}/output/se/{{sample}}/{{sample}}_mapping_log.txt"
    output:
        f"{config['output_dir']}/qc/alignment_stats/{{sample}}_stats.tsv"
    log:
        "logs/alignment_stats/{sample}.log"
    message:
        "Extracting alignment stats for sample {wildcards.sample}"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})

        awk -v samp="{wildcards.sample}" '
            /reads; of these:/          {{ total = $1 }}
            /aligned exactly 1 time/   {{ gsub(/[()%]/, "", $2); unique = $2 }}
            /aligned >1 times/         {{ gsub(/[()%]/, "", $2); multi  = $2 }}
            /overall alignment rate/   {{ gsub(/%/, "", $1);     overall = $1 }}
            END {{
                print "sample\ttotal_reads\tunique_mapped_pct\tmulti_mapped_pct\toverall_alignment_pct"
                printf "%s\t%d\t%.2f\t%.2f\t%.2f\n", samp, total, unique, multi, overall
            }}
        ' {input.mapping_log} > {output} 2>> {log}
        """


rule aggregate_alignment_stats:
    """
    Combine per-sample alignment stats into a single summary table.
    """
    input:
        expand(
            f"{config['output_dir']}/qc/alignment_stats/{{sample}}_stats.tsv",
            sample=SAMPLES
        )
    output:
        f"{config['output_dir']}/qc/alignment_stats_summary.tsv"
    message:
        "Aggregating alignment stats across all samples"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """


# ---------------------------------------------------------------------------
# Positive-control alignment-% QC (any positive control, not just p4p6 --
# see workflow/rules/positive_control.smk for the p4p6-only structure plot)
# ---------------------------------------------------------------------------

if config.get("positive_control_fasta"):

    rule alignment_pct_positive_control:
        """
        Calculate per-sample alignment percentage to the positive-control
        sequence, identified by config['positive_control_name'] (must match
        that sequence's FASTA header exactly, e.g. "p4p6"). Uses
        samtools view -F 4 to select only mapped reads, then counts those
        mapping to the positive control vs total mapped reads.
        """
        input:
            sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"
        output:
            f"{config['output_dir']}/qc/positive_control_alignment/{{sample}}_pct.tsv"
        params:
            pc_name=config.get("positive_control_name", ""),
        log:
            "logs/positive_control_alignment/{sample}.log"
        message:
            "Calculating positive-control alignment percentage for sample {wildcards.sample}"
        conda:
            "../envs/structurefold.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output})

            samtools view -F 4 {input.sam} | \
            awk -v samp="{wildcards.sample}" -v pc="{params.pc_name}" '
                {{
                    total++
                    if ($3 == pc) pc_reads++
                }}
                END {{
                    pct = (total > 0) ? pc_reads / total * 100 : 0
                    print "sample\ttotal_mapped\tpositive_control_reads\tpositive_control_pct"
                    printf "%s\t%d\t%d\t%.4f\n", samp, total, pc_reads, pct
                }}
            ' > {output} 2>> {log}
            """


    rule aggregate_positive_control_alignment:
        """
        Combine per-sample positive-control alignment percentages into a
        single summary table.
        """
        input:
            expand(
                f"{config['output_dir']}/qc/positive_control_alignment/{{sample}}_pct.tsv",
                sample=SAMPLES
            )
        output:
            f"{config['output_dir']}/qc/positive_control_alignment_summary.tsv"
        message:
            "Aggregating positive-control alignment percentages across all samples"
        shell:
            """
            head -1 {input[0]} > {output}
            for f in {input}; do tail -n +2 "$f" >> {output}; done
            """
