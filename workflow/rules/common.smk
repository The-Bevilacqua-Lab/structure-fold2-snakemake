############################################################################
# Alignment and QC rules shared by all downstream tools
############################################################################

import os

# ---------------------------------------------------------------------------
# Transcriptome preparation
# ---------------------------------------------------------------------------
TRANSCRIPT_SOURCE = config.get("transcript_source", "file")
TMP = config.get("tmp_dir", "tmp")

# Bases excluded from the 3' end of each transcript when generating the
# 2-8% normalization scale (see config.yaml's trim3 comment).
TRIM3 = config.get("trim3", 0)

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule below.
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"

# Uppercase-only copy of TRANSCRIPTOME, produced by uppercase_transcriptome
# below. Used (instead of TRANSCRIPTOME) by every rule downstream of
# mapping that checks base identity -- see that rule's docstring for why.
TRANSCRIPTOME_UPPER = f"{TMP}/resources/transcriptome_upper.fa"


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
        log:
            "logs/filter_extracted_seqs/filter.log"
        shell:
            "python3 workflow/scripts/filter_extracted_seqs.py --input {input} --output {output} > {log} 2>&1"


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


rule uppercase_transcriptome:
    """
    Write an uppercase-only copy of TRANSCRIPTOME (TRANSCRIPTOME_UPPER).

    Bowtie2 alignment is case-insensitive, so a source transcriptome FASTA
    with soft-masked (lowercase) bases still maps correctly -- but every
    base-identity check downstream of mapping (reactivity/specificity "is
    this position an A or C" tests, UTR/CDS annotation) compares
    case-sensitively against an uppercase specificity set, so those
    lowercase bases would otherwise be silently skipped. Kept as a separate
    file (rather than uppercasing TRANSCRIPTOME itself) so the already-built
    Bowtie2 index/alignments never need to be redone for this.
    """
    input:
        TRANSCRIPTOME
    output:
        TRANSCRIPTOME_UPPER
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/uppercase_transcriptome/uppercase.log"
    message:
        "Uppercasing transcriptome for downstream base-identity checks"
    shell:
        "python3 workflow/scripts/uppercase_fasta.py --input {input} --output {output} > {log} 2>&1"


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
        f"{TRANSCRIPTOME}.1.bt2"
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
    Align one sequencing run's trimmed reads to the reference transcriptome
    using Bowtie2 via StructureFold2 mapper. Runs per (sample, run) -- see
    the "Read preparation" section below for why -- landing in a scratch
    se_runs/ directory rather than the final per-sample se/ directory;
    combine_sam_runs merges every run's SAM back into the single per-sample
    SAM everything downstream expects.
    """
    input:
        transcriptome_index=f"{TRANSCRIPTOME}.1.bt2",
        fastq1=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"
    output:
        sam=f"{TMP}/output/se_runs/{{sample}}_{{run}}/{{sample}}_{{run}}_trimmed_mapped.sam",
        mapping_log=f"{TMP}/output/se_runs/{{sample}}_{{run}}/{{sample}}_{{run}}_mapping_log.txt"
    params:
        directory=lambda wc: os.path.abspath(f"{TMP}/output/se_runs/{wc.sample}_{wc.run}"),
        abs_index=os.path.abspath(TRANSCRIPTOME),
        abs_log=lambda wc, output: os.path.abspath(output.mapping_log),
        script="scripts/StructureFold2/fastq_mapper.py",
        workdir=f"{workflow.basedir}/workflow"
    log:
        "logs/align_with_bowtie2/{sample}_{run}.log"
    conda:
        "../envs/structurefold.yaml"
    message:
        "Aligning trimmed reads for sample {wildcards.sample}, run {wildcards.run} with Bowtie2"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        mkdir -p {params.directory}

        cp {input.fastq1} {params.directory}/

        cd {params.directory}
        python {params.workdir}/{params.script} {params.abs_index} 2> {params.abs_log}
        cd {params.workdir}
        """


def get_run_sams_for_sample(wildcards):
    runs = get_runs_for_sample(wildcards)
    return [
        f"{TMP}/output/se_runs/{wildcards.sample}_{r}/{wildcards.sample}_{r}_trimmed_mapped.sam"
        for r in runs
    ]


def get_run_mapping_logs_for_sample(wildcards):
    runs = get_runs_for_sample(wildcards)
    return [
        f"{TMP}/output/se_runs/{wildcards.sample}_{r}/{wildcards.sample}_{r}_mapping_log.txt"
        for r in runs
    ]


rule combine_sam_runs:
    """
    Merge one sample's per-run SAMs (each independently trimmed and
    aligned -- see rename_fastq/trim_reads/align_with_bowtie2 above) into
    the single per-sample SAM every downstream rule expects. All per-run
    SAMs share an identical header (same Bowtie2 index for every run), so
    it's safe to take the header from the first run and append every run's
    alignment records after it -- same approach as add_header_to_filtered_sam
    below. Per-run mapping logs are combined the same way (raw counts
    summed, percentages recomputed -- see combine_mapping_logs.py) into one
    Bowtie2-log-formatted {sample}_mapping_log.txt, so alignment_stats needs
    no changes.
    """
    input:
        sams=get_run_sams_for_sample,
        logs=get_run_mapping_logs_for_sample,
    output:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam",
        mapping_log=f"{TMP}/output/se/{{sample}}/{{sample}}_mapping_log.txt",
    conda:
        "../envs/samtools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
    log:
        "logs/combine_sam_runs/{sample}.log"
    message:
        "Combining sequencing runs of sample {wildcards.sample} into one SAM"
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1
        mkdir -p $(dirname {output.sam})

        first=1
        for f in {input.sams}; do
            if [ "$first" -eq 1 ]; then
                samtools view -H "$f" > {output.sam}
                first=0
            fi
            grep -v '^@' "$f" >> {output.sam}
        done

        python3 {params.workdir}/scripts/combine_mapping_logs.py {input.logs} -o {output.mapping_log}
        """


# ---------------------------------------------------------------------------
# Read preparation: rename -> trim
#
# Each sequencing run of a sample is renamed, trimmed, and aligned
# separately (rename_fastq / trim_reads / align_with_bowtie2 above all
# carry a {run} wildcard alongside {sample}) rather than concatenating raw
# reads across runs up front -- combine_sam_runs above merges the
# resulting per-run SAMs back into one per-sample SAM afterward, so every
# rule downstream of alignment is unaffected and still sees a single
# {sample}_trimmed_mapped.sam.
# ---------------------------------------------------------------------------

def get_runs_for_sample(wildcards):
    """Return every sequencing run identifier for a sample, sorted."""
    rows = samples_runs.loc[samples_runs["sample"] == wildcards.sample]
    if rows.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in samplesheet")
    return sorted(rows["run"].unique().tolist())


def get_fastq_for_sample_run(wildcards):
    """Return the single raw fastq for one (sample, run) pair."""
    rows = samples_runs.loc[
        (samples_runs["sample"] == wildcards.sample) & (samples_runs["run"] == wildcards.run)
    ]
    if rows.empty:
        raise ValueError(f"No row for sample '{wildcards.sample}', run '{wildcards.run}' in samplesheet")
    return rows["r1"].iloc[0]


def get_samples_by_condition(condition):
    rows = samples.loc[samples["condition"] == condition]
    if rows.empty:
        raise ValueError(f"No samples found for condition '{condition}'")
    return rows["sample"].tolist()


def get_samples_by_condition_and_temperature(condition, temperature):
    rows = samples.loc[(samples["condition"] == condition) & (samples["temperature"] == temperature)]
    if rows.empty:
        raise ValueError(f"No samples found for condition '{condition}' and temperature '{temperature}'")
    return rows["sample"].tolist()


rule rename_fastq:
    """
    Copy one sequencing run's raw FASTQ for a sample into the canonical
    per-(sample, run) location for downstream rules.
    """
    input:
        get_fastq_for_sample_run
    output:
        f"{TMP}/reads/renamed/{{sample}}_{{run}}.fastq"
    log:
        "logs/rename_fastq/{sample}_{run}.log"
    message:
        "Renaming FASTQ for sample {wildcards.sample}, run {wildcards.run}"
    shell:
        """
        cat {input} > {output} 2> {log}
        """


rule trim_reads:
    """
    Trim reads using the trimming script from StructureFold2. Runs per
    (sample, run) -- see the section comment above.
    """
    input:
        reads=f"{TMP}/reads/renamed/{{sample}}_{{run}}.fastq"
    output:
        fastq=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"
    threads: 4
    conda:
        "../envs/structurefold.yaml"
    params:
        tmpdir=f"{TMP}/structurefold2/{{sample}}_{{run}}",
        workdir=f"{workflow.basedir}/workflow"
    log:
        "logs/trim_reads/{sample}_{run}.log"
    message:
        "Trimming reads for sample {wildcards.sample}, run {wildcards.run}"
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
        mv {params.tmpdir}/{wildcards.sample}_{wildcards.run}_trimmed.fastq {output.fastq}
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
    log:
        "logs/filter_sam/{sample}.log"
    shell:
        """
        mkdir -p $(dirname {log}) && \
        LOG_ABS=$(readlink -f {log}) && \
        cd {params.sample_dir} &&
        python {params.workdir}/{params.script} -max_mismatch 3 -logname {params.log} > "$LOG_ABS" 2>&1
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
    log:
        "logs/add_header_to_filtered_sam/{sample}.log"
    shell:
        """
        (
        samtools view -H {input.original} > {output}
        grep -v '^@' {input.filtered} >> {output}
        ) > {log} 2>&1
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
    log:
        "logs/filtered_sam_to_sorted_bam/{sample}.log"
    shell:
        """
        (
        samtools view -bS {input} | samtools sort -o {output.bam}
        samtools index {output.bam}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

rule fastqc:
    """
    Run FastQC on trimmed reads for each sequencing run (trimming now
    happens per (sample, run) -- see the "Read preparation" comment above).
    """
    input:
        f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"
    output:
        html=f"{config['output_dir']}/qc/fastqc/{{sample}}_{{run}}_trimmed_fastqc.html",
        zip=f"{config['output_dir']}/qc/fastqc/{{sample}}_{{run}}_trimmed_fastqc.zip"
    params:
        outdir=f"{config['output_dir']}/qc/fastqc"
    log:
        "logs/fastqc/{sample}_{run}.log"
    message:
        "Running FastQC on trimmed reads for sample {wildcards.sample}, run {wildcards.run}"
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p {params.outdir} && fastqc {input} --outdir {params.outdir} > {log} 2>&1"


rule multiqc:
    """
    Aggregate FastQC reports across every (sample, run) with MultiQC
    """
    input:
        expand(
            f"{config['output_dir']}/qc/fastqc/{{sample}}_{{run}}_trimmed_fastqc.zip",
            zip,
            sample=samples_runs["sample"].tolist(),
            run=samples_runs["run"].tolist(),
        )
    output:
        report(
            f"{config['output_dir']}/qc/multiqc/multiqc_report.html",
            category="QC",
            subcategory="Alignment",
            caption="../report_captions/multiqc.rst",
        )
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
        files=expand(
            f"{config['output_dir']}/qc/alignment_stats/{{sample}}_stats.tsv",
            sample=SAMPLES
        )
    output:
        report(
            f"{config['output_dir']}/qc/alignment_stats_summary.tsv",
            category="QC",
            subcategory="Alignment",
            caption="../report_captions/alignment_stats_summary.rst",
        )
    log:
        "logs/aggregate_alignment_stats.log"
    message:
        "Aggregating alignment stats across all samples"
    shell:
        """
        (
        head -1 {input.files[0]} > {output}
        for f in {input.files}; do tail -n +2 "$f" >> {output}; done
        ) > {log} 2>&1
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
            files=expand(
                f"{config['output_dir']}/qc/positive_control_alignment/{{sample}}_pct.tsv",
                sample=SAMPLES
            )
        output:
            report(
                f"{config['output_dir']}/qc/positive_control_alignment_summary.tsv",
                category="QC",
                subcategory="Positive control",
                caption="../report_captions/positive_control_alignment_summary.rst",
            )
        log:
            "logs/aggregate_positive_control_alignment.log"
        message:
            "Aggregating positive-control alignment percentages across all samples"
        shell:
            """
            (
            head -1 {input.files[0]} > {output}
            for f in {input.files}; do tail -n +2 "$f" >> {output}; done
            ) > {log} 2>&1
            """
