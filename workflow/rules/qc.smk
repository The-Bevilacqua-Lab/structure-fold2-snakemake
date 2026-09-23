############################################################################
# Quality control
############################################################################

rule fastqc:
    """
    Run FastQC on trimmed reads for each sequencing run (trimming now
    happens per (sample, run) -- see the comment at the top of read_prep.smk).
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


if TRIMMER == "fastp":

    rule multiqc_fastp:
        """
        Aggregate fastp's own trimming stats (reads before/after filtering,
        adapter content, quality -- MultiQC reads these straight out of
        fastp's <.json>, not the <.html>) across every (sample, run) with
        MultiQC. Only defined when config['trimmer'] is "fastp" -- the
        cutadapt trimmer has no equivalent per-run json report.
        """
        input:
            expand(
                f"{TMP}/reads/trimmed/{{sample}}_{{run}}_fastp.json",
                zip,
                sample=samples_runs["sample"].tolist(),
                run=samples_runs["run"].tolist(),
            )
        output:
            report(
                f"{config['output_dir']}/qc/multiqc_fastp/multiqc_report.html",
                category="QC",
                subcategory="Trimming",
                caption="../report_captions/multiqc_fastp.rst",
            )
        params:
            outdir=f"{config['output_dir']}/qc/multiqc_fastp"
        log:
            "logs/multiqc_fastp.log"
        message:
            "Running MultiQC on fastp trimming reports"
        conda:
            "../envs/qc.yaml"
        shell:
            "multiqc --force {input} --outdir {params.outdir} > {log} 2>&1"


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
            "../envs/samtools.yaml"
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
