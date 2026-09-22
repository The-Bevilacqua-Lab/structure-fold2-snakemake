############################################################################
# Read preparation: rename -> trim
#
# Each sequencing run of a sample is renamed, trimmed, and aligned
# separately (rename_fastq / trim_reads / align_with_bowtie2 all carry a
# {run} wildcard alongside {sample}) rather than concatenating raw reads
# across runs up front -- combine_sam_runs (alignment.smk) merges the
# resulting per-run SAMs back into one per-sample SAM afterward, so every
# rule downstream of alignment is unaffected and still sees a single
# {sample}_trimmed_mapped.sam.
############################################################################

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

if TRIMMER == "fastp":
    rule trim_reads_fastp:
        """
        Trim reads using fastp. Runs per (sample, run) -- see the section
        comment at the top of this file.
        """
        input:
            sample=[f"{TMP}/reads/renamed/{{sample}}_{{run}}.fastq"]
        output:
            trimmed=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq",
            failed=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.failed.fastq",
            html=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_fastp.html",
            json=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_fastp.json"
        log:
            "logs/trim_reads_fastp/{sample}_{run}.log"
        params:
            adapters="--adapter_sequence GATCGGAAGAGCACACGTCTG",
            extra="--trim_poly_g"
        threads: 4
        wrapper:
            "v3.3.3/bio/fastp"

if TRIMMER == "cutadapt":
    rule trim_reads_cutadapt:
        """
        Trim reads with the Python 3 port of StructureFold2's fastq_trimmer.py. Runs per
        (sample, run) -- see the section comment at the top of this file.
        """
        input:
            reads=f"{TMP}/reads/renamed/{{sample}}_{{run}}.fastq"
        output:
            fastq=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"
        threads: 4
        conda:
            "../envs/cutadapt.yaml"
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
            python3 {params.workdir}/scripts/StructureFold3/fastq_trimmer.py
            cd {params.workdir}
            cd ..
            mkdir -p $(dirname {output.fastq})
            mv {params.tmpdir}/{wildcards.sample}_{wildcards.run}_trimmed.fastq {output.fastq}
            """
