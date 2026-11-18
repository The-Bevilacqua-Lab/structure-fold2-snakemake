############################################################################
# Read trimming
#
# Each sequencing run of a sample is trimmed and aligned on its own ({run} wildcard);
# combine_sam_from_multiple_runs (alignment.smk) then merges the runs into one SAM per sample.
# The trimmers read the samplesheet's FASTQ (plain or gzipped) directly.
############################################################################

if TRIMMER == "fastp":

    rule trim_reads_fastp:
        """Trim the 3' adapter, poly-G tails and short reads with fastp."""
        input:
            sample=lambda wc: [get_fastq_for_sample_run(wc)],   # the wrapper wants a list
        output:
            trimmed=intermediate(f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"),
            failed=temp(f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.failed.fastq"),
            html=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_fastp.html",
            json=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_fastp.json",
        log:
            "logs/trim_reads_fastp/{sample}_{run}.log",
        params:
            adapters=f"--adapter_sequence {FASTP_ADAPTER}",
            extra=f"--trim_poly_g --length_required {FASTP_MIN_LENGTH}",
        threads: 4
        wrapper:
            "v3.3.3/bio/fastp"

else:

    rule trim_reads_cutadapt:
        """
        Trim with StructureFold2's fastq_trimmer.py (Python 3 port): cutadapt removes the 5' and
        3' Structure-seq2 adapters, then bases below Q30 and reads shorter than 20 nt.
        """
        input:
            get_fastq_for_sample_run,
        output:
            intermediate(f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"),
        log:
            "logs/trim_reads/{sample}_{run}.log",
        threads: 3   # the three cutadapt steps run at the same time, connected by pipes
        conda:
            "../envs/cutadapt.yaml"
        message:
            "Trimming reads for sample {wildcards.sample}, run {wildcards.run}"
        shell:
            "python3 workflow/scripts/StructureFold3/fastq_trimmer.py -input {input} -output {output} > {log} 2>&1"
