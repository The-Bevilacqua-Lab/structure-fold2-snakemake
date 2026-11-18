############################################################################
# Alignment with bowtie2
############################################################################

rule build_bowtie2_index:
    """
    Build a Bowtie2 index from the canonical transcriptome FASTA.
    """
    input:
        TRANSCRIPTOME
    output:
        f"{TRANSCRIPTOME}.1.bt2"
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2_index/build.log"
    message:
        "Building Bowtie2 index for the transcriptome"
    shell:
        "bowtie2-build {input} {input} > {log} 2>&1"


rule align_with_bowtie2:
    """
    Align one sequencing run's trimmed reads to the reference transcriptome
    using Bowtie2.
    """
    input:
        transcriptome_index=get_bowtie2_index_for_sample,
        fastq1=f"{TMP}/reads/trimmed/{{sample}}_{{run}}_trimmed.fastq"
    output:
        sam=intermediate(f"{TMP}/output/se_runs/{{sample}}_{{run}}/{{sample}}_{{run}}_trimmed_mapped.sam"),
        mapping_log=f"{TMP}/output/se_runs/{{sample}}_{{run}}/{{sample}}_{{run}}_mapping_log.txt"
    params:
        # Bowtie2 index prefix: the .1.bt2 input minus its suffix
        index=lambda wc, input: input.transcriptome_index[: -len(".1.bt2")]
    threads: 4
    log:
        "logs/align_with_bowtie2/{sample}_{run}.log"
    conda:
        "../envs/bowtie2.yaml"
    message:
        "Aligning trimmed reads for sample {wildcards.sample}, run {wildcards.run} with Bowtie2"
    shell:
        # SAM goes to the output file; bowtie2's alignment summary (stderr) is
        # kept as the mapping_log output (read by the combine/stats rules)
        # and copied to the rule log.
        "bowtie2 -a -p {threads} -x {params.index} -q {input.fastq1} "
        "2>&1 >{output.sam} | tee {output.mapping_log} > {log}"


rule combine_sam_from_multiple_runs:
    """
    One SAM per sample: the runs' alignments one after another under the first run's header (a
    single run is hard-linked, not copied), and the runs' Bowtie2 summaries added up.
    """
    input:
        sams=get_run_sams_for_sample,
        logs=get_run_mapping_logs_for_sample,
    output:
        sam=intermediate(f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"),
        mapping_log=f"{TMP}/output/se/{{sample}}/{{sample}}_mapping_log.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/combine_sam_runs/{sample}.log"
    message:
        "Combining sequencing runs of sample {wildcards.sample} into one SAM"
    shell:
        r"""
        exec > {log} 2>&1
        set -- {input.sams}
        if [ $# -eq 1 ]; then
            ln -f "$1" {output.sam} || cp "$1" {output.sam}
        else
            samtools view -H "$1" > {output.sam}
            for f in "$@"; do grep -v '^@' "$f" >> {output.sam}; done
        fi
        python3 workflow/scripts/alignment/combine_mapping_logs.py {input.logs} -o {output.mapping_log}
        """
