############################################################################
# Coverage QC
############################################################################

rule samtools_depth_by_condition:
    """samtools depth (all positions) for one condition."""
    input:
        lambda wc: get_sams_by_condition_and_id(wc.condition, wc.id),
    output:
        f"{config['output_dir']}/{{id}}/depth_{{condition}}.txt",
    log:
        "logs/samtools_depth_by_condition/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
        condition="plus|minus",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && samtools depth -a {input} >{output} 2>{log}
        """


rule calculate_stop_coverage:
    """Per-transcript RT-stop coverage of the +DMS counts."""
    input:
        f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
    output:
        coverage=f"{config['output_dir']}/{{id}}/coverage.csv",
    log:
        "logs/calculate_stop_coverage/{id}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        coverage_name=f"{config['output_dir']}/{{id}}/coverage",
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_coverage.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        mkdir -p $(dirname {output.coverage}) \
            && python3 {params.workdir}/{params.script} -f {input} -name {params.coverage_name} {params.transcriptome} >{log} 2>&1
        """


rule calculate_nucleotide_coverage:
    """Per-nucleotide RT-stop coverage."""
    input:
        rtsc=f"{TMP}/output/se/{{id}}/combined_{{condition}}.rtsc",
    output:
        f"{config['output_dir']}/{{id}}/nucleotide_coverage_{{condition}}.csv",
    log:
        "logs/calculate_nucleotide_coverage/{id}_{condition}.log",
    wildcard_constraints:
        id="[^/]+",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 workflow/scripts/rtsc_to_nucleotide_coverage.py \
            --rtsc {input.rtsc} \
            --fasta {params.transcriptome} \
            --output {output} >{log} 2>&1
        """


rule filter_covered_transcripts:
    """List transcripts with mean RT-stop coverage >= 1 in every replicate."""
    input:
        expand(
            "{out}/{id}/coverage.csv",
            out=config["output_dir"],
            id=IDS,
        ),
    output:
        f"{config['output_dir']}/qc/covered_transcripts.txt",
    log:
        "logs/filter_covered_transcripts.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/filter_covered_transcripts.py \
                --coverages {input} \
                --threshold 1 \
                --output {output} >{log} 2>&1
        """


rule upset_covered_transcripts_replicates:
    """UpSet plot of covered transcripts across the true (unpooled) replicates."""
    input:
        expand(f"{config['output_dir']}/{{id}}/coverage.csv", id=CORRELATION_PLUS_IDS),
    output:
        f"{config['output_dir']}/qc/covered_transcripts_upset_replicates.png",
    log:
        "logs/upset_covered_transcripts_replicates.log",
    conda:
        "../envs/upset.yaml"
    params:
        labels=CORRELATION_PLUS_IDS,
        threshold=1.0,
    shell:
        """
        mkdir -p $(dirname {output}) \
            && python3 workflow/scripts/plot_covered_transcripts_upset.py \
                --coverages {input} \
                --labels {params.labels} \
                --threshold {params.threshold} \
                --output {output} >{log} 2>&1
        """


# Only when pooling leaves more than one set to compare (not "both").
if MERGED_IDS_DIFFER and len(IDS) > 1:

    rule upset_covered_transcripts_merged:
        """UpSet plot of covered transcripts across the pooled grouping (IDS)."""
        input:
            expand(f"{config['output_dir']}/{{id}}/coverage.csv", id=IDS),
        output:
            f"{config['output_dir']}/qc/covered_transcripts_upset_merged.png",
        log:
            "logs/upset_covered_transcripts_merged.log",
        conda:
            "../envs/upset.yaml"
        params:
            labels=IDS,
            threshold=1.0,
        shell:
            """
            mkdir -p $(dirname {output}) \
                && python3 workflow/scripts/plot_covered_transcripts_upset.py \
                    --coverages {input} \
                    --labels {params.labels} \
                    --threshold {params.threshold} \
                    --output {output} >{log} 2>&1
            """
