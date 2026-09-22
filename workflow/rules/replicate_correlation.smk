############################################################################
# Replicate correlation of RT-stop counts
############################################################################

rule rtsc_stop_correlation:
    """RT-stop correlation across the true +DMS replicates (unpooled)."""
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_plus.rtsc",
            id=CORRELATION_PLUS_IDS,
        ),
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv",
    log:
        "logs/rtsc_stop_correlation.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_correlation.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 {params.workdir}/{params.script} \
            {input.rtsc} \
            -spec AC \
            -fasta {params.transcriptome} \
            -name {output} >{log} 2>&1
        """


rule plot_rt_stop_correlation:
    """Pairwise replicate correlation plot (+DMS)."""
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation.csv",
    output:
        report(
            f"{config['output_dir']}/qc/rt_stop_replicate_correlation.png",
            category="QC",
            subcategory="Replicate correlation",
            caption="../report_captions/rt_stop_replicate_correlation.rst",
        ),
    log:
        "logs/plot_rt_stop_correlation.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output} >{log} 2>&1
        """


rule rtsc_stop_correlation_minus:
    """RT-stop correlation across the true -DMS replicates (unpooled)."""
    input:
        rtsc=expand(
            f"{TMP}/output/se/{{id}}/combined_minus.rtsc",
            id=CORRELATION_MINUS_IDS,
        ),
    output:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv",
    log:
        "logs/rtsc_stop_correlation_minus.log",
    conda:
        "../envs/rtsc_tools.yaml"
    params:
        workdir=f"{workflow.basedir}/workflow",
        script="scripts/StructureFold3/rtsc_correlation.py",
        transcriptome=TRANSCRIPTOME_UPPER,
    shell:
        """
        python3 {params.workdir}/{params.script} \
            {input.rtsc} \
            -spec AC \
            -fasta {params.transcriptome} \
            -name {output} >{log} 2>&1
        """


rule plot_rt_stop_correlation_minus:
    """Pairwise replicate correlation plot (-DMS)."""
    input:
        f"{config['output_dir']}/qc/rt_stop_correlation_minus.csv",
    output:
        report(
            f"{config['output_dir']}/qc/rt_stop_replicate_correlation_minus.png",
            category="QC",
            subcategory="Replicate correlation",
            caption="../report_captions/rt_stop_replicate_correlation_minus.rst",
        ),
    log:
        "logs/plot_rt_stop_correlation_minus.log",
    conda:
        "../envs/plotting.yaml"
    shell:
        """
        python3 workflow/scripts/plot_rt_stop_correlation.py \
            --input {input} \
            --output {output} >{log} 2>&1
        """
