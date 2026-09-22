############################################################################
# SAM/BAM post-processing
############################################################################

rule filter_sam:
    """
    Filter mapped reads by mismatch count using the Python 3 port of
    StructureFold2's sam_filter (scripts/StructureFold3/sam_filter.py).
    """
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"
    output:
        sam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
    threads: 4
    conda:
        "../envs/samtools.yaml"
    params:
        script="workflow/scripts/StructureFold3/sam_filter.py",
        filter_log=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.log",
        max_mismatch=3,
    log:
        "logs/filter_sam/{sample}.log"
    message:
        "Filtering mapped reads for sample {wildcards.sample}"
    shell:
        "python3 {params.script} -sam {input} -max_mismatch {params.max_mismatch} "
        "-logname {params.filter_log} -threads {threads} > {log} 2>&1"


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
