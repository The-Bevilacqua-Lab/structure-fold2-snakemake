############################################################################
# SAM filtering and the sorted BAM
############################################################################

rule filter_sam:
    """
    StructureFold2's sam_filter (Python 3 port): keep forward-strand alignments with at most
    max_mismatch mismatches/indels and no mismatch at the read's first base.
    """
    input:
        f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam"
    output:
        sam=intermediate(f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    params:
        script="workflow/scripts/StructureFold3/sam_filter.py",
        filter_log=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.log",
        max_mismatch=MAX_MISMATCH,
    log:
        "logs/filter_sam/{sample}.log"
    message:
        "Filtering mapped reads for sample {wildcards.sample}"
    shell:
        "python3 {params.script} -sam {input} -max_mismatch {params.max_mismatch} "
        "-logname {params.filter_log} -threads {threads} > {log} 2>&1"


rule filtered_sam_to_sorted_bam:
    """
    Sorted, indexed BAM of the filtered alignments (used for read depth: metagene and 3' bias
    correction). sam_filter.py drops the SAM header, so the original header is put back first.
    """
    input:
        filtered=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.sam",
        original=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped.sam",
    output:
        bam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_sorted.bam",
        bai=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_sorted.bam.bai",
    threads: 4
    resources:
        mem_mb=8000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/filtered_sam_to_sorted_bam/{sample}.log"
    shell:
        """
        exec 2> {log}
        {{ samtools view -H {input.original}; grep -v '^@' {input.filtered} || true; }} \
            | samtools sort -@ {threads} -m 1500M -T {output.bam}.tmp -o {output.bam} -
        samtools index {output.bam}
        """
