############################################################################
# Metagene coverage plots (one line per sample, all samples on one plot)
#
# Per transcript, coverage c becomes ln(c + 1) / mean(ln(c + 1)) before being
# binned by relative position and averaged across transcripts -- see
# workflow/scripts/metagene_profile.py. Two plots, one per coverage source:
#   depth -- samtools depth (-a) on each sample's filtered, sorted BAM
#   rtsc  -- each sample's RT-stop counts (.rtsc)
############################################################################

METAGENE_BINS = config.get("metagene_bins", 100)
METAGENE_DIR = f"{config['output_dir']}/qc/metagene"


rule metagene_profile_depth:
    """
    Metagene profile of one sample's samtools-depth coverage. The depth
    stream is piped straight into the script, so the (transcriptome-sized)
    per-position depth table is never written to disk.
    """
    input:
        bam=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered_sorted.bam",
    output:
        f"{METAGENE_DIR}/profiles/{{sample}}_depth.tsv",
    params:
        bins=METAGENE_BINS,
    log:
        "logs/metagene_profile_depth/{sample}.log",
    conda:
        "../envs/metagene.yaml"
    message:
        "Metagene depth profile for sample {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output})
        samtools depth -a {input.bam} 2>{log} \
            | python3 workflow/scripts/metagene_profile.py --mode depth --bins {params.bins} --output {output} 2>>{log}
        """


rule metagene_profile_rtsc:
    """Metagene profile of one sample's RT-stop counts."""
    input:
        rtsc=f"{TMP}/output/se/{{sample}}/{{sample}}_trimmed_mapped_filtered.rtsc",
    output:
        f"{METAGENE_DIR}/profiles/{{sample}}_rtsc.tsv",
    params:
        bins=METAGENE_BINS,
    log:
        "logs/metagene_profile_rtsc/{sample}.log",
    conda:
        "../envs/metagene.yaml"
    message:
        "Metagene RT-stop profile for sample {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output})
        python3 workflow/scripts/metagene_profile.py --mode rtsc --input {input.rtsc} \
            --bins {params.bins} --output {output} >{log} 2>&1
        """


rule plot_metagene_depth:
    """Overlay every sample's samtools-depth metagene profile on one plot."""
    input:
        expand(f"{METAGENE_DIR}/profiles/{{sample}}_depth.tsv", sample=SAMPLES),
    output:
        report(
            f"{METAGENE_DIR}/metagene_depth.png",
            category="QC",
            subcategory="Metagene coverage",
            caption="../report_captions/metagene_depth.rst",
        ),
    params:
        samples=SAMPLES,
        conditions=samples["condition"].tolist(),
    log:
        "logs/plot_metagene_depth.log",
    conda:
        "../envs/metagene.yaml"
    shell:
        """
        python3 workflow/scripts/plot_metagene.py \
            --profiles {input} --samples {params.samples} --conditions {params.conditions} \
            --title "Metagene coverage (samtools depth)" --output {output} >{log} 2>&1
        """


rule plot_metagene_rtsc:
    """Overlay every sample's RT-stop metagene profile on one plot."""
    input:
        expand(f"{METAGENE_DIR}/profiles/{{sample}}_rtsc.tsv", sample=SAMPLES),
    output:
        report(
            f"{METAGENE_DIR}/metagene_rtsc.png",
            category="QC",
            subcategory="Metagene coverage",
            caption="../report_captions/metagene_rtsc.rst",
        ),
    params:
        samples=SAMPLES,
        conditions=samples["condition"].tolist(),
    log:
        "logs/plot_metagene_rtsc.log",
    conda:
        "../envs/metagene.yaml"
    shell:
        """
        python3 workflow/scripts/plot_metagene.py \
            --profiles {input} --samples {params.samples} --conditions {params.conditions} \
            --title "Metagene coverage (RT-stops)" --output {output} >{log} 2>&1
        """
