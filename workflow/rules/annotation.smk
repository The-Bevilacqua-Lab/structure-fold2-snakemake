############################################################################
# Transcript position annotation rules
#
# Produces a per-position CSV for the active TRANSCRIPTOME labelling each
# base as 5UTR, CDS, 3UTR, or noncoding, with codon position (1/2/3) for
# CDS bases.
#
# Requires the following config keys to be set:
#   genome         – path to the reference genome FASTA
#   annotation_gtf – path to a GTF/GFF annotation for that genome
#
# If either key is absent the rules are not defined and the output is not
# added to rule all (see Snakefile).
############################################################################

if config.get("annotation_gtf") and config.get("genome"):

    rule extract_cds_and_ref_transcriptome:
        """
        Use gffread to extract CDS sequences (-x) and full reference
        transcript sequences (-w) from the genome + annotation.
        These reference sequences are used to locate CDS positions within
        each transcript regardless of which transcriptome FASTA the pipeline
        is using (e.g. one carrying substituted variants vs. the reference).
        """
        input:
            genome=config["genome"],
            annotation=config["annotation_gtf"],
        output:
            cds=f"{TMP}/resources/cds.fa",
            ref_tx=f"{TMP}/resources/ref_transcriptome.fa",
        log:
            "logs/extract_cds_and_ref_transcriptome.log",
        conda:
            "../envs/gffread.yaml"
        message:
            "Extracting CDS and reference transcript sequences with gffread"
        shell:
            """
            mkdir -p $(dirname {output.cds})
            gffread -g {input.genome} {input.annotation} \
                -x {output.cds} \
                -w {output.ref_tx} \
                >{log} 2>&1
            """

    rule annotate_transcript_positions:
        """
        Annotate each position in each transcript as 5UTR, CDS, 3UTR, or
        noncoding (transcripts with no CDS match). For CDS positions,
        report the codon position (1, 2, or 3).

        The CDS location in transcript space is determined by finding the
        gffread CDS sequence as a substring of the gffread reference
        transcript, then applied to the actual TRANSCRIPTOME (which may
        carry substituted variants but has the same exon structure).
        """
        input:
            ref_tx=f"{TMP}/resources/ref_transcriptome.fa",
            cds=f"{TMP}/resources/cds.fa",
            tx=TRANSCRIPTOME,
        output:
            f"{config['output_dir']}/transcript_position_annotations.csv",
        log:
            "logs/annotate_transcript_positions.log",
        conda:
            "../envs/biopython.yaml"
        message:
            "Annotating per-position UTR/CDS regions for all transcripts"
        shell:
            """
            python3 workflow/scripts/annotate_transcript_positions.py \
                --ref_tx_fasta {input.ref_tx} \
                --cds_fasta {input.cds} \
                --tx_fasta {input.tx} \
                --output {output} \
                2>{log}
            """


# ---------------------------------------------------------------------------
# Optional: reactivity_region (see the Snakefile comment near
# REACTIVITY_REGION) restricts the main reactivity calculation to one mRNA
# region. These two rules build the region-sliced coordinates and
# transcriptome that region_rtsc (workflow/rules/structurefold2.smk) and
# rtsc_to_react then use in place of the full transcript. REACTIVITY_REGION
# being truthy already guarantees 'genome'/'annotation_gtf' are set (the
# Snakefile raises otherwise), so transcript_position_annotations above is
# always defined when this block is.
# ---------------------------------------------------------------------------
if REACTIVITY_REGION:

    rule region_coordinates:
        """
        Locate each transcript's reactivity_region boundaries (start/end, in
        the same transcript coordinates as TRANSCRIPTOME) from the
        per-position annotation CSV.
        """
        input:
            f"{config['output_dir']}/transcript_position_annotations.csv",
        output:
            f"{config['output_dir']}/region_coordinates.tsv",
        log:
            "logs/region_coordinates.log",
        conda:
            "../envs/biopython.yaml"
        params:
            region=REACTIVITY_REGION,
        message:
            "Locating reactivity_region boundaries for the region-restricted reactivity calculation"
        shell:
            """
            python3 workflow/scripts/region_coordinates.py \
                --annotations {input} --region {params.region} --output {output} \
                2>{log}
            """

    rule region_transcriptome:
        """
        Slice TRANSCRIPTOME_UPPER down to just reactivity_region per
        transcript. Used by rtsc_to_react/convert_react_to_csv instead of
        the full transcriptome when reactivity_region is set.
        """
        input:
            fasta=TRANSCRIPTOME_UPPER,
            coords=f"{config['output_dir']}/region_coordinates.tsv",
        output:
            f"{TMP}/resources/transcriptome_region.fa",
        log:
            "logs/region_transcriptome.log",
        conda:
            "../envs/biopython.yaml"
        shell:
            """
            python3 workflow/scripts/extract_region_fasta.py \
                --fasta {input.fasta} --coords {input.coords} --output {output} \
                2>{log}
            """
