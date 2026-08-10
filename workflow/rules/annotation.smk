############################################################################
# Transcript position annotation rules
#
# Produces a per-position TSV for the active TRANSCRIPTOME labelling each
# base as 5UTR, CDS, or 3UTR, with codon position (1/2/3) for CDS bases.
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
            genome     = config["genome"],
            annotation = config["annotation_gtf"],
        output:
            cds    = f"{TMP}/resources/cds.fa",
            ref_tx = f"{TMP}/resources/ref_transcriptome.fa",
        conda:
            "../envs/gffread.yaml"
        log:
            "logs/extract_cds_and_ref_transcriptome.log"
        message:
            "Extracting CDS and reference transcript sequences with gffread"
        shell:
            """
            mkdir -p $(dirname {output.cds})
            gffread -g {input.genome} {input.annotation} \
                -x {output.cds} \
                -w {output.ref_tx} \
                > {log} 2>&1
            """


    rule annotate_transcript_positions:
        """
        Annotate each position in each transcript as 5UTR, CDS, or 3UTR.
        For CDS positions, report the codon position (1, 2, or 3).

        The CDS location in transcript space is determined by finding the
        gffread CDS sequence as a substring of the gffread reference
        transcript, then applied to the actual TRANSCRIPTOME (which may
        carry substituted variants but has the same exon structure).
        """
        input:
            ref_tx = f"{TMP}/resources/ref_transcriptome.fa",
            cds    = f"{TMP}/resources/cds.fa",
            tx     = TRANSCRIPTOME,
        output:
            f"{config['output_dir']}/transcript_position_annotations.tsv"
        conda:
            "../envs/biopython.yaml"
        log:
            "logs/annotate_transcript_positions.log"
        message:
            "Annotating per-position UTR/CDS regions for all transcripts"
        shell:
            """
            python3 workflow/scripts/annotate_transcript_positions.py \
                --ref_tx_fasta {input.ref_tx} \
                --cds_fasta    {input.cds} \
                --tx_fasta     {input.tx} \
                --output       {output} \
                2> {log}
            """
