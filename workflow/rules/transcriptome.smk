############################################################################
# Transcriptome preparation
############################################################################


if TRANSCRIPT_SOURCE == "gffread":
    # gffread reads the GFF as given, or the longest-isoform GFF from AGAT
    GFFREAD_GFF = f"{TMP}/resources/longest_isoform.gff" if KEEP_LONGEST_ISOFORM else config["gff"]

    if KEEP_LONGEST_ISOFORM:

        rule agat_longest_isoform:
            """
            Use AGAT to retain only the longest isoform per gene (option keep_longest_isoform).
            """
            input:
                gff=config["gff"]
            output:
                f"{TMP}/resources/longest_isoform.gff"
            conda:
                "../envs/agat.yaml"
            log:
                "logs/agat_longest_isoform/agat.log"
            message:
                "Filtering GFF to longest isoform per gene with AGAT"
            shell:
                "agat_sp_keep_longest_isoform.pl --gff {input.gff} --out {output} > {log} 2>&1"

    rule gffread_extract_transcripts:
        """
        Extract transcript sequences from the genome FASTA using gffread.
        """
        input:
            genome=config["genome"],
            gff=GFFREAD_GFF
        output:
            f"{TMP}/resources/transcriptome_gffread.fa"
        conda:
            "../envs/gffread.yaml"
        log:
            "logs/gffread/gffread.log"
        message:
            "Extracting transcript sequences with gffread"
        shell:
            "gffread {input.gff} -g {input.genome} -w {output} > {log} 2>&1"

    rule filter_extracted_seqs:
        """
        Keep only those transcripts with a CDS
        """
        input:
            f"{TMP}/resources/transcriptome_gffread.fa"
        output:
            f"{TMP}/resources/transcriptome_gffread_filtered.fa"
        conda:
            "../envs/biopython.yaml"
        log:
            "logs/filter_extracted_seqs/filter.log"
        shell:
            "python3 workflow/scripts/filter_extracted_seqs.py --input {input} --output {output} > {log} 2>&1"


rule prepare_transcriptome:
    """
    Assemble the canonical transcriptome FASTA at resources/transcriptome.fa.

    Base sequences come from _source_transcriptome (either the gffread
    pipeline output or the user-supplied FASTA, depending on
    transcript_source). Any extras from _extra_transcriptome_fastas
    (positive controls, ncRNAs) are appended after the base sequences.
    """
    input:
        transcriptome=_source_transcriptome,
        extras=_extra_transcriptome_fastas
    output:
        TRANSCRIPTOME
    log:
        "logs/prepare_transcriptome/prepare.log"
    message:
        "Preparing final transcriptome FASTA"
    shell:
        "awk 1 {input.transcriptome} {input.extras} > {output} 2> {log}"


rule uppercase_transcriptome:
    """
    Write an uppercase-only copy of TRANSCRIPTOME (TRANSCRIPTOME_UPPER).

    Bowtie2 alignment is case-insensitive, so a source transcriptome FASTA
    with soft-masked (lowercase) bases still maps correctly -- but every
    base-identity check downstream of mapping (reactivity/specificity "is
    this position an A or C" tests, UTR/CDS annotation) compares
    case-sensitively against an uppercase specificity set, so those
    lowercase bases would otherwise be silently skipped. Kept as a separate
    file (rather than uppercasing TRANSCRIPTOME itself) so the already-built
    Bowtie2 index/alignments never need to be redone for this.
    """
    input:
        TRANSCRIPTOME
    output:
        TRANSCRIPTOME_UPPER
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/uppercase_transcriptome/uppercase.log"
    message:
        "Uppercasing transcriptome for downstream base-identity checks"
    shell:
        "python3 workflow/scripts/uppercase_fasta.py --input {input} --output {output} > {log} 2>&1"
