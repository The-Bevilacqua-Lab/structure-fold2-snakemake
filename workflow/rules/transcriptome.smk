############################################################################
# Transcriptome preparation
############################################################################


if TRANSCRIPT_SOURCE == "gffread" and not GENOTYPE_MODE:   # genotype mode: genotypes.smk
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

    if CDS_ONLY:

        rule filter_extracted_seqs:
            """
            Keep only those transcripts with a CDS (option cds_only).
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
                "python3 workflow/scripts/transcriptome/filter_extracted_seqs.py --input {input} --output {output} > {log} 2>&1"


_EXTRA_RNA_FROM_ANNOTATION = [k for k, v in EXTRA_RNA.items() if "annotation" in v]
if _EXTRA_RNA_FROM_ANNOTATION:

    rule gffread_extract_extra_rna:
        """
        Extract ncRNA / lncRNA sequences from a genome with gffread, when they are given as a
        GFF3/GTF (ncrna_gff/gtf, lncrna_gff/gtf) instead of a FASTA. With <kind>_noncoding_only
        (default true) transcripts with a CDS are discarded (gffread --nc).
        """
        input:
            genome=lambda wc: EXTRA_RNA[wc.kind]["genome"],
            annotation=lambda wc: EXTRA_RNA[wc.kind]["annotation"],
        output:
            f"{TMP}/resources/{{kind}}_gffread.fa"
        wildcard_constraints:
            kind="|".join(_EXTRA_RNA_FROM_ANNOTATION)
        params:
            nc=lambda wc: "--nc" if EXTRA_RNA[wc.kind]["noncoding_only"] else "",
        conda:
            "../envs/gffread.yaml"
        log:
            "logs/gffread/{kind}.log"
        message:
            "Extracting {wildcards.kind} sequences with gffread"
        shell:
            "gffread {input.annotation} -g {input.genome} {params.nc} -w {output} > {log} 2>&1"


rule prepare_transcriptome:
    """
    Assemble the canonical transcriptome FASTA at resources/transcriptome.fa.

    Base sequences come from _source_transcriptome (either the gffread
    pipeline output or the user-supplied FASTA, depending on
    transcript_source). Any extras from _extra_transcriptome_fastas
    (positive controls, ncRNAs, lncRNAs) are appended after the base sequences;
    a record whose ID already appeared earlier is dropped (logged).
    """
    input:
        transcriptome=_source_transcriptome,
        extras=_extra_transcriptome_fastas
    output:
        TRANSCRIPTOME
    params:
        awk="workflow/scripts/transcriptome/concat_fasta_unique.awk"
    log:
        "logs/prepare_transcriptome/prepare.log"
    message:
        "Preparing final transcriptome FASTA"
    shell:
        "awk -f {params.awk} {input.transcriptome} {input.extras} > {output} 2> {log}"


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
        "python3 workflow/scripts/transcriptome/uppercase_fasta.py --input {input} --output {output} > {log} 2>&1"
