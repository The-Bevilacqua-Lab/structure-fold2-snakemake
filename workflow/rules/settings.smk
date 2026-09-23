############################################################################
# Shared settings: paths and options read from config.yaml
############################################################################

import os

TRANSCRIPT_SOURCE = config.get("transcript_source", "file")

# gffread source only: run AGAT first to keep just the longest isoform per gene
# (off by default -- gffread then reads the GFF as given).
KEEP_LONGEST_ISOFORM = bool(config.get("keep_longest_isoform", False))
TMP = config.get("tmp_dir", "tmp")

# Read-trimming tool: "cutadapt" (default, StructureFold2's fastq_trimmer.py)
# or "fastp" -- see config.yaml's trimmer comment.
TRIMMER = config.get("trimmer", "cutadapt")
if TRIMMER not in ("cutadapt", "fastp"):
    raise ValueError(
        f"config['trimmer'] must be 'cutadapt' or 'fastp', got {TRIMMER!r}"
    )

# Optional 3' coverage-bias correction of the +DMS RT-stops before the
# reactivity calculation -- see config.yaml's 3_prime_bias_correction comment
# and workflow/rules/bias_correction.smk.
BIAS_CORRECTION = bool(config.get("3_prime_bias_correction", False))
BIAS_REFERENCE = config.get("3_prime_bias_reference", "pooled")
BIAS_N_TRANSCRIPTS = int(config.get("3_prime_bias_n_transcripts", 2500))
BIAS_BINS = int(config.get("3_prime_bias_bins", 50))
BIAS_ANNOTATION = config.get(
    "3_prime_bias_annotation", config.get("annotation_gtf", config.get("gff"))
)
BIAS_TMP = f"{TMP}/bias_correction"
if BIAS_CORRECTION:
    if BIAS_REFERENCE not in ("pooled", "independent"):
        raise ValueError(
            "config['3_prime_bias_reference'] must be 'pooled' or 'independent', "
            f"got {BIAS_REFERENCE!r}"
        )
    if not BIAS_ANNOTATION:
        raise ValueError(
            "3_prime_bias_correction needs a GTF/GFF3 to find protein-coding "
            "transcripts: set 3_prime_bias_annotation (or annotation_gtf/gff)."
        )

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule (transcriptome.smk).
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"

# Uppercase-only copy of TRANSCRIPTOME, produced by uppercase_transcriptome
# (transcriptome.smk). Used (instead of TRANSCRIPTOME) by every rule downstream of
# mapping that checks base identity -- see that rule's docstring for why.
TRANSCRIPTOME_UPPER = f"{TMP}/resources/transcriptome_upper.fa"
