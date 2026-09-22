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

# Bases excluded from the 3' end of each transcript when generating the
# 2-8% normalization scale (see config.yaml's trim3 comment).
TRIM3 = config.get("trim3", 0)

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule (transcriptome.smk).
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"

# Uppercase-only copy of TRANSCRIPTOME, produced by uppercase_transcriptome
# (transcriptome.smk). Used (instead of TRANSCRIPTOME) by every rule downstream of
# mapping that checks base identity -- see that rule's docstring for why.
TRANSCRIPTOME_UPPER = f"{TMP}/resources/transcriptome_upper.fa"
