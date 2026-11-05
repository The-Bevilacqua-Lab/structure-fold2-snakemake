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
# reactivity calculation
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

# Optional depth-matched reactivity (workflow/rules/depth_matching.smk) -- see
# config.yaml's depth_matching comment.
DEPTH_MATCHING = bool(config.get("depth_matching", False))
DEPTH_MATCHING_DRAWS = int(config.get("depth_matching_draws", 100))
DEPTH_MATCHING_MATCH_MINUS = bool(config.get("depth_matching_match_minus", True))
DEPTH_MATCHING_ROC_TRIM3 = int(config.get("depth_matching_roc_trim3", 50))

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule (transcriptome.smk).
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"

# Uppercase-only copy of TRANSCRIPTOME, produced by uppercase_transcriptome
# (transcriptome.smk).
TRANSCRIPTOME_UPPER = f"{TMP}/resources/transcriptome_upper.fa"

# ---------------------------------------------------------------------------
# Optional: estimate per-transcript relative abundance (RPKM and/or TPM)
# from the untreated (-DMS) RT-stop counts, via StructureFold2's own
# rtsc_abundances.py 
# ---------------------------------------------------------------------------
_ALLOWED_ABUNDANCE_MODES = {"RPKM", "TPM"}
ABUNDANCE_MODES = config.get("transcript_abundance", [])
if isinstance(ABUNDANCE_MODES, str):
    ABUNDANCE_MODES = [ABUNDANCE_MODES]
ABUNDANCE_MODES = sorted({m.upper() for m in ABUNDANCE_MODES})
_invalid_abundance_modes = set(ABUNDANCE_MODES) - _ALLOWED_ABUNDANCE_MODES
if _invalid_abundance_modes:
    raise ValueError(
        f"config['transcript_abundance'] must only contain {sorted(_ALLOWED_ABUNDANCE_MODES)}, "
        f"got {sorted(_invalid_abundance_modes)}"
    )


# The p4p6 structure plot rules depend on a hand-annotated coordinate map
# specific to that one construct (see workflow/rules/positive_control.smk) --
# only wired in when the positive control actually is p4p6. Any other
# positive control still gets the generic alignment-% QC from qc.smk.
POSITIVE_CONTROL_IS_P4P6 = config.get("positive_control_name") == "p4p6"
