############################################################################
# Shared settings: paths and options read from config.yaml
############################################################################

import os


def _flag(value, default=False):
    """A config on/off value: YAML booleans, or strings like "false"/"no"/"0"/"off" (from --config)."""
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() not in ("false", "no", "0", "off", "")
    return bool(value)


TRANSCRIPT_SOURCE = config.get("transcript_source", "file")
if TRANSCRIPT_SOURCE not in ("file", "gffread"):
    raise ValueError(f"config['transcript_source'] must be 'file' or 'gffread', got {TRANSCRIPT_SOURCE!r}")

# gffread source only: run AGAT first to keep just the longest isoform per gene
# (off by default -- gffread then reads the GFF as given).
KEEP_LONGEST_ISOFORM = _flag(config.get("keep_longest_isoform"), False)

# gffread source only: keep just the transcripts with a CDS (on by default --
# set false to keep every transcript gffread extracts, ncRNAs included).
CDS_ONLY = _flag(config.get("cds_only"), True)

TMP = config.get("tmp_dir", "tmp")

# The large per-read intermediates (trimmed FASTQ, SAMs) are deleted once every rule that reads
# them has run, unless keep_intermediates is true. Deleted files are only rebuilt if a later run
# needs them again (e.g. after changing max_mismatch, or turning on genome_check).
KEEP_INTERMEDIATES = _flag(config.get("keep_intermediates"), False)


def intermediate(path):
    """temp() unless keep_intermediates: true."""
    return path if KEEP_INTERMEDIATES else temp(path)


# Read-trimming tool: "cutadapt" (default, StructureFold2's fastq_trimmer.py)
# or "fastp" -- see config.yaml's trimmer comment.
TRIMMER = config.get("trimmer", "cutadapt")
if TRIMMER not in ("cutadapt", "fastp"):
    raise ValueError(
        f"config['trimmer'] must be 'cutadapt' or 'fastp', got {TRIMMER!r}"
    )

# 3' adapter trimmed by fastp (trimmer: fastp); default is the Structure-seq
# 3' adapter, overridable with config key fastp_adapter.
FASTP_ADAPTER = str(config.get("fastp_adapter") or "GATCGGAAGAGCACACGTCTG").strip().upper()
if not FASTP_ADAPTER or set(FASTP_ADAPTER) - set("ACGTN"):
    raise ValueError(
        f"config['fastp_adapter'] must be a DNA sequence (A/C/G/T/N), got {FASTP_ADAPTER!r}"
    )

# Minimum read length kept by fastp (trimmer: fastp). Default 20, matching the
# cutadapt trimmer's -m 20 (fastq_trimmer.py); fastp's own default is 15.
try:
    FASTP_MIN_LENGTH = int(config.get("fastp_min_length", 20))
except (TypeError, ValueError):
    raise ValueError(
        f"config['fastp_min_length'] must be an integer, got {config.get('fastp_min_length')!r}"
    )
if FASTP_MIN_LENGTH < 1:
    raise ValueError(f"config['fastp_min_length'] must be >= 1, got {FASTP_MIN_LENGTH}")

# Maximum mismatches/indels per alignment kept by filter_sam (sam_filter.py
# -max_mismatch, sam_processing.smk). Default 3.
try:
    MAX_MISMATCH = int(config.get("max_mismatch", 3))
except (TypeError, ValueError):
    raise ValueError(
        f"config['max_mismatch'] must be an integer, got {config.get('max_mismatch')!r}"
    )
if MAX_MISMATCH < 0:
    raise ValueError(f"config['max_mismatch'] must be >= 0, got {MAX_MISMATCH}")

# Optional 3' coverage-bias correction of the +DMS RT-stops before the
# reactivity calculation
BIAS_CORRECTION = _flag(config.get("3_prime_bias_correction"), False)
BIAS_REFERENCE = config.get("3_prime_bias_reference", "pooled")
BIAS_N_TRANSCRIPTS = int(config.get("3_prime_bias_n_transcripts", 2500))
# Training set: "top_n" (top 3_prime_bias_n_transcripts) or "all_covered"
# (every protein-coding transcript with coverage >= 3_prime_bias_min_coverage
# in every ID's +DMS).
BIAS_TRAINING = config.get("3_prime_bias_training", "top_n")
BIAS_MIN_COVERAGE = float(config.get("3_prime_bias_min_coverage", 1))
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
    if BIAS_TRAINING not in ("top_n", "all_covered"):
        raise ValueError(
            "config['3_prime_bias_training'] must be 'top_n' or 'all_covered', "
            f"got {BIAS_TRAINING!r}"
        )
    if not BIAS_ANNOTATION and not config.get("genotypes"):   # genotype mode: each genotype's own annotation
        raise ValueError(
            "3_prime_bias_correction needs a GTF/GFF3 to find protein-coding "
            "transcripts: set 3_prime_bias_annotation (or annotation_gtf/gff)."
        )

# How each channel's ln(RT-stops + 1) is scaled before the -DMS subtraction (config.yaml's
# reactivity_normalization comment): "all" = mean over all positions (StructureFold2), "gu" = mean
# over the G/U positions. Used by the standard, raw and heat-corrected reactivity
# (rtsc_to_react.py -gu_norm) and, unless depth_matching_normalization overrides it, depth matching.
REACTIVITY_NORMALIZATION = str(config.get("reactivity_normalization", "all"))
if REACTIVITY_NORMALIZATION not in ("all", "gu"):
    raise ValueError("reactivity_normalization must be 'all' or 'gu'")
REACTIVITY_MIN_GU = int(config.get("reactivity_min_gu", 20))
REACTIVITY_NORM_ARGS = f"-gu_norm -min_gu {REACTIVITY_MIN_GU}" if REACTIVITY_NORMALIZATION == "gu" else ""

# Optional depth-matched reactivity (workflow/rules/depth_matching.smk) -- see
# config.yaml's depth_matching comment.
DEPTH_MATCHING = _flag(config.get("depth_matching"), False)
DEPTH_MATCHING_DRAWS = int(config.get("depth_matching_draws", 100))
DEPTH_MATCHING_MATCH_MINUS = _flag(config.get("depth_matching_match_minus"), True)
DEPTH_MATCHING_ROC_TRIM3 = int(config.get("depth_matching_roc_trim3", 50))
DEPTH_MATCHING_SHARDS = int(config.get("depth_matching_shards", 1))
DEPTH_MATCHING_NORMALIZATION = str(config.get("depth_matching_normalization", REACTIVITY_NORMALIZATION))
if DEPTH_MATCHING_NORMALIZATION not in ("all", "gu"):
    raise ValueError("depth_matching_normalization must be 'all' or 'gu'")
DEPTH_MATCHING_MIN_GU = int(config.get("depth_matching_min_gu", REACTIVITY_MIN_GU))
if DEPTH_MATCHING_SHARDS < 1:
    raise ValueError("depth_matching_shards must be >= 1")

# Optional RNA Framework XML copies of the reactivities (rf-norm format, one <transcript>.xml per
# transcript; reactivity.smk) -- see config.yaml's rnaframework_xml comment.
RNAFRAMEWORK_XML = _flag(config.get("rnaframework_xml"), False)

# Canonical transcriptome FASTA used by all downstream rules.
# It is produced by the prepare_transcriptome rule (transcriptome.smk).
TRANSCRIPTOME = f"{TMP}/resources/transcriptome.fa"

# Uppercase-only copy of TRANSCRIPTOME, produced by uppercase_transcriptome
# (transcriptome.smk).
TRANSCRIPTOME_UPPER = f"{TMP}/resources/transcriptome_upper.fa"

# ---------------------------------------------------------------------------
# Optional extra RNA sets appended to the transcriptome before mapping (shared by both genotypes
# in genotype mode): rRNAs, other ncRNAs (tRNA/snoRNA/...) and lncRNAs. Each comes from a FASTA (<kind>_fasta)
# or from a GFF3/GTF annotation (<kind>_gff or <kind>_gtf, either key takes either format) whose
# transcripts gffread extracts from <kind>_genome (default: genome). With an annotation,
# <kind>_noncoding_only (default true) keeps only transcripts without a CDS (gffread --nc), so a
# full annotation can be given. Their reads count as rRNA / ncRNA / lncRNA in the read-composition QC
# (rRNAs inside the ncRNA set are labelled rRNA when its annotation says so).
# ---------------------------------------------------------------------------
# order = order appended to the transcriptome (a duplicated ID keeps its first copy)
EXTRA_RNA_CATEGORY = {"rrna": "rRNA", "ncrna": "ncRNA", "lncrna": "lncRNA"}
EXTRA_RNA = {}          # kind -> {"fasta": FASTA appended, "source": what the user gave, ...}

# Read-composition pie (report.smk): a read whose forward-strand alignments hit several types goes
# to the first of these; reads hitting none of them are "other". Default rRNA > mRNA > lncRNA.
READ_COMPOSITION_PRIORITY = config.get("read_composition_priority") or ["rRNA", "mRNA", "lncRNA"]
if isinstance(READ_COMPOSITION_PRIORITY, str):
    READ_COMPOSITION_PRIORITY = READ_COMPOSITION_PRIORITY.replace(",", " ").split()
for _kind in EXTRA_RNA_CATEGORY:
    _fa = config.get(f"{_kind}_fasta")
    _ann = [k for k in (f"{_kind}_gff", f"{_kind}_gtf") if config.get(k)]
    if _fa and _ann:
        raise ValueError(f"give {_kind}_fasta or {_ann[0]}, not both")
    if len(_ann) > 1:
        raise ValueError(f"give one of {_kind}_gff / {_kind}_gtf, not both")
    if _fa:
        EXTRA_RNA[_kind] = {"fasta": _fa, "source": _fa}
    elif _ann:
        _genome = config.get(f"{_kind}_genome") or config.get("genome")
        if not _genome:
            raise ValueError(
                f"{_ann[0]} needs a genome FASTA to extract sequences from: set {_kind}_genome "
                "(or genome)")
        EXTRA_RNA[_kind] = {
            "fasta": f"{TMP}/resources/{_kind}_gffread.fa",
            "annotation": config[_ann[0]],
            "genome": _genome,
            "noncoding_only": _flag(config.get(f"{_kind}_noncoding_only"), True),
            "source": f"{config[_ann[0]]} (gffread from {_genome})",
        }

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
# p4p6 ROC (plot_p4p6_roc): leave out the last N (3') positions. Default 50.
P4P6_ROC_TRIM3 = int(config.get("p4p6_roc_trim3", 50))

# ---------------------------------------------------------------------------
# Optional: two genotypes (species, or individuals of one species), each with its own genome and
# GFF/GTF (config key genotypes). Each genotype's transcripts are extracted with gffread from its
# own genome and renamed "<genotype>__<transcript_id>", so both live in one combined transcriptome
# (TRANSCRIPTOME) downstream without ID clashes; reads map only to their own genotype's transcripts
# (a per-genotype Bowtie2 index). The positive control and ncrna_fasta are shared and keep their
# names. Transcripts are matched between genotypes by identical ID, or by transcript_matching
# (TSV, one column per genotype name). See workflow/rules/genotypes.smk.
# ---------------------------------------------------------------------------
GENOTYPES = config.get("genotypes") or {}
GENOTYPE_MODE = bool(GENOTYPES)
GENOTYPE_NAMES = list(GENOTYPES)
GENOTYPE_SEP = "__"
GENOTYPE_DIR = f"{TMP}/resources/genotypes"
TRANSCRIPT_MATCHING = config.get("transcript_matching")
POOL_MINUS_ACROSS_GENOTYPES = _flag(config.get("pool_minus_across_genotypes"), False)
# Combined -DMS for every reactivity when pool_minus_across_genotypes is on (genotypes.smk)
CROSS_GENOTYPE_MINUS = f"{TMP}/output/se/pooled/combined_minus_cross_genotype.rtsc"
GENOTYPE_LIFTOVER = f"{config['output_dir']}/genotype_comparison/matched_transcripts.tsv"
# Optional variant-derived position map for the cross-genotype -DMS pooling (replaces the alignment)
CROSS_GENOTYPE_COORDINATE_MAP = config.get("cross_genotype_coordinate_map")
CROSS_GENOTYPE_MAP_REF = config.get("cross_genotype_coordinate_map_ref")
if GENOTYPE_MODE:
    if len(GENOTYPE_NAMES) != 2:
        raise ValueError(f"config['genotypes'] must name exactly 2 genotypes, got {GENOTYPE_NAMES}")
    for _g, _spec in GENOTYPES.items():
        if GENOTYPE_SEP in _g or not _g.replace("_", "").replace("-", "").isalnum():
            raise ValueError(f"genotype name {_g!r} must be letters/digits/-/_ and must not contain '{GENOTYPE_SEP}'")
        if not isinstance(_spec, dict):
            raise ValueError(f"config['genotypes'][{_g!r}] needs 'genome' and 'gtf' or 'gff'")
        if _spec.get("gff") and _spec.get("gtf") and _spec["gff"] != _spec["gtf"]:
            raise ValueError(f"config['genotypes'][{_g!r}] sets both 'gff' and 'gtf'; give one annotation")
        _spec["gff"] = _spec.get("gff") or _spec.get("gtf")   # 'gtf' is an alias: gffread reads GFF3 and GTF
        if not _spec.get("genome") or not _spec["gff"]:
            raise ValueError(f"config['genotypes'][{_g!r}] needs both 'genome' and an annotation ('gtf' or 'gff')")
    if TRANSCRIPT_SOURCE != "gffread":
        TRANSCRIPT_SOURCE = "gffread"   # each genotype's transcripts always come from its genome + GFF
    if CROSS_GENOTYPE_COORDINATE_MAP:
        CROSS_GENOTYPE_MAP_REF = CROSS_GENOTYPE_MAP_REF or GENOTYPE_NAMES[0]
        if CROSS_GENOTYPE_MAP_REF not in GENOTYPE_NAMES:
            raise ValueError(f"cross_genotype_coordinate_map_ref {CROSS_GENOTYPE_MAP_REF!r} must be one of {GENOTYPE_NAMES}")
        if not POOL_MINUS_ACROSS_GENOTYPES:
            raise ValueError("cross_genotype_coordinate_map is only used with pool_minus_across_genotypes: true")
    if POOL_MINUS_ACROSS_GENOTYPES and config.get("pool_replicates", "none") != "minus":
        raise ValueError(
            "pool_minus_across_genotypes needs pool_replicates: minus (every -DMS sample of both "
            "genotypes is pooled, then each matched transcript also gets the other genotype's -DMS)."
        )
elif TRANSCRIPT_MATCHING or POOL_MINUS_ACROSS_GENOTYPES:
    raise ValueError("transcript_matching / pool_minus_across_genotypes need config['genotypes']")


def genotype_annotation(genotype):
    """Annotation used for genotype's protein-coding list / position annotations (annotation_gtf, else gff)."""
    spec = GENOTYPES[genotype]
    return spec.get("annotation_gtf") or spec["gff"]


# ---------------------------------------------------------------------------
# Optional: map the reads Bowtie2 could not align to the transcriptome to the genome (HISAT2,
# splice-aware) and report where they land (config key genome_check; workflow/rules/genome_check.smk).
# One reference per genome: "genome" (config genome), or each genotype's own genome in genotype
# mode. Annotation used to place the reads: genome_check_annotation (or a genotype's
# genome_check_annotation), else annotation_gtf, else gff; none = chromosome-level only.
# genome_check_hisat2_index (or a genotype's hisat2_index): an existing HISAT2 index prefix to use
# instead of building one.
# ---------------------------------------------------------------------------
GENOME_CHECK = _flag(config.get("genome_check"), False)
GENOME_CHECK_DIR = f"{config['output_dir']}/genome_check"
GENOME_CHECK_TMP = f"{TMP}/genome_check"
# Chromosomes counted as organelle genomes (case-insensitive regex on the whole name)
GENOME_CHECK_ORGANELLES = config.get(
    "genome_check_organelles",
    r"(chr)?(Mt|M|MT|mito.*|mitochondri.*|Pt|C|chloroplast.*|plastid.*|plastome.*)",
)
# Read 5' ends of at most this many nt apart (same chromosome and strand) form one locus
GENOME_CHECK_LOCUS_GAP = int(config.get("genome_check_locus_gap", 50))
GENOME_CHECK_TOP_LOCI = int(config.get("genome_check_top_loci", 200))
GENOME_CHECK_REFS = {}   # reference name -> {"genome", "annotation", "index"}
if GENOME_CHECK:
    if GENOTYPE_MODE:
        for _g, _spec in GENOTYPES.items():
            GENOME_CHECK_REFS[_g] = {
                "genome": _spec["genome"],
                "annotation": _spec.get("genome_check_annotation") or genotype_annotation(_g),
                "index": _spec.get("hisat2_index", ""),
            }
    else:
        if not config.get("genome") and not config.get("genome_check_hisat2_index"):
            raise ValueError("genome_check needs config['genome'] (or genome_check_hisat2_index)")
        GENOME_CHECK_REFS["genome"] = {
            "genome": config.get("genome", ""),
            "annotation": config.get("genome_check_annotation")
            or config.get("annotation_gtf") or config.get("gff") or "",
            "index": config.get("genome_check_hisat2_index", ""),
        }
