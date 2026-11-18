###############################################################################
# structure-fold2-snakemake
#
# StructureFold2 as a Snakemake workflow: Structure-seq / Structure-seq2 reads (DMS, RT stops)
# -> per-transcript, per-position reactivity, plus QC. Setup and options: README.md.
#
#   snakemake -n --configfile config/my_config.yaml                 # dry run: list the jobs
#   snakemake --sdm conda --cores 8 --configfile config/my_config.yaml
###############################################################################

from snakemake.shell import shell

shell.prefix(
    "set -euo pipefail; "
    "module load gcc >/dev/null 2>&1 || true; "
    'if [ -n "${{CONDA_PREFIX:-}}" ]; then '
    'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"; '
    "fi; "
)


configfile: "config/config.yaml"


# Used with --sdm conda apptainer: each conda environment is built inside this image
container: "docker://continuumio/miniconda3"


# Config options, samplesheet, and shared helpers (no rules)
include: "workflow/rules/settings.smk"
include: "workflow/rules/samples.smk"
include: "workflow/rules/helpers.smk"
include: "workflow/rules/targets.smk"


# Defined before the other rules so it is the default target.
rule all:
    input:
        get_all_targets,


# Pipeline rules
include: "workflow/rules/transcriptome.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/read_prep.smk"
include: "workflow/rules/sam_processing.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/rtsc.smk"
include: "workflow/rules/reactivity.smk"
include: "workflow/rules/heat_correction.smk"
include: "workflow/rules/coverage.smk"
include: "workflow/rules/abundance.smk"
include: "workflow/rules/specificity.smk"
include: "workflow/rules/replicate_correlation.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/metagene.smk"


# Two genotypes (config key genotypes): per-genotype transcriptomes/indexes, transcript matching,
# cross-genotype -DMS pooling and the reactivity comparison
if GENOTYPE_MODE:
    include: "workflow/rules/genotypes.smk"


# Optional modules
if BIAS_CORRECTION:
    include: "workflow/rules/bias_correction.smk"

if POSITIVE_CONTROL_IS_P4P6:
    include: "workflow/rules/positive_control.smk"

if DEPTH_MATCHING:
    include: "workflow/rules/depth_matching.smk"

if GENOME_CHECK:
    include: "workflow/rules/genome_check.smk"


# Read composition by target type + the self-contained HTML QC report (qc/report.html)
include: "workflow/rules/report.smk"
