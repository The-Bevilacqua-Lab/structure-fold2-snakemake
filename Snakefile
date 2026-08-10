###############################################################################
# structure-fold2-snakemake
#
# A generalized Snakemake wrapper around StructureFold2 for DMS-MaPseq /
# Structure-seq-style +DMS/-DMS reactivity calculation on any organism.
#
# Minimal usage:
#   snakemake --sdm conda -j 8 --configfile config/config.yaml
#
# On a SLURM cluster:
#   snakemake --sdm conda -j 100 --executor cluster-generic \
#     --cluster-generic-submit-cmd 'sbatch --time=4:00:00 --ntasks=10 --mem=40gb \
#       --output=slurm_logs/%j.out' \
#     --configfile config/config.yaml --latency-wait 60
#
# See config/config.yaml for every available config key, and README.md for
# an overview of what this pipeline does and does not generalize.
###############################################################################

from snakemake.utils import min_version
import pandas as pd
from snakemake.shell import shell

# On HPC systems that use Environment Modules, a loaded module (e.g. an R or
# Python module) can prepend its own lib dir onto LD_LIBRARY_PATH, shadowing
# the newer libstdc++ shipped inside conda envs and breaking compiled
# extensions. `module load gcc` first (a no-op, and harmless, on any system
# without an Environment Modules setup -- the `|| true` swallows the "module:
# command not found" case) gets a modern libstdc++ onto LD_LIBRARY_PATH, then
# the conda env's own lib dir (once activated) takes final precedence.
# Re-declaring "set -euo pipefail;" here is required because shell.prefix()
# overrides Snakemake's default rather than extending it.
shell.prefix(
    "set -euo pipefail; "
    "module load gcc >/dev/null 2>&1 || true; "
    'if [ -n "${{CONDA_PREFIX:-}}" ]; then '
    'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"; '
    "fi; "
)

configfile: "config/config.yaml"

container: "docker://continuumio/miniconda3"

# ---------------------------------------------------------------------------
# Samplesheet
#
# One row per (sample, sequencing run) pair, tab-separated, with columns:
#   sample     – sample identifier (matches fastq naming)
#   condition  – exactly "plus" (+DMS/treated) or "minus" (-DMS/untreated)
#   r1         – path to that (sample, run)'s raw FASTQ
#   ID         – optional; groups samples into biological replicates (e.g.
#                "rep1") -- the reactivity calculation pairs one condition's
#                ID against the other's, so a given ID needs both a "plus"
#                and a "minus" row somewhere in the sheet. If the column is
#                absent, ID defaults to each sample's row position *within
#                its condition* ("rep1", "rep2", ... independently for plus
#                and minus), i.e. the 1st minus row is paired with the 1st
#                plus row as one replicate, the 2nd with the 2nd, and so on
#                -- this assumes equal plus/minus sample counts in that
#                order; give an explicit ID column if that's not your layout.
#   run        – optional; distinguishes multiple sequencing runs of the same
#                sample (resequenced on a second flow cell, etc.), which get
#                concatenated before trimming. Defaults to a single run if
#                the column is absent.
# ---------------------------------------------------------------------------
samples_runs = pd.read_csv(config["samplesheet"], sep="\t")
if "run" not in samples_runs.columns:
    samples_runs["run"] = "run1"
if "ID" not in samples_runs.columns:
    # Compute once per unique sample (not per run-row, which would
    # over-count multi-run samples), ranked by row order within each
    # condition, then broadcast back across that sample's run-rows.
    sample_condition = samples_runs.drop_duplicates(subset="sample")[["sample", "condition"]].copy()
    sample_condition["ID"] = "rep" + (sample_condition.groupby("condition").cumcount() + 1).astype(str)
    samples_runs = samples_runs.merge(sample_condition[["sample", "ID"]], on="sample", how="left")

samples = samples_runs.drop_duplicates(subset="sample", keep="first")

SAMPLES = samples["sample"].tolist()
MINUS_SAMPLES = samples.loc[samples["condition"] == "minus", "sample"].tolist()

# See workflow/rules/structurefold2.smk for the full explanation of what
# each pool_replicates value does to the reactivity calculation.
POOL_REPLICATES = config.get("pool_replicates", "none")

# "both" collapses every replicate's -DMS and +DMS into one pooled pair, so
# there is exactly one reactivity output instead of one per replicate ID --
# every per-ID rule in the pipeline (coverage, specificity, QC plots...)
# then runs once, on that single pooled dataset, under the synthetic ID
# "pooled". Any other pool_replicates value leaves IDS as the real
# replicate-group list.
if POOL_REPLICATES == "both":
    IDS = ["pooled"]
else:
    IDS = samples["ID"].unique().tolist()

include: "workflow/rules/common.smk"
include: "workflow/rules/structurefold2.smk"
include: "workflow/rules/annotation.smk"

# The p4p6 structure plot rules depend on a hand-annotated coordinate map
# specific to that one construct (see workflow/rules/positive_control.smk) --
# only wired in when the positive control actually is p4p6. Any other
# positive control still gets the generic alignment-% QC from common.smk.
POSITIVE_CONTROL_IS_P4P6 = config.get("positive_control_name") == "p4p6"
if POSITIVE_CONTROL_IS_P4P6:
    include: "workflow/rules/positive_control.smk"


def get_all_targets(wildcards):
    out = config["output_dir"]
    targets = (
        expand("{out}/{id}/reactivity.react", out=out, id=IDS) +
        expand("{out}/{id}/reactivity.csv", out=out, id=IDS) +
        expand("{out}/{id}/coverage.csv", out=out, id=IDS) +
        expand("{out}/{id}/counts_minus.csv", out=out, id=IDS) +
        expand("{out}/{id}/specificity_plus.csv", out=out, id=IDS) +
        expand("{out}/{id}/specificity_minus.csv", out=out, id=IDS) +
        expand("{out}/specificity_{sample}.csv", out=out, sample=SAMPLES) +
        [
            f"{out}/qc/alignment_stats_summary.tsv",
            f"{out}/qc/specificity_plot.png",
        ] +
        (
            [
                f"{out}/qc/rt_stop_replicate_correlation.png",
                f"{out}/qc/rt_stop_replicate_correlation_minus.png",
            ]
            if len(IDS) > 1
            else []
        ) +
        (
            [f"{out}/transcript_position_annotations.tsv"]
            if config.get("annotation_gtf") and config.get("genome")
            else []
        )
    )

    if config.get("positive_control_fasta"):
        targets += [f"{out}/qc/positive_control_alignment_summary.tsv"]

    if POSITIVE_CONTROL_IS_P4P6:
        targets += (
            expand("{out}/{id}/p4p6_react.png", out=out, id=IDS) +
            expand("{out}/{id}/p4p6_react_plus_only.png", out=out, id=IDS) +
            [
                f"{out}/qc/p4p6_all_samples.png",
                f"{out}/qc/p4p6_all_samples_plus_only.png",
            ]
        )

    if config.get("specificity_comparison_runs"):
        targets += [f"{out}/qc/specificity_comparison.png"]

    return targets


rule all:
    input:
        get_all_targets,
