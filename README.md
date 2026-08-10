# structure-fold2-snakemake

## Overview
A generalized Snakemake pipeline around [StructureFold2](https://github.com/StructureFold2/StructureFold2)
to get reactivities from a Structure-seq experiment. Given
raw/trimmed single-end reads and a transcriptome (or a genome + annotation to
build one), it aligns with Bowtie2, counts RT-stops per replicate, and
calculates per-transcript, per-position reactivity. 

# Running the pipeline
## Setting up the environment

To run this pipeline, you need Snakemake installed. You can create and activate a conda environment 
with containing Snakemake with the following command. 

```
conda env create -f environment.yml
conda activate structure-fold2-snakemake
```
The pipeline itself handles all of the other software dependencies so don't worry about that. 


## Quick start

1. Copy `config/config.yaml`, fill in your samplesheet/transcriptome paths
   (every key is documented inline).
2. Build a samplesheet like `config/samplesheet_example.tsv` (or the
   3-column `config/samplesheet_minimal_example.tsv` if you don't need
   replicate grouping or multi-run concatenation -- see below).
3. `snakemake --sdm conda -j 8 --configfile config/your_config.yaml`


## Samplesheet

Tab-separated, one row per (sample, sequencing run):

| column | required | meaning |
|---|---|---|
| `sample` | yes | sample identifier, matches your FASTQ naming |
| `condition` | yes | exactly `plus` (+DMS/treated) or `minus` (-DMS/untreated) |
| `r1` | yes | path to that (sample, run)'s raw FASTQ |
| `ID` | no | groups samples into biological replicates (e.g. `rep1`) -- a given ID needs both a `plus` and `minus` row somewhere in the sheet, since the reactivity calculation pairs one condition's ID against the other's. **If omitted**, defaults to each sample's row position *within its condition* (`rep1`, `rep2`, ... independently for `plus` and `minus`) -- i.e. the 1st `minus` row pairs with the 1st `plus` row as one replicate, the 2nd with the 2nd, and so on. This assumes equal `plus`/`minus` counts given in matching order; give an explicit `ID` column if that's not your layout. |
| `run` | no | distinguishes multiple sequencing runs of the same sample (e.g. resequenced on a second flow cell); same-sample rows get concatenated before trimming. **Defaults to a single run** if omitted. |


## Outputs (under `output_dir`)

Per replicate ID (or the single `pooled` ID under `pool_replicates: both`):
- `{id}/reactivity.react`, `.csv` -- +DMS/-DMS subtracted, 2-8%-normalized reactivity
- `{id}/coverage.csv`, `{id}/specificity_{plus,minus}.csv`, `{id}/counts_minus.csv` -- per-replicate QC
- `{id}/p4p6_react*.png` -- only when `positive_control_name: p4p6`

Pipeline-wide, under `qc/`: alignment stats, replicate correlation (when
there's more than one replicate), base specificity, and (when a positive
control is configured) its alignment-%.

## Testing

`tests/` has a separate, self-contained Snakemake pipeline that downloads a
real (subsampled) public Structure-seq2 dataset and transcriptome, runs
this pipeline against them end-to-end, and sanity-checks the results:

```
snakemake --sdm conda -j 8 -s tests/Snakefile
```

Run from the repo root. See `tests/README.md` for what data it uses, what
"passing" actually checks, and expected runtime/disk usage (it's a real
integration test against ~12 million reads, not a quick unit test).
