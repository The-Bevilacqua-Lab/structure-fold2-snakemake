# Choosing options

The defaults reproduce a standard StructureFold2 analysis. This page says when to change them.
Every key is listed, with its default, in `config/config.yaml`.

## Coming from StructureFold2

Each StructureFold2 step is a rule here. The scripts are Python 3 ports with the same defaults and
the same output (`workflow/scripts/StructureFold3/`).

| StructureFold2 step | Rule(s) | Notes |
|---|---|---|
| `fastq_trimmer.py` | `trim_reads_cutadapt` | Or `trimmer: fastp` |
| `fastq_mapper.py` (Bowtie2 `-a`) | `align_with_bowtie2` | Unmapped reads are kept in the SAM for QC; `sam_filter` drops them |
| `sam_filter.py` | `filter_sam` | `max_mismatch` (default 3) |
| `sam_to_rtsc.py` | `sam_to_rtsc` | One `.rtsc` per sample |
| `rtsc_combine.py` | `combine_rtsc` | Per replicate ID, or pooled (`pool_replicates`) |
| `rtsc_to_react.py` | `rtsc_to_react` | One `reactivity.react` + `.csv` per ID |
| `rtsc_coverage.py`, `rtsc_specificity.py`, `rtsc_correlation.py`, `rtsc_abundances.py` | QC rules | All in `qc/report.html` |
| `react_heat_correct.py` | heat-correction rules | Turned on by a `temperature` samplesheet column |

What you no longer do by hand: put files in the right folder, rename outputs, pair each +DMS file
with its −DMS file. The samplesheet does the pairing: samples with the same `ID` form one
replicate, and each ID gets one reactivity.

## Snakemake in five commands

Run from the repository folder.

| Task | Command |
|---|---|
| See what would run, without running it | `snakemake -n --configfile my.yaml` |
| Run on this machine with 8 cores | `snakemake --sdm conda --cores 8 --configfile my.yaml` |
| Resume after a crash or a cancelled job | the same command plus `--rerun-incomplete` |
| "Directory cannot be locked" after a crash | the same command plus `--unlock`, then run again |
| Rerun only for changed input files, not for changed code or settings | add `--rerun-triggers mtime` |

Snakemake works backwards from the final files: it runs a step only if its output is missing or
older than its inputs (and, by default, if the step's code or settings changed). So changing one
option reruns only the steps it affects. Each job writes a log under `logs/<rule>/`; when a job
fails, the error message names that log.

On a SLURM cluster, see [README §3](../README.md#3-run).

## Decisions, in order

### Reference: `transcript_source`

- **`file`** if you already have a transcriptome FASTA.
- **`gffread`** to build it from a genome + GFF3/GTF. By default only transcripts with a CDS are kept
  (`cds_only: true`); add `keep_longest_isoform: true` for one isoform per gene, which reduces
  multimapping between isoforms.
- **Add an ncRNA set** (`ncrna_fasta` or `ncrna_gtf`) when the transcriptome is coding-only.
  Otherwise rRNA reads can land on similar mRNA sequence.
- **Add the spike-in** (`positive_control_fasta`) if the library has one. With p4p6 you get a ROC
  against the known structure, the best single check that the reactivities are right.

### Trimming: `trimmer`

- **`cutadapt`** (default) for Structure-seq2 libraries: StructureFold2's trimming, both adapters.
- **`fastp`** for a different 3′ adapter (`fastp_adapter`), or for NextSeq/NovaSeq data with poly-G
  tails. fastp does not trim a 5′ adapter.

### Replicates: `pool_replicates`

| You want | Use |
|---|---|
| One reactivity per replicate, e.g. to test differences between conditions (dStruct, deltaSHAPE) or to check reproducibility | `none` (default) |
| The same, but the −DMS libraries are shallow or noisy | `minus`: one pooled −DMS background for every replicate |
| One best reactivity profile from all the data | `both` |
| One profile per temperature (with a `temperature` column) | `both_by_temperature` |

Replicate correlation QC always uses the unpooled replicates, whatever you choose. To drop one bad
library without deleting it from the samplesheet, use `exclude_samples_plus` / `exclude_samples_minus`.

### Reactivity: `reactivity_normalization`

- **`all`** (default) is the StructureFold2 calculation. Use it to match earlier StructureFold2
  results.
- **`gu`** scales each channel by its mean at G/U, which DMS does not modify, instead of at all
  positions. The all-position mean of +DMS includes the A/C signal, so `all` subtracts more −DMS
  background from more reactive libraries. Prefer `gu` when you compare libraries, conditions or
  genotypes. On the rice data it gave better replicate and genotype agreement and the same p4p6
  AUC. Details: [methods.md §8](methods.md#gu-normalization-reactivity_normalization-gu).

### Optional corrections

| Turn on | When | Read more |
|---|---|---|
| `3_prime_bias_correction: true` | Poly(A)-selected or oligo-dT-primed RNA, and the report's depth metagene shows a different 5′→3′ slope for +DMS than for −DMS | [3_prime_end_bias.md](3_prime_end_bias.md) |
| `depth_matching: true` | You compare libraries (conditions, genotypes) whose depth differs; depth alone changes reactivity | [depth_matching.md](depth_matching.md) |
| `temperature` samplesheet column | The same RNA probed at two temperatures | [pipeline.md §10](pipeline.md#10-optional-heat-correction) |
| `genotypes` | Two genomes (species or cultivars), each read mapped to its own genotype | [pipeline.md §11G](pipeline.md#11g-optional-two-genotypes) |
| `genome_check: true` | Many reads fail to map to the transcriptome and you want to know why | [pipeline.md §12.10](pipeline.md#1210-optional-transcriptome-unaligned-reads-on-the-genome--workflowrulesgenome_checksmk) |

## Which reactivity file to use

| Enabled | Use |
|---|---|
| Nothing optional, or 3′ bias correction | `{id}/reactivity.csv` |
| Depth matching | `depth_matched/{id}/reactivity.csv` to compare libraries; `{id}/reactivity.csv` for one library alone |
| Heat correction | `{id}/heat_correction/reactivity_corrected.csv` |

Trust a transcript only if its coverage (`{id}/coverage.csv`, +DMS A/C stops per A/C base) is at
least 1; that is the threshold the heat correction and most StructureFold2 analyses use.

## Example configs

A standard run, rice in vitro with the p4p6 spike-in:

```yaml
samplesheet: config/my_samples.tsv
output_dir: results
tmp_dir: /scratch/me/sf2_tmp
transcript_source: gffread
genome: data/genome.fa
gff: data/annotation.gff3
ncrna_gtf: data/full_annotation.gtf
positive_control_fasta: resources/p4p6/p4p6.fa
positive_control_name: p4p6
```

Comparing two conditions of different depth: add

```yaml
reactivity_normalization: gu
depth_matching: true
```
