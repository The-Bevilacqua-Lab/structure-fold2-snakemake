# Pipeline walkthrough

This document goes through every step of structure-fold2-snakemake in the
order data flows through it, from raw FASTQ to reactivity and QC. For each
step it lists the Snakemake rule(s), the rule file they live in, the
script or tool each one runs, and the files that go in and come out.
What the steps compute (RT stops, filtering, the reactivity formula) is in
[methods.md](methods.md); setup and running are in the [README](../README.md).

- Rule files: `workflow/rules/*.smk` (pulled in by `Snakefile`)
- Scripts: `workflow/scripts/` (StructureFold2's Python 3 ports are in
  `workflow/scripts/StructureFold3/`)
- Conda environments: `workflow/envs/*.yaml`

## Contents

0. [Setup: config, samplesheet, targets](#0-setup-config-samplesheet-targets)
1. [Transcriptome preparation](#1-transcriptome-preparation)
2. [Bowtie2 index](#2-bowtie2-index)
3. [Read trimming](#3-read-trimming)
4. [Alignment and merging runs](#4-alignment-and-merging-runs)
5. [SAM filtering and BAM conversion](#5-sam-filtering-and-bam-conversion)
6. [RT-stop counting](#6-rt-stop-counting)
7. [Combining replicates / pooling](#7-combining-replicates--pooling)
8. [Optional: 3′ coverage-bias correction](#8-optional-3-coverage-bias-correction)
9. [Reactivity calculation](#9-reactivity-calculation)
10. [Optional: heat correction](#10-optional-heat-correction)
11. [Optional: depth-matched reactivity](#11-optional-depth-matched-reactivity)
    - [11G. Optional: two genotypes](#11g-optional-two-genotypes)
12. [QC and summary outputs](#12-qc-and-summary-outputs)
    - [12.10 Optional: transcriptome-unaligned reads on the genome](#1210-optional-transcriptome-unaligned-reads-on-the-genome--workflowrulesgenome_checksmk)
13. [Final outputs (`rule all`)](#13-final-outputs-rule-all)
14. [Script index](#14-script-index)
15. [Known issues](#15-known-issues)

---

## Overview

```
                  samplesheet.tsv + config.yaml
                              │
 genome+GFF ──(AGAT)──► gffread ──► filter ──┐
 or transcriptome.fa ────────────────────────┤
 + positive control / ncRNA FASTAs ──────────┤
                                             ▼
                              prepare_transcriptome ──► uppercase_transcriptome
                                             │
                                    build_bowtie2_index
                                             │
 raw FASTQ (per sample, per run)             │
   └► trim (cutadapt|fastp) ─────────────► align_with_bowtie2
                                                     │
                               combine_sam_from_multiple_runs   (per sample)
                                                     │
                            filter_sam ──────► sorted BAM
                           │                          │
                      sam_to_rtsc            (depth, metagene, bias model)
                           │
              combine_rtsc   (+DMS and −DMS, per replicate ID, or pooled)
                           │
            [optional 3′ bias correction of +DMS]
                           │
                     rtsc_to_react ─► reactivity.react / reactivity.csv
                           │
     [optional heat correction]   [optional depth matching]   QC / plots
```

### Wildcards and directories

| Wildcard      | Meaning                                                                 |
|---------------|-------------------------------------------------------------------------|
| `{sample}`    | One row group in the samplesheet's `sample` column (one library)        |
| `{run}`       | Sequencing run (samplesheet `run` column); a sample may have several    |
| `{id}`        | Reactivity ID: a replicate group (`ID` column), or `pooled` / `pooled_<temperature>` when pooling |
| `{condition}` / `{treatment}` | `plus` (+DMS) or `minus` (-DMS)                         |

- `{TMP}` = `config["tmp_dir"]` (default `tmp`): intermediate files
  (trimmed reads, SAM/BAM, per-sample and combined `.rtsc`). Can be scratch. The trimmed FASTQs
  and the SAMs are deleted as soon as every rule that reads them has finished (Snakemake
  `temp()`), unless `keep_intermediates: true`; a later run rebuilds them only if it needs them.
- `{out}` = `config["output_dir"]`: final reactivity, QC tables and plots.
- `logs/<rule>/...`: one log per rule invocation.

---

## 0. Setup: config, samplesheet, targets

No jobs run here; these files are read when Snakemake parses the workflow
and decide which rules exist and which outputs are requested.

| File | What it does |
|------|--------------|
| `Snakefile` | Sets a shell prefix (`set -euo pipefail`, `module load gcc`, conda `LD_LIBRARY_PATH`), loads `config/config.yaml`, includes every rule file, and includes the optional modules (`bias_correction.smk`, `positive_control.smk`, `depth_matching.smk`) only when enabled. |
| `workflow/rules/settings.smk` | Reads and checks the config options (`transcript_source`, `tmp_dir`, `keep_intermediates`, `trimmer`, `reactivity_normalization`, `3_prime_bias_*`, `depth_matching*`, `genotypes`, `genome_check` …). Defines `TRANSCRIPTOME` (`{TMP}/resources/transcriptome.fa`) and `TRANSCRIPTOME_UPPER`. |
| `workflow/rules/samples.smk` | Parses and checks the samplesheet (required columns, `plus`/`minus`, the `ID` and `run` defaults) and builds the sample/ID lists: `SAMPLES`, `IDS`, `REPLICATE_IDS`, `IDS_WITH_PLUS/MINUS`, `CORRELATION_PLUS/MINUS_IDS`, the `exclude_samples_plus/minus` sets, and the heat-correction temperature groups (enabled when the sheet has a `temperature` column, which must have exactly 2 values). Sets `IDS` to `["pooled"]` or `["pooled_<t>", ...]` for the pooling modes that collapse replicates. |
| `workflow/rules/helpers.smk` | Input functions used by rules: which FASTQ/runs belong to a sample, which `.rtsc` files feed each combine step, which `.rtsc` (raw, pooled, or bias-corrected) feeds reactivity, which BAMs feed the bias model, and the extra FASTAs appended to the transcriptome. |
| `workflow/rules/targets.smk` | `get_all_targets()`: the list of final files `rule all` asks for (see [§13](#13-final-outputs-rule-all)). |

See [samplesheet.md](samplesheet.md) for the samplesheet columns.

---

## 1. Transcriptome preparation

**Rule file:** `workflow/rules/transcriptome.smk`

Builds the single reference FASTA that everything maps to. There are two ways
to supply it (`transcript_source`):

- `file`: use `config["transcriptome"]` directly; skip to `prepare_transcriptome`.
- `gffread`: build it from `genome` + `gff`.

| Order | Rule | Runs when | Tool / script | Input → Output |
|------|------|-----------|---------------|----------------|
| 1a | `agat_longest_isoform` | `transcript_source: gffread` and `keep_longest_isoform: true` | `agat_sp_keep_longest_isoform.pl` (env `agat.yaml`) | `gff` → `{TMP}/resources/longest_isoform.gff` |
| 1b | `gffread_extract_transcripts` | `transcript_source: gffread` | `gffread -w` (env `gffread.yaml`) | `genome` + GFF (original or longest-isoform) → `{TMP}/resources/transcriptome_gffread.fa` |
| 1c | `filter_extracted_seqs` | `transcript_source: gffread` and `cds_only: true` (default) | `workflow/scripts/transcriptome/filter_extracted_seqs.py` (env `biopython.yaml`) | Keeps only records whose header mentions a `CDS` → `{TMP}/resources/transcriptome_gffread_filtered.fa` |
| 1d | `prepare_transcriptome` | always | `awk 1` concatenation | Base transcriptome (user FASTA; filtered gffread FASTA; or unfiltered `transcriptome_gffread.fa` when `cds_only: false`) + `positive_control_fasta` + ncRNA set + lncRNA set (whichever are set; a FASTA, or `{TMP}/resources/{ncrna,lncrna}_gffread.fa` from rule `gffread_extract_extra_rna` when given as GFF/GTF) → `{TMP}/resources/transcriptome.fa` (`TRANSCRIPTOME`). A record whose ID already appeared is dropped (`concat_fasta_unique.awk`, logged) |
| 1e | `uppercase_transcriptome` | always | `workflow/scripts/transcriptome/uppercase_fasta.py` (env `biopython.yaml`) | `TRANSCRIPTOME` → `{TMP}/resources/transcriptome_upper.fa` (`TRANSCRIPTOME_UPPER`) |

The uppercase copy matters: Bowtie2 doesn't care about case, but every
downstream "is this base A or C?" check (reactivity, specificity,
annotation) compares against uppercase letters, so soft-masked (lowercase)
bases would otherwise be skipped silently. Bowtie2 maps to `TRANSCRIPTOME`;
every rule after mapping reads `TRANSCRIPTOME_UPPER`.

---

## 2. Bowtie2 index

**Rule file:** `workflow/rules/alignment.smk`

| Rule | Tool | Input → Output |
|------|------|----------------|
| `build_bowtie2_index` | `bowtie2-build` (env `bowtie2.yaml`) | `TRANSCRIPTOME` → `TRANSCRIPTOME.*.bt2` (the rule tracks `.1.bt2`) |

---

## 3. Read trimming

**Rule file:** `workflow/rules/read_prep.smk`

Each (sample, run) pair is processed **separately** through trimming and
alignment. Runs are only merged after alignment (step 4), so a sample
sequenced on several runs ends up as one SAM. Both trimmers read the samplesheet's `r1`
directly, plain or gzipped.

| Rule | Runs when | Tool / script | Input → Output |
|------|-----------|---------------|----------------|
| `trim_reads_cutadapt` | `trimmer: cutadapt` (default) | `workflow/scripts/StructureFold3/fastq_trimmer.py -input -output` (Python 3 port of StructureFold2's trimmer, cutadapt-based; env `cutadapt.yaml`) | `r1` → `{TMP}/reads/trimmed/{sample}_{run}_trimmed.fastq` |
| `trim_reads_fastp` | `trimmer: fastp` | snakemake-wrapper `v3.3.3/bio/fastp` | `r1` → trimmed FASTQ + `_fastp.html`, `_fastp.json` |

Trimmer settings:

- **cutadapt / `fastq_trimmer.py`**: 3′ adapter
  `GATCGGAAGAGCACACGTCTG`, 5′ adapter `TGAACAGCGACTAGGCTCTTCA`, minimum length
  20, minimum quality 30.
- **fastp**: `--adapter_sequence <fastp_adapter> --trim_poly_g --length_required
  <fastp_min_length>` (`fastp_adapter` defaults to `GATCGGAAGAGCACACGTCTG`,
  `fastp_min_length` to 20, matching cutadapt), plus
  fastp's default quality filtering.

---

## 4. Alignment and merging runs

**Rule file:** `workflow/rules/alignment.smk`

| Order | Rule | Tool / script | Input → Output |
|------|------|---------------|----------------|
| 4a | `align_with_bowtie2` | `bowtie2 -a -p 4 -x <index> -q` (env `bowtie2.yaml`) | Trimmed FASTQ → `{TMP}/output/se_runs/{sample}_{run}/{sample}_{run}_trimmed_mapped.sam` + `..._mapping_log.txt` (bowtie2's stderr summary) |
| 4b | `combine_sam_from_multiple_runs` | `samtools view -H` + `grep` (one run: a hard link); `workflow/scripts/alignment/combine_mapping_logs.py` | All runs' SAMs for one sample → `{TMP}/output/se/{sample}/{sample}_trimmed_mapped.sam`; all runs' mapping logs → `{sample}_mapping_log.txt` (summed into one bowtie2-style summary) |

`-a` reports every alignment for each read, so multimappers are kept at this
stage.

---

## 5. SAM filtering and BAM conversion

**Rule file:** `workflow/rules/sam_processing.smk`

| Order | Rule | Tool / script | Input → Output |
|------|------|---------------|----------------|
| 5a | `filter_sam` | `workflow/scripts/StructureFold3/sam_filter.py -max_mismatch <max_mismatch>` (env `samtools.yaml`; config `max_mismatch`, default 3) | Merged SAM → `{sample}_trimmed_mapped_filtered.sam` (+ filter log `.log` in the same dir). With its defaults this drops reverse-strand mappings, reads with a mismatch at the first base, and reads with more than `max_mismatch` mismatches/indels. |
| 5b | `filtered_sam_to_sorted_bam` | the merged SAM's header + the filtered SAM (`sam_filter.py` drops the header) piped into `samtools sort -@ 4`, then `samtools index` | → `{sample}_trimmed_mapped_filtered_sorted.bam` + `.bai` |

Downstream:

- The **filtered SAM** (5a) is what gets RT-stop counted (step 6).
- The **sorted BAM** (5b) is used for read-depth work: `samtools depth`
  coverage, the depth metagene, and the 3′ bias model.

---

## 6. RT-stop counting

**Rule file:** `workflow/rules/rtsc.smk`

| Rule | Script | Input → Output |
|------|--------|----------------|
| `sam_to_rtsc` | `workflow/scripts/StructureFold3/sam_to_rtsc.py -single` (env `rtsc_tools.yaml`) | Filtered SAM + `TRANSCRIPTOME_UPPER` → `{TMP}/output/se/{sample}/{sample}_trimmed_mapped_filtered.rtsc` |

A `.rtsc` file has, for each transcript, a name line and then one RT-stop count
per position: the number of reads whose 5′ end sits just downstream of that
nucleotide, i.e. where reverse transcription stopped.

---

## 7. Combining replicates / pooling

**Rule file:** `workflow/rules/rtsc.smk`

One rule, `combine_rtsc` (`workflow/scripts/StructureFold3/rtsc_combine.py`, env
`rtsc_tools.yaml`), sums per-sample `.rtsc` files into
`{TMP}/output/se/{id}/combined_{condition}.rtsc`. Which samples go in is decided by
`samples_to_combine()` in `helpers.smk`:

| `{id}` | Samples |
|--------|---------|
| a replicate `ID` | that ID's samples of the condition |
| `pooled` | every sample of the condition |
| `pooled_<temperature>` | every sample of the condition at that temperature |

Samples listed in `exclude_samples_plus` / `exclude_samples_minus` are left out here, but still
go through every per-sample step and per-sample QC.

Which combined files feed the reactivity calculation depends on
`pool_replicates` (see `get_minus_rtsc_for_reactivity` /
`get_raw_plus_rtsc` in `helpers.smk`):

| `pool_replicates` | Reactivity IDs | +DMS used | -DMS used |
|-------------------|----------------|-----------|-----------|
| `none` (default) | each replicate `ID` | that ID's own | that ID's own |
| `minus` | each replicate `ID` | that ID's own | `pooled` (all -DMS) |
| `plus` | each replicate `ID` | `pooled` (all +DMS) | that ID's own |
| `both` | `pooled` | `pooled` | `pooled` |
| `both_by_temperature` | `pooled_<t>` per temperature | `pooled_<t>` | `pooled_<t>` |
| `minus_all_plus_by_temperature` | `pooled_<t>` per temperature | `pooled_<t>` | `pooled` (all -DMS, every temperature) |

The replicate-comparison QC (correlation, UpSet across replicates) always
reads the per-replicate `combined_*.rtsc`, whatever `pool_replicates` is set
to.

---

## 8. Optional: 3′ coverage-bias correction

**Rule file:** `workflow/rules/bias_correction.smk`; included only when
`3_prime_bias_correction: true`. Background and model details are in
[3_prime_end_bias.md](3_prime_end_bias.md).

This step corrects the +DMS RT-stops for 3′ coverage bias (for example from
poly(A) selection), using -DMS read depth as the reference. It runs after
step 7 and before step 9.

| Order | Rule | Script (env `metagene.yaml`) | Input → Output |
|------|------|------------------------------|----------------|
| 8a | `bias_protein_coding_transcripts` | `three_prime_end_bias/list_protein_coding_transcripts.py` | `3_prime_bias_annotation` (falls back to `annotation_gtf`, then `gff`) + `TRANSCRIPTOME_UPPER` → `{TMP}/bias_correction/protein_coding_transcripts.txt` |
| 8b | `bias_training_transcripts` | `three_prime_end_bias/select_bias_training_transcripts.py` | Raw +DMS `.rtsc` of every ID → one training set shared by all IDs, from each protein-coding transcript's lowest +DMS coverage across IDs. `3_prime_bias_training: top_n` (default) keeps the top `3_prime_bias_n_transcripts` (default 2500); `all_covered` keeps every transcript with coverage ≥ `3_prime_bias_min_coverage` (default 1) in all IDs → `{out}/qc/3prime_bias_correction/training_transcripts.tsv` + `.bed` |
| 8c | `bias_training_depth` | `samtools depth -a -b <bed>` + `awk` | BAMs for one ID and channel (`plus`, or the `minus` reference: all -DMS samples if `3_prime_bias_reference: pooled`, otherwise the ID's own) → per-nt summed depth `{TMP}/bias_correction/{id}/depth_{condition}.tsv.gz` (temp) |
| 8d | `fit_3prime_bias_model` | `three_prime_end_bias/fit_3prime_bias_model.py` (uses `bias_model.py`) | Fits `log(+DMS depth / -DMS depth)` = natural cubic spline in distance from the 3′ end + relative position + natural cubic spline in distance from the 5′ end ([3_prime_end_bias.md §5](3_prime_end_bias.md#5-the-model)) → `{out}/{id}/3prime_bias_correction/model.tsv`, `coverage_residual_metagene.{tsv,png}` |
| 8e | `apply_3prime_bias_correction` | `three_prime_end_bias/apply_3prime_bias_correction.py` | Divides every transcript's +DMS RT-stops (mRNA, ncRNA and the positive control; the model is trained on protein-coding only) by the fitted bias, then rescales so the total matches the raw total. With correction on, every reactivity, p4p6 plot and ROC downstream uses these corrected stops → `{TMP}/bias_correction/{id}/combined_plus_corrected.rtsc` + `rtsc_metagene_before_after.{tsv,png}` |

When this module is on, `get_plus_rtsc_for_id()` hands
`combined_plus_corrected.rtsc` to every reactivity rule (steps 9–11) in place
of the raw `combined_plus.rtsc`. Every transcript is corrected, including ncRNAs and
the positive control, so the p4p6 plots and ROCs (steps 11d, 12.8) reflect the correction.

---

## 9. Reactivity calculation

**Rule file:** `workflow/rules/reactivity.smk`

| Order | Rule | Script | Input → Output |
|------|------|--------|----------------|
| 9a | `rtsc_to_react` | `workflow/scripts/StructureFold3/rtsc_to_react.py <minus> <plus> <fasta>` (env `rtsc_tools.yaml`) | +DMS `.rtsc` (raw or bias-corrected) + -DMS `.rtsc` (per §7) + `TRANSCRIPTOME_UPPER` → **`{out}/{id}/reactivity.react`** |
| 9b | `react_to_csv` | `workflow/scripts/reactivity/react_to_csv.py` (env `biopython.yaml`) | `reactivity.react` + FASTA → **`{out}/{id}/reactivity.csv`** (one row per transcript position, with the nucleotide). The same rule converts every other `.react` (raw, heat correction, norm comparison) |
| 9c | `react_to_rnaframework_xml` (`rnaframework_xml: true`) | `workflow/scripts/reactivity/react_to_rnaframework_xml.py` (env `biopython.yaml`) | `reactivity.react` + FASTA → `{out}/{id}/rnaframework_xml/<transcript>.xml`: rf-norm XML (attributes, order and 60-per-line wrapping as rf-norm 2.9.6; `scoring="Ding" pseudo="1" norm="2-8%" win=offset=1e9 reactive="AC" max="7"`; NA → NaN; transcripts with no value left out). The same script makes `depth_matched/{id}/rnaframework_xml/` (rule `depth_matched_to_rnaframework_xml`, from the CSV) and `{id}/heat_correction/rnaframework_xml_<suffix>/` (rule `heat_corrected_to_rnaframework_xml`) |

`rtsc_to_react.py` uses its defaults here: natural log of the stop counts,
-DMS subtracted from +DMS, A/C specificity (`-bases AC`), 2–8% normalization,
and reactivities capped at 7. The formula is in [methods.md §8](methods.md#8-reactivity).
With `reactivity_normalization: gu`, `rtsc_to_react`, `rtsc_to_raw_react` and both heat-correction
reactivity rules add `-gu_norm -min_gu {reactivity_min_gu}`: each channel is scaled by its mean over the
G/U positions instead of all positions before the −DMS subtraction
([methods.md §8, G/U normalization](methods.md#gu-normalization-reactivity_normalization-gu)).
The on-demand `rtsc_to_react_comparison` rule always uses the StructureFold2 scaling.

**On-demand rules** (not requested by `rule all`; ask for the output file by
name to build it):

| Rule | What it makes |
|------|---------------|
| `rtsc_to_raw_react` → `react_to_csv` | `{out}/{id}/raw_reactivity.{react,csv}`: no 2–8% normalization (`-nrm_off`) and effectively uncapped (`-threshold 1000000`) |
| `rtsc_to_react_comparison` → `react_to_csv` | `{out}/{id}/norm_comparison/{trim}_{nlog}_{norm}/reactivity.{react,csv}` for every combination of `notrim`/`trim125` (`-trim3 125`), `noln`/`ln` (`-ln_off`), `nonrm`/`nrm` (`-nrm_off`). Uses each replicate's unpooled `.rtsc`; trimmed 3′ positions are masked with `workflow/scripts/reactivity/mask_trim3_react.py` |

---

## 10. Optional: heat correction

**Rule file:** `workflow/rules/heat_correction.smk`; turned on by adding a
`temperature` column to the samplesheet (exactly 2 distinct values). It
follows Su et al. 2018 PNAS SI, "Determination of DMS reactivity", steps
3a/3b. All scripts are in `workflow/scripts/StructureFold3/`, env
`rtsc_tools.yaml` unless noted.

| Order | Rule | Script | Input → Output |
|------|------|--------|----------------|
| 10a | `heat_correction_plus_coverage` | `rtsc_coverage.py` | Pooled +DMS `.rtsc` at both temperatures → `{out}/qc/heat_correction_plus_coverage.csv` |
| 10b | `heat_correction_shared_transcripts` | `coverage_overlap.py -n 1.0` | Transcripts with coverage ≥ 1 at **both** temperatures → `{out}/qc/heat_correction_shared_transcripts.txt` |
| 10c | `heat_correction_lower_temperature_scale` | `rtsc_to_react.py -restrict` | Pooled lower-temperature +DMS and the −DMS background the reactivities use (all −DMS with `pool_replicates: minus` / `minus_all_plus_by_temperature`, else the pooled lower-temperature −DMS), restricted to the shared transcripts → `{out}/qc/heat_correction_lower_temperature.{react,scale}` (the 2–8% scale, step 3a) |
| 10d | `rtsc_to_react_heat_shared` | `rtsc_to_react.py -restrict -scale` | Each ID's reactivity, computed with the **lower-temperature scale** so the higher temperature doesn't renormalize away its own heat signal → `{out}/{id}/heat_correction/reactivity_precorrection.react` |
| 10e | `heat_correction_shared_react` | `react_intersect_transcripts.py` | Keeps only the transcripts present in every ID → `{out}/{id}/heat_correction/reactivity.react` |
| 10f | `react_to_csv` | `workflow/scripts/reactivity/react_to_csv.py` | → `{out}/{id}/heat_correction/reactivity.csv` |
| 10g | `heat_correct_reactivity` | `react_heat_correct.py -lower ... -higher ... -suffix` | Rescales lower- and higher-temperature reactivities to a common overall signal → `{out}/{id}/heat_correction/reactivity_{suffix}.react`; scale factors in `{out}/qc/heat_correction_scale_factors.log` |
| 10h | `react_to_csv` | `workflow/scripts/reactivity/react_to_csv.py` | → `{out}/{id}/heat_correction/reactivity_{suffix}.csv` |

`{suffix}` = `heat_correction_suffix` (default `corrected`).

**Comparison run** (`heat_correction_compare_all_transcripts: true`): the same
correction with no coverage threshold, starting from the plain
`{out}/{id}/reactivity.react`:

`heat_correction_all_transcripts_shared_react` (`react_intersect_transcripts.py`)
→ `heat_correct_reactivity_all_transcripts` (`react_heat_correct.py`)
→ `react_to_csv`.
The outputs go to `{out}/{id}/heat_correction_all_transcripts/` and
`{out}/qc/heat_correction_all_transcripts_scale_factors.log`.

---

## 11. Optional: depth-matched reactivity

**Rule file:** `workflow/rules/depth_matching.smk`; included only when
`depth_matching: true`, and requires at least 2 reactivity IDs. Details are
in [depth_matching.md](depth_matching.md).

| Order | Rule | Script | Input → Output |
|------|------|--------|----------------|
| 11a | `depth_matched_reactivity` | `workflow/scripts/depth_matching/depth_matched_reactivity.py` (env `metagene.yaml`, 8 threads, 32 GB) | Every ID's +DMS `.rtsc` (bias-corrected if §8 is on) and -DMS `.rtsc`. Each transcript's +DMS is thinned to the smallest total among the IDs (and the -DMS too if `depth_matching_match_minus`), reactivity is recomputed, and the mean over `depth_matching_draws` (default 100) draws is taken → `{out}/depth_matched/{id}/reactivity.csv`, `matched_plus_draw1.rtsc`, `{out}/depth_matched/summary.tsv` |
| 11b | `depth_matched_rtsc_correlation` | `StructureFold3/rtsc_correlation.py -spec AC` | First-draw matched +DMS → `{out}/qc/depth_matched/rt_stop_correlation.csv` |
| 11c | `plot_depth_matched_rt_stop_correlation` | `workflow/scripts/qc_plots/plot_rt_stop_correlation.py` | → `{out}/qc/depth_matched/rt_stop_replicate_correlation.png` |
| 11d′ | `plot_depth_matched_p4p6_all_samples` (p4p6 only) | `plot_p4p6_structure.py` | Depth-matched `reactivity.csv` → `{out}/qc/depth_matched/p4p6_all_samples.png` (P4-P6 structure, one panel per ID) |
| 11a′ | `depth_matched_reactivity_shard` + `depth_matched_reactivity` (`depth_matching_shards: N` > 1) | same script with `--shard k --n-shards N` (1 thread each), then `depth_matching/merge_depth_matched_shards.py` | N jobs, each one contiguous block of the sorted transcript list → `{tmp}/depth_matched_shards/{k}/` (temp); the merge concatenates them in shard order into the same outputs as 11a, byte-identical (draws are seeded per transcript index) |
| 11d | `plot_depth_matched_p4p6_roc` (p4p6 only) | `qc_plots/plot_p4p6_roc.py` | Depth-matched `reactivity.csv` → `{out}/qc/depth_matched/p4p6_roc.png` + `p4p6_auc.csv` (last `depth_matching_roc_trim3` nt left out) |

In genotype mode (§11G) the matched transcripts of the two genotypes share a key
(`depth_matching_key_map`), so each one is thinned to the same +DMS total in every library of
both genotypes; transcripts without a match in the other genotype get no depth-matched
reactivity. Steps 11b–11c are skipped (the genotypes' RT-stops sit under different names).

### 11G. Optional: two genotypes

**Rule file:** `workflow/rules/genotypes.smk`; included only when `genotypes` is set in the
config (exactly two entries, each with its own `genome` FASTA and annotation, given as `gtf` or
`gff`; both keys accept GTF and GFF3). Each genotype's transcripts are extracted from its own
genome the same way as `transcript_source: gffread` (§1), with `keep_longest_isoform` and
`cds_only` applied. The samplesheet needs a
`genotype` column. Heat correction is not supported, and `pool_replicates` must be `none` or
`minus`.

Each genotype's transcripts get the prefix `<genotype>__`, so the two genotypes share one
combined transcriptome downstream: steps 6–12 run unchanged on it. Reads map only to their own
genotype (step 2/4 use the genotype's own index), and each sample's `.rtsc` lists only its
genotype's transcripts. The positive control and `ncrna_fasta` are shared and keep their names.

| Order | Rule | Script | Input → Output |
|------|------|--------|----------------|
| G1 | `agat_longest_isoform_genotype` (if `keep_longest_isoform`), `gffread_extract_transcripts_genotype`, `filter_extracted_seqs_genotype` (if `cds_only`) | AGAT / gffread / `filter_extracted_seqs.py` | Each genotype's genome + GFF → `{TMP}/resources/genotypes/{genotype}/transcriptome_gffread[_filtered].fa` |
| G2 | `prefix_genotype_transcripts` | inline awk | → `{genotype}/transcripts.fa` (IDs `<genotype>__<id>`) |
| G3 | `combine_genotype_transcripts` | inline awk | Both → `all_genotypes_transcripts.fa`, the base of `TRANSCRIPTOME` (§1, plus the shared extras) |
| G4 | `genotype_mapping_fasta`, `uppercase_genotype_mapping_fasta`, `build_bowtie2_index_genotype` | awk / `uppercase_fasta.py` / bowtie2-build | Genotype transcripts + shared extras → `{genotype}/mapping.fa` (+ index, + uppercase copy used by `sam_to_rtsc`) |
| G5 | `match_genotype_transcripts` | `genotypes/match_transcripts.py` (env `genotypes.yaml`) | Matched pairs (identical IDs, or `transcript_matching`) and their position map: equal-length pairs ≥ 80% identical map 1:1, others are aligned (Biopython PairwiseAligner, global, end gaps free) → `{out}/genotype_comparison/matched_transcripts.tsv` (blocks as `start_a,start_b,length`, 0-based) |
| G6 | `pool_minus_across_genotypes` (if `pool_minus_across_genotypes`) | `genotypes/pool_minus_across_genotypes.py` | The -DMS pooled within each genotype (`pool_replicates: minus`) + the other genotype's stops at matching positions: from `cross_genotype_coordinate_map` when set (variant-derived; every pair sharing an ID), else the aligned blocks for pairs with identity ≥ `cross_genotype_min_identity` (0.9) → `{TMP}/output/se/pooled/combined_minus_cross_genotype.rtsc`, every reactivity's -DMS |
| G7 | `depth_matching_key_map` | inline awk | Pair map → `name, key` for §11 |
| G8 | `compare_genotypes` | `genotypes/compare_genotypes.py` (plotnine) | Final reactivity of every ID (depth-matched if §11 is on) → `{out}/genotype_comparison/reactivity_comparison.tsv` (per pair: n positions, Pearson, Spearman, mean abs Δ), `_positions.tsv.gz` (every shared A/C position), `_summary.tsv` and `.png` (per-transcript Pearson r: genotype means, one replicate each, and within each genotype as a baseline) |
| G9 | `rtsc_stop_correlation_genotype`, `plot_rt_stop_correlation_genotype`, `upset_covered_transcripts_genotype` | as §12.3–12.4 | Replicate QC within each genotype → `{out}/qc/rt_stop_replicate_correlation_{plus,minus}_{genotype}.png`, `covered_transcripts_upset_replicates_{genotype}.png` (replace the run-wide versions) |
| G10 | `bias_protein_coding_transcripts_genotype`/`_combined` (if §8), `extract_cds_and_ref_transcriptome_genotype`, `annotate_transcript_positions_genotype`, `combine_transcript_position_annotations` | as §8 / §12.7 | Each genotype's own annotation, renamed → `{TMP}/bias_correction/protein_coding_transcripts.txt`, `{out}/transcript_position_annotations.csv` |

3′ bias correction (§8) in genotype mode: the training set is chosen within each genotype (each
transcript scored over its own genotype's IDs; `top_n` keeps the top N of each genotype), and
with `3_prime_bias_reference: pooled` the -DMS reference depth comes from the ID's own
genotype's -DMS samples.

---

## 12. QC and summary outputs

These run beside the main path. They're grouped by the stage whose output
they read.

Every plot is set in Arial by `workflow/scripts/plot_style.py`, which each plotting script imports;
change the font there.

### 12.1 Read / alignment QC — `workflow/rules/qc.smk`

| Rule | Tool / script | Input → Output | In `rule all`? |
|------|---------------|----------------|----------------|
| `fastqc` | FastQC (env `qc.yaml`) | Trimmed FASTQ per (sample, run) → `{out}/qc/fastqc/{sample}_{run}_trimmed_fastqc.{html,zip}` | via `multiqc` only |
| `multiqc` | MultiQC | All FastQC zips → `{out}/qc/multiqc/multiqc_report.html` | No (on demand) |
| `multiqc_fastp` | MultiQC | All fastp JSONs → `{out}/qc/multiqc_fastp/multiqc_report.html` | Yes, when `trimmer: fastp` |
| `alignment_stats` | `awk` over the merged bowtie2 log | → `{out}/qc/alignment_stats/{sample}_stats.tsv` (total reads, unique %, multi %, overall %) | via aggregate |
| `aggregate_alignment_stats` | shell | → **`{out}/qc/alignment_stats_summary.tsv`** | Yes |
| `alignment_pct_positive_control` | `samtools view -F 0x904` + `awk` | Merged (unfiltered) SAM, primary alignments only (one per mapped read) → `{out}/qc/positive_control_alignment/{sample}_pct.tsv`: share of mapped reads on `positive_control_name` | via aggregate |
| `aggregate_positive_control_alignment` | shell | → `{out}/qc/positive_control_alignment_summary.tsv` | Yes, when `positive_control_fasta` is set |

### 12.2 RT-stop specificity — `workflow/rules/specificity.smk`

| Rule | Script | Input → Output |
|------|--------|----------------|
| `calculate_specificity_by_sample` | `StructureFold3/rtsc_specificity.py` | Per-sample `.rtsc` → `{out}/specificity_{sample}.csv` (base composition of the RT-stops) |
| `calculate_specificity` | `StructureFold3/rtsc_specificity.py` | Per-ID `combined_{plus,minus}.rtsc` → `{out}/{id}/specificity_{treatment}.csv` |
| `plot_specificity` | `qc_plots/plot_specificity.py` (plotnine) | All per-sample CSVs → `{out}/qc/specificity_plot.png` |
| `plot_specificity_comparison` | `qc_plots/plot_specificity_comparison.py` (plotnine) | Per-sample specificity CSVs of every run in `specificity_comparison_runs` → `{out}/qc/specificity_comparison.png` (one panel per run) + `specificity_comparison_manifest.tsv` |

For DMS, a good library has +DMS stops enriched at A and C, and -DMS stops
without that enrichment.

### 12.3 Coverage — `workflow/rules/coverage.smk`

| Rule | Script | Input → Output | In `rule all`? |
|------|--------|----------------|----------------|
| `calculate_stop_coverage` | `StructureFold3/rtsc_coverage.py` | Per-ID `combined_plus.rtsc` → **`{out}/{id}/coverage.csv`** (per-transcript RT-stop coverage) | Yes |
| `upset_covered_transcripts_replicates` | `plot_covered_transcripts_upset.py` (env `upset.yaml`) | `coverage.csv` of each true replicate (threshold 1) → `{out}/qc/covered_transcripts_upset_replicates.png` | Yes, if ≥ 2 replicates |
| `upset_covered_transcripts_merged` | `plot_covered_transcripts_upset.py` | `coverage.csv` of each pooled ID → `{out}/qc/covered_transcripts_upset_merged.png` | Yes, if pooling changed `IDS` and left ≥ 2 |
| `samtools_depth_by_condition` | `samtools depth -a` | BAMs of one ID/condition → `{out}/{id}/depth_{condition}.txt` | No (on demand) |
| `calculate_nucleotide_coverage` | `rtsc_to_nucleotide_coverage.py` | `combined_{condition}.rtsc` → `{out}/{id}/nucleotide_coverage_{condition}.csv` | No (on demand) |
| `filter_covered_transcripts` | `filter_covered_transcripts.py` | All `coverage.csv` → `{out}/qc/covered_transcripts.txt` (coverage ≥ 1 in every ID) | No (on demand) |

### 12.4 Replicate correlation — `workflow/rules/replicate_correlation.smk`

| Rule | Script | Input → Output |
|------|--------|----------------|
| `rtsc_stop_correlation` | `StructureFold3/rtsc_correlation.py -spec AC` | Each true replicate's `combined_plus.rtsc` → `{out}/qc/rt_stop_correlation.csv` |
| `plot_rt_stop_correlation` | `plot_rt_stop_correlation.py` | → `{out}/qc/rt_stop_replicate_correlation.png` (only if ≥ 2 +DMS replicates) |
| `rtsc_stop_correlation_minus` | `StructureFold3/rtsc_correlation.py -spec AC` | Each true replicate's `combined_minus.rtsc` → `{out}/qc/rt_stop_correlation_minus.csv` |
| `plot_rt_stop_correlation_minus` | `plot_rt_stop_correlation.py` | → `{out}/qc/rt_stop_replicate_correlation_minus.png` (only if ≥ 2 -DMS replicates) |

### 12.5 Metagene coverage — `workflow/rules/metagene.smk`

For each transcript, coverage `c` becomes `ln(c+1) / mean(ln(c+1))`, which is
then binned by relative position (`metagene_bins`, default 100) and averaged
across transcripts. The profile rules use env `metagene.yaml`; the plot rules
use `plotnine.yaml`.

| Rule | Script | Input → Output |
|------|--------|----------------|
| `metagene_profile_depth` | `samtools depth -a` piped into `metagene_profile.py --mode depth` | Sorted BAM → `{out}/qc/metagene/profiles/{sample}_depth.tsv` |
| `metagene_profile_rtsc` | `metagene_profile.py --mode rtsc` | Per-sample `.rtsc` → `{out}/qc/metagene/profiles/{sample}_rtsc.tsv` |
| `plot_metagene_depth` | `qc_plots/plot_metagene.py` (plotnine, env `plotnine.yaml`) | All depth profiles → `{out}/qc/metagene/metagene_depth.png` |
| `plot_metagene_rtsc` | `qc_plots/plot_metagene.py` (plotnine, env `plotnine.yaml`) | All RT-stop profiles → `{out}/qc/metagene/metagene_rtsc.png` |

### 12.6 Abundance from -DMS — `workflow/rules/abundance.smk`

| Rule | Script | Input → Output | In `rule all`? |
|------|--------|----------------|----------------|
| `rtsc_total_stops` | `rtsc_total_stops.py` | `combined_minus.rtsc` → `{out}/{id}/counts_minus.csv` (total -DMS stops per transcript) | Yes |
| `calculate_transcript_abundance` | `StructureFold3/rtsc_abundances.py {RPKM\|TPM}` | `combined_minus.rtsc` → `{out}/{id}/abundance_{mode}.csv` | Yes, for each mode in `transcript_abundance` |

### 12.7 Transcript position annotation — `workflow/rules/annotation.smk`

Defined only when both `annotation_gtf` and `genome` are set.

| Order | Rule | Tool / script | Input → Output |
|------|------|---------------|----------------|
| a | `extract_cds_and_ref_transcriptome` | `gffread -x -w` (env `gffread.yaml`) | `genome` + `annotation_gtf` → `{TMP}/resources/cds.fa`, `ref_transcriptome.fa` |
| b | `annotate_transcript_positions` | `annotation/annotate_transcript_positions.py` (env `biopython.yaml`) | Finds the CDS inside each reference transcript and maps it onto `TRANSCRIPTOME` → `{out}/transcript_position_annotations.csv` (5UTR / CDS / 3UTR / noncoding, plus codon position for CDS bases) |

### 12.8 Positive control (p4p6) — `workflow/rules/positive_control.smk`

Included only when `positive_control_name: p4p6`. Draws reactivity onto the
Tetrahymena P4-P6 secondary structure (`resources/p4p6/p4p6.xrna`) with
`plot_p4p6_structure.py` (matplotlib, env `plotnine.yaml`): residues as circles coloured by reactivity bin (black < 0.3, gold 0.3–0.7, brick red ≥ 0.7, grey = no data) with contrasting letters, a backbone line, G·C / A·U / G·U pair symbols, and numbering every 10 nt (transcript positions). The construct has 2 extra 5′ nt, so
reactivity position `p` is drawn at XRNA nucleotide `p − 2`.

| Rule | Script | Input → Output |
|------|--------|----------------|
| `plot_p4p6_react` | `plot_p4p6_structure.py` | One ID's `reactivity.csv` → `{out}/{id}/p4p6_react.png` |
| `plot_p4p6_all_samples` | `plot_p4p6_structure.py` | Every ID's `reactivity.csv` → `{out}/qc/p4p6_all_samples.png` |
| `plot_p4p6_roc` | `qc_plots/plot_p4p6_roc.py` (plotnine) | Every ID's `reactivity.csv` → `{out}/qc/p4p6_roc.png` + `p4p6_auc.csv`: ROC over A/C, positive = unpaired, negative = Watson–Crick–Franklin-paired (other pairs left out); last `p4p6_roc_trim3` nt (default 50) left out |

### 12.9 Read composition and the HTML report — `workflow/rules/report.smk`

| Rule | Tool / script | Input → Output |
|------|---------------|----------------|
| `transcript_categories` | `qc_plots/transcript_categories.py` | Base transcriptome → `mRNA` (records with an Ensembl `transcript_biotype` other than `protein_coding` → `ncRNA`), rRNA biotypes / rRNA set / rRNAs the ncRNA annotation marks → `rRNA`, other ncRNA set → `ncRNA`, lncRNA set → `lncRNA` (its own type only when set), `positive_control_fasta` → `positive_control` → `{TMP}/resources/transcript_categories.tsv` |
| `read_composition` | `samtools view` + `workflow/scripts/qc_plots/read_composition.awk` | Merged, unfiltered SAM (unmapped reads included) → `{out}/qc/read_composition/{sample}.tsv`. Only forward-strand alignments give a type, as in the reactivity: **unique** = reads with exactly one alignment, by type; **by type** = reads whose forward alignments all hit one type (keeps reads shared between isoforms / rRNA copies); **type-ambiguous** = reads hitting more than one type; **antisense-only** = mapped reads with no forward alignment; **pie_*** = every read in one slice (first `read_composition_priority` type hit, else other; unmapped) |
| `aggregate_read_composition` | inline pandas | → `{out}/qc/read_composition_summary.tsv` (counts and %, one row per sample) |
| `plot_read_composition_pie` | `workflow/scripts/qc_plots/plot_read_composition_pie.py` (plotnine) | → `{out}/qc/read_composition_pie.png` (one pie per sample + all samples: mRNA / rRNA / lncRNA / other / unmapped) and `read_composition_pie_mapped.png` (the same, mapped reads only), and `.tsv` |
| `qc_report` | `workflow/scripts/qc_plots/make_report.py` (env `plotnine.yaml`) | Every QC table/figure above → **`{out}/qc/report.html`**: one self-contained file (figures embedded as base64 PNG) with run settings and samples, mapping, read composition, positive control, specificity, coverage, replicate correlation, metagene, and — when enabled — 3′ bias correction, depth matching (per-ID thinning summary, matched RT-stop correlation). The positive-control section shows the p4p6 structure figure, AUC table and ROC for the final reactivities only — the depth-matched ones when depth matching is on (built from 3′ bias-corrected stops when that is on too). Heat correction and, in genotype mode, a genotype-comparison section (matched pairs, per-transcript agreement table and plot) are added when enabled. Links to the fastp MultiQC report when `trimmer: fastp`. Title from `report_title` |

### 12.10 Optional: transcriptome-unaligned reads on the genome — `workflow/rules/genome_check.smk`

Included when `genome_check: true` (needs `genome`; in genotype mode each sample uses its own
genotype's genome). The reads Bowtie2 left unaligned (flag 4) are aligned to the genome with HISAT2
(splice-aware; `--no-softclip` so a hit is end-to-end like the transcriptome alignment; up to 10
loci per read) and placed by their 5′ end (the base after the RT stop; reads are sense to the RNA)
on `genome_check_annotation` (default `annotation_gtf`, then `gff`). Categories, first that applies:
`organelle (<chrom>)` (chromosome matches `genome_check_organelles`), `transcriptome exon` (exon of
a transcript the reads were aligned to: usually an unannotated junction, indels, or a base-1
mismatch), `<biotype> exon` (an annotated transcript left out of the transcriptome: other isoforms,
ncRNAs dropped by `cds_only`, rRNA, tRNA …), `intron`, `antisense` (inside a gene on the other
strand only), `intergenic`.

| Rule | Tool / script | Input → Output |
|------|---------------|----------------|
| `genome_check_hisat2_index` | `hisat2-build` (env `hisat2.yaml`, 8 threads, 16 GB) | Genome → `{TMP}/genome_check/index/{ref}/genome.*.ht2` (skipped when `genome_check_hisat2_index` / a genotype's `hisat2_index` is given) |
| `genome_check_unaligned_reads` | `samtools fastq -f 4` | Merged transcriptome SAM → `{TMP}/genome_check/{sample}/unaligned.fastq.gz` (temp) |
| `genome_check_align` | `hisat2` + `samtools` (env `hisat2.yaml`) | → `{out}/genome_check/{sample}.bam(.bai)`, HISAT2 summary, primary alignments table (temp) |
| `genome_check_annotation` | `genome_check/annotation_intervals.py` (plotnine env) | GTF/GFF3 + `TRANSCRIPTOME` → exon and gene intervals flagged by whether the transcript is in the transcriptome |
| `genome_check_classify` | `genome_check/classify_genome_hits.py` | One sample's primary alignments → reads per (chrom, strand, 5′ end, category, gene) and per category |
| `genome_check_summary` | `genome_check/summarize_genome_check.py` | Every sample → `{out}/genome_check/summary.tsv` (unaligned reads, % on genome, % per category), `categories.png`, `top_loci.tsv` (top `genome_check_top_loci` loci, 5′ ends ≤ `genome_check_locus_gap` nt apart, per-sample counts and main gene); also a report section |

---

## 13. Final outputs (`rule all`)

Built by `get_all_targets()` in `workflow/rules/targets.smk`.

**Always:**

| Output | Produced by |
|--------|-------------|
| `{out}/{id}/reactivity.react` | `rtsc_to_react` |
| `{out}/{id}/reactivity.csv` | `react_to_csv` |
| `{out}/{id}/coverage.csv` | `calculate_stop_coverage` |
| `{out}/{id}/specificity_plus.csv`, `specificity_minus.csv` | `calculate_specificity` |
| `{out}/{id}/counts_minus.csv` | `rtsc_total_stops` |
| `{out}/specificity_{sample}.csv` | `calculate_specificity_by_sample` |
| `{out}/qc/alignment_stats_summary.tsv` | `aggregate_alignment_stats` |
| `{out}/qc/specificity_plot.png` | `plot_specificity` |
| `{out}/qc/metagene/metagene_depth.png`, `metagene_rtsc.png` | `plot_metagene_depth`, `plot_metagene_rtsc` |
| `{out}/qc/read_composition_summary.tsv` | `aggregate_read_composition` |
| `{out}/qc/read_composition_pie.png`, `_mapped.png`, `.tsv` | `plot_read_composition_pie` |
| `{out}/qc/report.html` | `qc_report` |

The per-ID QC files (`coverage`, `specificity_*`, `counts_minus`, `abundance_*`)
are only requested for IDs whose own channel still has at least one sample
after exclusions (`IDS_WITH_PLUS` / `IDS_WITH_MINUS`).

**Conditional:**

| Condition | Outputs |
|-----------|---------|
| ≥ 2 +DMS replicates | `qc/rt_stop_replicate_correlation.png`, `qc/covered_transcripts_upset_replicates.png` |
| ≥ 2 -DMS replicates | `qc/rt_stop_replicate_correlation_minus.png` |
| Pooling changed `IDS` and left ≥ 2 | `qc/covered_transcripts_upset_merged.png` |
| `annotation_gtf` + `genome` | `transcript_position_annotations.csv` |
| `transcript_abundance` | `{id}/abundance_{RPKM,TPM}.csv` |
| `positive_control_fasta` | `qc/positive_control_alignment_summary.tsv` |
| `positive_control_name: p4p6` | `{id}/p4p6_react.png`, `qc/p4p6_all_samples.png`, `qc/p4p6_roc.png`, `qc/p4p6_auc.csv` |
| `specificity_comparison_runs` | `qc/specificity_comparison.png` |
| `trimmer: fastp` | `qc/multiqc_fastp/multiqc_report.html` |
| `3_prime_bias_correction` | `{id}/3prime_bias_correction/{model.tsv, coverage_residual_metagene.png, rtsc_metagene_before_after.png}` |
| `depth_matching` | `depth_matched/{id}/reactivity.csv`, `qc/depth_matched/rt_stop_replicate_correlation.png` (+ `qc/depth_matched/p4p6_roc.png` and `p4p6_all_samples.png` with p4p6) |
| `genotypes` | `genotype_comparison/matched_transcripts.tsv`, `genotype_comparison/reactivity_comparison.png` (+ its TSVs), `transcript_position_annotations.csv`, per-genotype `qc/rt_stop_replicate_correlation_{plus,minus}_{genotype}.png` and `qc/covered_transcripts_upset_replicates_{genotype}.png` (in place of the run-wide ones) |
| `genome_check` | `genome_check/{summary.tsv, categories.png, top_loci.tsv}`, `genome_check/{sample}.bam(.bai)` |
| `temperature` column | `{id}/heat_correction/reactivity.csv`, `reactivity_{suffix}.{react,csv}` |
| `rnaframework_xml` | `{id}/rnaframework_xml/` (+ `depth_matched/{id}/rnaframework_xml/` with depth matching, `{id}/heat_correction/rnaframework_xml_{suffix}/` with heat correction) |
| …plus `heat_correction_compare_all_transcripts` | `{id}/heat_correction_all_transcripts/reactivity_{suffix}.{react,csv}` |

---

## 14. Script index

Paths are relative to `workflow/scripts/`.

| Script | Used by rule(s) | Step |
|--------|-----------------|------|
| `transcriptome/filter_extracted_seqs.py` | `filter_extracted_seqs` | 1 |
| `transcriptome/uppercase_fasta.py` | `uppercase_transcriptome` | 1 |
| `StructureFold3/fastq_trimmer.py` | `trim_reads_cutadapt` | 3 |
| `alignment/combine_mapping_logs.py` | `combine_sam_from_multiple_runs` | 4 |
| `StructureFold3/sam_filter.py` | `filter_sam` | 5 |
| `StructureFold3/sam_to_rtsc.py` | `sam_to_rtsc` | 6 |
| `StructureFold3/rtsc_combine.py` | `combine_rtsc` | 7 |
| `three_prime_end_bias/list_protein_coding_transcripts.py` | `bias_protein_coding_transcripts` | 8 |
| `three_prime_end_bias/select_bias_training_transcripts.py` | `bias_training_transcripts` | 8 |
| `three_prime_end_bias/fit_3prime_bias_model.py`, `three_prime_end_bias/bias_model.py` | `fit_3prime_bias_model` | 8 |
| `three_prime_end_bias/apply_3prime_bias_correction.py` | `apply_3prime_bias_correction` | 8 |
| `three_prime_end_bias/plot_bias.py` | `fit_3prime_bias_model`, `apply_3prime_bias_correction` (plotnine figures, imported); run directly to redraw an ID's figures from its saved `3prime_bias_correction/` TSVs | 8 |
| `StructureFold3/rtsc_to_react.py` | `rtsc_to_react`, `rtsc_to_raw_react`, `rtsc_to_react_comparison`, `heat_correction_lower_temperature_scale`, `rtsc_to_react_heat_shared` | 9, 10 |
| `reactivity/react_to_csv.py` | `react_to_csv` | 9, 10 |
| `reactivity/mask_trim3_react.py` | `rtsc_to_react_comparison` | 9 |
| `StructureFold3/rtsc_coverage.py` | `calculate_stop_coverage`, `heat_correction_plus_coverage` | 10, 12.3 |
| `StructureFold3/coverage_overlap.py` | `heat_correction_shared_transcripts` | 10 |
| `StructureFold3/react_intersect_transcripts.py` | `heat_correction_shared_react`, `heat_correction_all_transcripts_shared_react` | 10 |
| `StructureFold3/react_heat_correct.py` | `heat_correct_reactivity`, `heat_correct_reactivity_all_transcripts` | 10 |
| `depth_matching/depth_matched_reactivity.py` | `depth_matched_reactivity` | 11 |
| `genotypes/match_transcripts.py` | `match_genotype_transcripts` | 11G |
| `genotypes/pool_minus_across_genotypes.py` | `pool_minus_across_genotypes` | 11G |
| `genotypes/compare_genotypes.py` | `compare_genotypes` | 11G |
| `StructureFold3/rtsc_correlation.py` | `rtsc_stop_correlation`, `rtsc_stop_correlation_minus`, `depth_matched_rtsc_correlation` | 11, 12.4 |
| `qc_plots/plot_rt_stop_correlation.py` | `plot_rt_stop_correlation`, `plot_rt_stop_correlation_minus`, `plot_depth_matched_rt_stop_correlation` | 11, 12.4 |
| `StructureFold3/rtsc_specificity.py` | `calculate_specificity`, `calculate_specificity_by_sample` | 12.2 |
| `qc_plots/plot_specificity.py` | `plot_specificity` | 12.2 |
| `qc_plots/plot_specificity_comparison.py` | `plot_specificity_comparison` | 12.2 |
| `qc_plots/plot_covered_transcripts_upset.py` | `upset_covered_transcripts_*` | 12.3 |
| `coverage/rtsc_to_nucleotide_coverage.py` | `calculate_nucleotide_coverage` | 12.3 |
| `coverage/filter_covered_transcripts.py` | `filter_covered_transcripts` | 12.3 |
| `coverage/metagene_profile.py` | `metagene_profile_depth`, `metagene_profile_rtsc` | 12.5 |
| `qc_plots/plot_metagene.py` | `plot_metagene_depth`, `plot_metagene_rtsc` | 12.5 |
| `coverage/rtsc_total_stops.py` | `rtsc_total_stops` | 12.6 |
| `StructureFold3/rtsc_abundances.py` | `calculate_transcript_abundance` | 12.6 |
| `annotation/annotate_transcript_positions.py` | `annotate_transcript_positions` | 12.7 |
| `qc_plots/plot_p4p6_structure.py` | `plot_p4p6_react`, `plot_p4p6_all_samples` | 12.8 |
| `genome_check/annotation_intervals.py` | `genome_check_annotation` | 12.10 |
| `genome_check/classify_genome_hits.py` | `genome_check_classify` | 12.10 |
| `genome_check/summarize_genome_check.py` | `genome_check_summary` | 12.10 |

Folders group scripts by pipeline stage: `transcriptome/`, `alignment/`, `reactivity/`,
`three_prime_end_bias/`, `depth_matching/`, `genotypes/`, `genome_check/`, `coverage/`, `annotation/` and
`qc_plots/` (every plotting and report script). `StructureFold3/` holds the Python 3 ports of
the StructureFold2 tools.

Not used by any rule: `regions/` (`extract_region_fasta.py`, `extract_region_rtsc.py`,
`region_coordinates.py`; standalone helpers for slicing transcripts to a region),
`deprecated/` (superseded scripts: `rtsc_to_react_plus_only.py`, the old matplotlib
`plot_p4p6_roc.py`, `plot_p4p6_rnavigate.py`, and an older copy of
`react_intersect_transcripts.py`), `StructureFold3/fastq_mapper.py`, and the original Python 2
`StructureFold2/` scripts (kept for reference).

---

## 15. Known issues

None open. (`plot_specificity_comparison` used to call the deleted `plot_specificity_combined.R`; it
now runs `qc_plots/plot_specificity_comparison.py`, a plotnine port.)
