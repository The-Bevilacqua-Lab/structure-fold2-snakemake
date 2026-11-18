## Samplesheet

Tab-separated, one row per (sample, sequencing run). Examples: `config/samplesheet_minimal_example.tsv`
(the three required columns) and `config/samplesheet_example.tsv`.

| column | required | meaning |
|---|---|---|
| `sample` | yes | Library name, used in file names. One library = one `sample`. |
| `condition` | yes | `plus` (+DMS) or `minus` (−DMS) |
| `r1` | yes | That (sample, run)'s FASTQ, plain or gzipped |
| `ID` | no | Replicate: the +DMS and −DMS samples with the same `ID` give one reactivity (e.g. `rep1`). **Default:** the n-th `minus` sample and the n-th `plus` sample, in sheet order, form `rep<n>`. `pooled` and `pooled_*` are reserved. |
| `run` | no | Sequencing run, for a library sequenced more than once (e.g. `run1`, `run2`). Each run is trimmed and aligned on its own, then the runs are merged into one SAM per sample. **Default:** one run. |
| `temperature` | no | Probing temperature, exactly two distinct values (e.g. `23` and `42`, or `ambient` and `heat`). Turns on heat correction between them ([pipeline.md §10](pipeline.md#10-optional-heat-correction)). Every row of an `ID` must have the same temperature. |
| `genotype` | with `genotypes` in the config | Which of the two `genotypes` the row belongs to. Every row of an `ID` must have the same genotype, and each genotype needs `plus` and `minus` samples. Not combinable with `temperature`. See [pipeline.md §11G](pipeline.md#11g-optional-two-genotypes). |

The pipeline stops with an error naming the problem if a required column is missing, a
`condition` is not `plus`/`minus`, a (sample, run) pair appears twice, or one sample's rows
disagree on `condition`, `ID` or `temperature`.

Example with two replicates, one library sequenced twice:

```
sample	ID	condition	run	r1
S1	rep1	minus	run1	reads/S1_run1.fastq.gz
S2	rep2	minus	run1	reads/S2_run1.fastq.gz
S3	rep1	plus	run1	reads/S3_run1.fastq.gz
S3	rep1	plus	run2	reads/S3_run2.fastq.gz
S4	rep2	plus	run1	reads/S4_run1.fastq.gz
```
