# Methods

How the pipeline turns reads into DMS reactivities, step by step and with the actual
calculations. [pipeline.md](pipeline.md) lists the rules and files for each step; this page
explains what they compute and why. The optional corrections have their own pages:

- [3_prime_end_bias.md](3_prime_end_bias.md): 3′ coverage-bias correction of the +DMS RT stops
- [depth_matching.md](depth_matching.md): depth-matched reactivity
- Heat correction: [pipeline.md §10](pipeline.md#10-optional-heat-correction) (Su et al. 2018)

A methods paragraph you can adapt for a paper is at the end.

## Contents

1. [What an RT stop is](#1-what-an-rt-stop-is)
2. [Reads and trimming](#2-reads-and-trimming)
3. [Reference transcriptome](#3-reference-transcriptome)
4. [Alignment](#4-alignment)
5. [Read filtering](#5-read-filtering)
6. [RT-stop counting](#6-rt-stop-counting)
7. [Combining samples](#7-combining-samples)
8. [Reactivity](#8-reactivity)
9. [Coverage](#9-coverage)
10. [Methods paragraph](#10-methods-paragraph)

---

## 1. What an RT stop is

In Structure-seq, DMS methylates unpaired A (N1) and C (N3). Reverse transcriptase stops one
nucleotide before a methylated base, so the cDNA's 3′ end, and the 5′ end of a read in the sense
orientation, lies on the nucleotide just 3′ of the modification:

```
RNA      5′ ─ ─ ─ ─ ─ A* x ─ ─ ─ ─ ─ ─ 3′        A* = methylated base at position i
read                     x ─ ─ ─ ─ ─ ─ ─        read 5′ end (SAM POS) = i + 1
```

So a read whose alignment starts at transcript position `POS` (1-based) is an RT stop at position
`POS − 1`. Every count in the pipeline is of read 5′ ends. Coverage across a read's body only
appears in the read-depth steps (metagene depth plots and the 3′ bias model).

The −DMS library is made the same way without DMS. Its stops come from RT pausing at structure,
natural modifications and RNA ends, and it's subtracted as background (§8).

## 2. Reads and trimming

The pipeline takes **single-end** reads (one FASTQ per sample and sequencing run). Reads must be
in the **sense orientation of the RNA**: §5 drops every reverse-strand alignment, so a library
whose reads are antisense would lose almost all of them.

Trimming, per (sample, run):

- **cutadapt** (default; StructureFold2's `fastq_trimmer.py`): three chained cutadapt steps. They
  trim the 5′ adapter `TGAACAGCGACTAGGCTCTTCA`, then the 3′ adapter `GATCGGAAGAGCACACGTCTG`, then
  quality-trim at Q30 and drop reads shorter than 20 nt. These adapters are the Structure-seq2
  defaults; other library preps need their own adapters (see the README).
- **fastp**: trims `fastp_adapter` (default `GATCGGAAGAGCACACGTCTG`), poly-G tails, and fastp's
  default quality filter, and drops reads shorter than `fastp_min_length` (default 20). fastp
  does not trim a 5′ adapter.

Several runs of one sample are trimmed and aligned separately and merged after alignment.

## 3. Reference transcriptome

Reads are aligned to transcripts, not the genome. The reference is either a FASTA you provide
(`transcript_source: file`) or built from a genome and GFF/GTF with gffread
(`transcript_source: gffread`). The gffread build can first keep only the longest isoform per gene
(`keep_longest_isoform`, using AGAT) and by default keeps only transcripts with a CDS (`cds_only`).

`positive_control_fasta` (e.g. the Tetrahymena P4-P6 domain) and `ncrna_fasta` (rRNAs, tRNAs …)
are appended. Every read competes for the same reference, so adding rRNA keeps rRNA reads from
landing on similar mRNA sequence. An uppercase copy of the FASTA is used for every base check
afterwards, so soft-masked (lowercase) bases still count as A or C.

## 4. Alignment

Bowtie2 aligns each run's trimmed reads with `bowtie2 -a`: default `--sensitive` end-to-end mode
(no soft clipping, so the read's 5′ end is the true end) and **every** valid alignment reported,
not just the best. A read that matches several isoforms of one gene aligns to all of them.

## 5. Read filtering

StructureFold2's `sam_filter.py` (Python 3 port), with these settings:

| Removed | Why |
|---|---|
| Unaligned reads (flag 4) | — |
| Reverse-strand alignments (flag 16) | Sense reads only align forward to the transcript; reverse hits are spurious |
| Alignments with more than `max_mismatch` mismatches + indels (`NM` tag; default 3) | Poor alignments |
| Alignments with a mismatch at the read's first base (`MD:Z:0…`) | RT often adds a non-templated base at the cDNA end, so a first-base mismatch puts the stop in the wrong place |

Secondary alignments are **kept**: a read that aligns to several transcripts counts as an RT stop
on each of them. This is StructureFold2's behaviour, and it is what lets isoforms that share
sequence all get reactivity. Read depth from the sorted BAM (`samtools depth`) counts only
primary alignments.

## 6. RT-stop counting

`sam_to_rtsc.py` turns each filtered SAM into an `.rtsc` file: for every transcript, a vector of
length `L` (the transcript length) whose 0-based entry `i` counts alignments starting at 1-based
position `i + 1`. By §1, entry `i` is therefore the number of stops **at base `i`** (1-based).
Entry 0 (reads starting at base 1) has no base to stop on and is ignored downstream. The last base,
`L`, can't be measured because a stop there would need a read starting past the end.

## 7. Combining samples

Per reactivity ID, the per-sample `.rtsc` vectors are summed position by position: +DMS samples
into `combined_plus.rtsc`, −DMS samples into `combined_minus.rtsc`. Which samples are summed
depends on the samplesheet `ID` column and `pool_replicates`
([pipeline.md §7](pipeline.md#7-combining-replicates--pooling)). For example, `minus` sums every
−DMS sample into one shared background. `exclude_samples_plus` / `exclude_samples_minus` leave
samples out of these sums only.

## 8. Reactivity

StructureFold2's `rtsc_to_react.py` (Python 3 port, same numbers as the original), with its
defaults: natural log on, A/C specificity, 2–8% normalization, cap 7. For one transcript with
+DMS stops `p_i` and −DMS stops `m_i` (combined, §7) over entries `i = 0 … L−1`:

1. **Log.** `P_i = ln(p_i + 1)`, `M_i = ln(m_i + 1)`. The log dampens very high-count positions,
   for example PCR duplicates or strong RT pauses.
2. **Scale each channel to its mean.** `P_i · L / ΣP` and `M_i · L / ΣM` (sums over the whole
   transcript), so each channel averages 1 per position and library depth cancels out. A
   transcript without stops in either channel gets no reactivity.
3. **Subtract the background.**

   ```
   θ_i = max(0,  L·P_i / ΣP  −  L·M_i / ΣM)
   ```

   Negative values become 0. Only positions `1 … L−1` are kept (§6), and only those within both
   vectors' length.
4. **2–8% normalization.** Take θ at the A/C positions in `1 … L−1` (the counted set includes
   zeros) and sort high to low. With `n` of them, average the values ranked from `⌊0.02 n⌋` to
   `⌊0.10 n⌋` (0-based, end exclusive): the top 2% are skipped as outliers and the next 8% are
   averaged. That average, the transcript's normalizer, is written to the `.scale` file. A
   transcript whose normalizer is 0 is dropped. With fewer than 10 A/C positions the window is
   empty and the transcript is always dropped.
5. **Final value.** `reactivity_i = min(θ_i / normalizer, 7)`, rounded to 3 decimals. G and U
   positions are `NA`, as is the last base.

A reactivity of 1 is therefore "as reactive as a typical highly reactive base in this transcript".
Values are comparable between transcripts because each one is normalized to itself.

### G/U normalization (`reactivity_normalization: gu`)

Step 2 divides each channel by its mean over **all** positions. For −DMS that mean is background.
For +DMS it is background plus the DMS signal at A/C, so `L·P_i / ΣP` puts the +DMS background
*below* the −DMS background, and step 3 subtracts too much. The over-subtraction grows with the
transcript's (and the library's) DMS signal, so it differs between libraries and conditions: it
is a technical difference that can look like a structural one.

With `reactivity_normalization: gu`, step 2 instead divides each channel by its mean over the
transcript's G and U positions, which DMS doesn't modify on the Watson-Crick face:

```
θ_i = max(0,  P_i / mean_{G/U} P  −  M_i / mean_{G/U} M)
```

where the means are over the positions `1 … L−1` (within both vectors) whose base is G or U. Both
channels are then on the same background scale, and only the excess +DMS signal remains. Steps
1, 3 (clipping at 0), 4 and 5 are unchanged, so the result is still 2–8% normalized over A/C and
capped at 7. A transcript with fewer than `reactivity_min_gu` (default 20) G/U positions, or no
stops at G/U in either channel, gets no reactivity.

Equivalently, the default subtracts the background with weight
`mean_all(P) / mean_all(M)` and G/U normalization with `mean_GU(P) / mean_GU(M)`. On the rice
in vitro Structure-seq libraries, the ratio of the two (the G/U slope) was 0.71–0.87: the default
removed 1.15–1.4× too much background, more at 42 °C. G/U normalization raised replicate agreement
in every library pair and NB–IR64 agreement on identical-sequence transcripts, and left the
p4p6 and 18S rRNA AUCs unchanged. It applies to the standard, raw and heat-corrected reactivity
(`rtsc_to_react.py -gu_norm`) and to depth matching (`depth_matched_reactivity.py --normalization gu`;
[depth_matching.md](depth_matching.md)).

Heat correction swaps step 4's normalizer for the lower temperature's at both temperatures, then
rescales one temperature to the other
([pipeline.md §10](pipeline.md#10-optional-heat-correction)). The optional 3′ bias correction
changes `p_i` before step 1, and depth matching thins `p_i` before step 1. Both are described on
their own pages.

## 9. Coverage

StructureFold's coverage (in `coverage.csv`, and used to pick the 3′ bias training transcripts):

```
coverage = (RT stops on A and C bases) / (number of A and C bases)
```

computed on the +DMS stops (pairing entry `i + 1` with base `i`, 0-based). A transcript with
coverage ≥ 1 has on average one +DMS stop per A/C position. Coverage ≥ 1 is the threshold
StructureFold analyses commonly use for trusting a transcript's reactivities, and the one heat
correction uses to pick its shared transcripts ([pipeline.md §10](pipeline.md#10-optional-heat-correction)).

## 10. Methods paragraph

Adapt to the options you used; the bracketed parts are config-dependent.

> Structure-seq reads were processed with structure-fold2-snakemake, a Snakemake workflow around
> StructureFold2 (Tang et al. 2015). Adapters were removed with cutadapt [5′ adapter
> TGAACAGCGACTAGGCTCTTCA, 3′ adapter GATCGGAAGAGCACACGTCTG], bases below Q30 were trimmed, and
> reads shorter than 20 nt were discarded. Reads were aligned end to end to [the transcriptome,
> plus the P4-P6 positive control and rRNA sequences] with Bowtie2 (Langmead and Salzberg 2012),
> keeping every alignment (-a). Alignments on the reverse strand, with more than [3] mismatches or
> indels, or with a mismatch at the first base were removed. RT stops were counted at the
> nucleotide 5′ of each read's start, and the +DMS and −DMS libraries of [each replicate] were
> summed. Reactivities were calculated with StructureFold2's rtsc_to_react: natural-log
> transformed stop counts were scaled by the transcript length over the transcript's total, the
> −DMS values were subtracted, negative values were set to 0, and the result was normalized
> per transcript by the mean of the 2–10% most reactive A and C positions and capped at 7
> (Ding et al. 2014; Ritchey et al. 2017).

Add, if enabled:

> **3′ bias:** Before the reactivity calculation, +DMS stops were corrected for 3′ coverage bias.
> On [2,500] well-covered protein-coding transcripts, log(+DMS / −DMS read depth) was modelled as
> natural cubic splines in distance from the 3′ end and from the 5′ end plus relative position;
> each position's +DMS stops were divided by the fitted bias, and the total was rescaled to the
> raw total.
>
> **Depth matching:** To compare libraries of different depth, each transcript's +DMS stops were
> randomly thinned (multivariate hypergeometric) to the smallest total among the libraries,
> reactivity was recalculated, and the mean over [100] draws was reported.
>
> **G/U normalization:** Instead of scaling each channel by its mean over all positions, the
> natural-log transformed +DMS and −DMS stop counts were each divided by their mean over the
> transcript's G and U positions (which DMS does not modify on the Watson-Crick face; transcripts
> with ≥ [20] G/U positions), so that both channels were on their background scale before the
> −DMS values were subtracted.

Software to cite: StructureFold2 (Tang et al. 2015, *Bioinformatics* 31:2668), Bowtie2 (Langmead
and Salzberg 2012, *Nat Methods* 9:357), cutadapt (Martin 2011, *EMBnet.journal* 17:10) or fastp
(Chen et al. 2018, *Bioinformatics* 34:i884), SAMtools (Danecek et al. 2021, *GigaScience*
10:giab008), Snakemake (Mölder et al. 2021, *F1000Research* 10:33), and, when used, gffread
(Pertea and Pertea 2020, *F1000Research* 9:304), AGAT, HISAT2 (Kim et al. 2019, *Nat
Biotechnol* 37:907) and plotnine. For the method: Ding et al. 2014 (*Nature* 505:696), Ritchey
et al. 2017 (*Nucleic Acids Res* 45:e135, Structure-seq2), Su et al. 2018 (*PNAS* 115:12170,
heat correction).
