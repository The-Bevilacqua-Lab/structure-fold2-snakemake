"""
Self-contained HTML QC report for one structure-fold2-snakemake run (rule qc_report).

Every table is written inline and every figure is embedded as a base64 PNG, so the single
.html file can be opened, archived or e-mailed on its own. Optional sections (positive
control, 3' bias correction, depth matching, heat correction) appear only when their inputs
were given. Runs via Snakemake's `script:` directive (reads the `snakemake` object).
"""
import base64
import datetime
import html
import io
import os
import sys

import pandas as pd
from plotnine import (aes, coord_cartesian, element_blank, element_text, facet_wrap, geom_col, ggplot, guide_legend, guides, labs,
                      scale_fill_manual, scale_y_continuous, theme, theme_bw)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # workflow/scripts
import plot_style  # noqa: E402,F401  Arial in every figure


sm = snakemake  # noqa: F821  (injected by Snakemake)
sys.stderr = open(sm.log[0], "w")
inp = dict(sm.input.items())
P = sm.params
TYPES = list(P.composition_types)
DISPLAY_SCALE = 0.85   # on-screen size relative to each figure's saved physical size
TYPE_LABEL = {"mRNA": "mRNA", "rRNA": "rRNA", "ncRNA": "ncRNA", "lncRNA": "lncRNA", "positive_control": "positive control",
              "other": "other"}
TYPE_COLOR = {"mRNA": "#2a6fb0", "rRNA": "#d1495b", "ncRNA": "#e08a1e", "lncRNA": "#b5487a", "positive_control": "#3a9e5c",
              "other": "#9a9a9a"}

CSS = """
/* Layout: logo + title header over a green rule; sticky section index on the left, report on the right.
   One light theme on purpose (white, with the logo's green and red), so no dark-mode variant. */
:root { --bg:#ffffff; --fg:#1c2620; --muted:#5d6b62; --line:#dde5df; --panel:#eef7f1;
        --green:#14823c; --green-deep:#0d5c2a; --red:#d42a1e; --figbg:#ffffff;
        --font-body:"Source Sans 3", -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        --font-data:"IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace; color-scheme: light; }
* { box-sizing: border-box; }
body { margin:0; font:15px/1.55 var(--font-body); color:var(--fg); background:var(--bg);
       padding-left:max(16px, env(safe-area-inset-left, 0px)); padding-right:max(16px, env(safe-area-inset-right, 0px)); }
header { display:flex; align-items:center; gap:20px; flex-wrap:wrap; padding-block:18px 14px;
         border-bottom:3px solid var(--green); }
header img.logo { height:72px; width:auto; flex:none; }
header .titles { min-width:0; }
header h1 { margin:0 0 2px; font-size:24px; font-weight:600; color:var(--green-deep); text-wrap:balance; }
header p { margin:0; color:var(--muted); font-size:14px; }
.layout { display:flex; align-items:flex-start; gap:28px; }
nav { position:sticky; top:env(safe-area-inset-top, 0px); width:200px; flex:none; padding-block:18px;
      max-height:100vh; overflow:auto; }
nav a { display:block; padding:3px 0 3px 10px; color:var(--muted); text-decoration:none; font-size:14px;
        border-left:2px solid transparent; }
nav a:hover, nav a:focus-visible { color:var(--green-deep); border-left-color:var(--green); }
main { flex:1; min-width:0; padding-block:4px 56px; max-width:1180px; }
section { padding-block:18px 18px; border-bottom:1px solid var(--line); }
h2 { font-size:19px; font-weight:600; margin:8px 0 10px; color:var(--green-deep); text-wrap:balance; }
h2::after { content:""; display:block; width:36px; height:2px; margin-top:6px; background:var(--red); }
h3 { font-size:15px; font-weight:600; margin:16px 0 6px; color:var(--fg); }
p.note { color:var(--muted); margin:4px 0 10px; font-size:14px; max-width:75ch; }
.tablewrap { overflow-x:auto; margin:8px 0 12px; }
table { border-collapse:collapse; font-size:13px; font-variant-numeric:tabular-nums; }
th, td { border-bottom:1px solid var(--line); padding:4px 10px; text-align:right; white-space:nowrap; }
th { background:var(--panel); color:var(--green-deep); font-weight:600; }
th:first-child, td:first-child { text-align:left; }
td { font-family:var(--font-data); font-size:12.5px; }
td:first-child, .kv td { font-family:var(--font-body); font-size:13px; }
figure { margin:10px 0 16px; }
figure img { display:block; max-width:100%; height:auto; border:1px solid var(--line); border-radius:4px;
             background:var(--figbg); }
figcaption { color:var(--muted); font-size:13px; margin-top:4px; }
.grid { display:grid; grid-template-columns:repeat(auto-fill, minmax(min(380px, 100%), 1fr)); gap:14px; }
.grid > * { min-width:0; }
.kv td:first-child { color:var(--muted); }
.kv td { text-align:left; }
a { color:var(--green); }
a:focus-visible { outline:2px solid var(--green); outline-offset:2px; }
code { font-family:var(--font-data); font-size:0.92em; }
footer { color:var(--muted); font-size:12px; padding-block:20px; border-top:1px solid var(--line); }
@media (max-width: 900px) { nav { display:none; } header img.logo { height:56px; } }
"""


# ----------------------------------------------------------------- helpers
def esc(x):
    return html.escape(str(x))


def img_tag(path_or_plot, caption=None):
    """Embed a PNG file, or a plotnine ggplot rendered on the fly, as a base64 <img>."""
    if isinstance(path_or_plot, str):
        if not os.path.exists(path_or_plot):
            return f'<p class="note">Missing figure: {esc(path_or_plot)}</p>'
        with open(path_or_plot, "rb") as f:
            data = f.read()
    else:
        import matplotlib.pyplot as plt   # plotnine's backend, only to close the rendered figure
        fig = path_or_plot.draw()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
        plt.close(fig)
        data = buf.getvalue()
    cap = f"<figcaption>{esc(caption)}</figcaption>" if caption else ""
    # Show every figure at its physical size (pixels / saved DPI, at 96 CSS px per inch), capped by the
    # column, so figures saved at different resolutions don't appear at arbitrary relative sizes.
    style = ""
    try:
        from PIL import Image   # ships with matplotlib
        with Image.open(io.BytesIO(data)) as im:
            dpi = float((im.info.get("dpi") or (96, 96))[0]) or 96.0
            style = f' style="width:{im.width / dpi * 96 * DISPLAY_SCALE:.0f}px"'
    except Exception:
        pass
    return (f'<figure><img alt="{esc(caption or "")}"{style} '
            f'src="data:image/png;base64,{base64.b64encode(data).decode()}">{cap}</figure>')


def table(df, index=True, floatfmt="{:,.2f}", intcols=()):
    df = df.copy()
    for c in df.columns:
        if c in intcols or pd.api.types.is_integer_dtype(df[c]):
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{int(v):,}")
        elif pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    return '<div class="tablewrap">' + df.to_html(index=index, border=0, escape=True) + "</div>"


def kv_table(d):
    rows = "".join(f"<tr><td>{esc(k)}</td><td>{esc(v)}</td></tr>" for k, v in d.items())
    return f'<div class="tablewrap"><table class="kv">{rows}</table></div>'


sections = []   # (anchor, title, body html)


def add(anchor, title, body):
    sections.append((anchor, title, body))


# ----------------------------------------------------------------- run summary
samples = pd.DataFrame(P.samples)
cols = [c for c in ["sample", "condition", "ID", "temperature", "run", "r1"] if c in samples.columns]
samples = samples[cols]
body = "<h3>Settings</h3>" + kv_table(P.settings)
body += f"<h3>Samples ({samples['sample'].nunique()} libraries, {len(samples)} sequencing runs)</h3>"
body += table(samples, index=False)
body += f'<p class="note">Reactivity IDs: {esc(", ".join(P.ids))}</p>'
add("run", "Run summary", body)

# ----------------------------------------------------------------- mapping
aln = pd.read_csv(inp["alignment"], sep="\t").set_index("sample")
body = '<p class="note">From bowtie2 (all sequencing runs of a sample combined), before mismatch filtering.</p>'
body += table(aln, intcols=("total_reads",))
ALN_LEVELS = ["not aligned", "aligned >1 times", "aligned exactly once"]   # first level = top of the stack
aln_long = pd.concat([
    pd.DataFrame({"sample": aln.index, "category": "aligned exactly once", "pct": aln["unique_mapped_pct"].to_numpy()}),
    pd.DataFrame({"sample": aln.index, "category": "aligned >1 times", "pct": aln["multi_mapped_pct"].to_numpy()}),
    pd.DataFrame({"sample": aln.index, "category": "not aligned", "pct": (100 - aln["overall_alignment_pct"]).to_numpy()}),
], ignore_index=True)
aln_long["sample"] = pd.Categorical(aln_long["sample"], categories=list(aln.index))
aln_long["category"] = pd.Categorical(aln_long["category"], categories=ALN_LEVELS)
p_aln = (ggplot(aln_long, aes("sample", "pct", fill="category"))
         + geom_col(width=0.75)
         + scale_fill_manual(values={"aligned exactly once": "#2a6fb0", "aligned >1 times": "#8fb6dd",
                                     "not aligned": "#d9dee5"})
         + guides(fill=guide_legend(reverse=True))
         + scale_y_continuous(expand=(0, 0)) + coord_cartesian(ylim=(0, 100))   # zoom, never drop bars
         + labs(x="", y="% of reads", fill="")
         + theme_bw(base_size=10)
         + theme(figure_size=(max(5, 0.55 * len(aln) + 2.5), 3.8), legend_position="top",
                 axis_text_x=element_text(rotation=45, ha="right"),
                 panel_grid_major_x=element_blank(), panel_grid_minor=element_blank()))
body += img_tag(p_aln, "Read alignment per sample")
add("mapping", "Read mapping", body)

# ----------------------------------------------------------------- composition
comp = pd.read_csv(inp["composition"], sep="\t").set_index("sample")
PIE = list(P.pie_slices)
body = ('<p class="note"><b>Pie</b>: every read that entered Bowtie2, in one slice. A mapped read goes to the '
        f'first of {" &gt; ".join(TYPE_LABEL[t] for t in PIE[:-2])} that its forward-strand alignments hit; '
        'reads hitting none of those (other ncRNA, positive control, unannotated, or antisense-only) are '
        '<i>other</i>. The last pie pools all samples. The second set of pies shows the mapped reads only '
        '(antisense-only reads are mapped, and count as other).</p>')
body += img_tag(inp["composition_pie"], "Where the reads go: all reads")
body += img_tag(inp["composition_pie_mapped"], "Where the reads go: mapped reads only")
body += ('<p class="note">Only forward-strand alignments give a type, as in the reactivity (reverse-strand '
         'alignments are dropped by the mismatch filter); mapped reads with none are <i>antisense-only</i>. '
         '<b>Unique</b>: reads with exactly one alignment, as % of uniquely mapped reads. '
         '<b>By type</b>: reads whose forward alignments all fall in one target type (reads shared between isoforms '
         'of an mRNA or between rRNA copies are kept), as % of all mapped reads; reads hitting more than one '
         'type are <i>type-ambiguous</i>. Multicopy RNAs (rRNA, tRNA) are mostly multimapped, so the unique '
         'measure under-counts ncRNA. Counts are from the merged alignments before mismatch filtering.</p>')
show = comp[["total_reads", "mapped_pct", "mapped_reads", "unique_reads", "unique_pct_of_mapped"]
            + [f"unique_{t}_pct" for t in TYPES] + [f"typed_{t}_pct" for t in TYPES]
            + ["type_ambiguous_pct", "antisense_only_pct"]].copy()
show.columns = (["reads", "% mapped", "mapped reads", "unique reads", "% unique"]
                + [f"unique: % {TYPE_LABEL[t]}" for t in TYPES]
                + [f"by type: % {TYPE_LABEL[t]}" for t in TYPES] + ["% type-ambiguous", "% antisense-only"])
body += table(show, intcols=("reads", "mapped reads", "unique reads"))
PANELS = ["Uniquely mapped reads (% of unique)", "All mapped reads, by type (% of mapped)"]
rows = []
for prefix, panel in zip(["unique", "typed"], PANELS):
    for t in TYPES:
        v = comp[f"{prefix}_{t}_pct"].fillna(0)
        if t == "other" and v.sum() == 0:
            continue
        rows.append(pd.DataFrame({"sample": comp.index, "type": TYPE_LABEL[t], "pct": v.to_numpy(), "panel": panel}))
rows.append(pd.DataFrame({"sample": comp.index, "type": "type-ambiguous",
                          "pct": comp["type_ambiguous_pct"].fillna(0).to_numpy(), "panel": PANELS[1]}))
rows.append(pd.DataFrame({"sample": comp.index, "type": "antisense-only",
                          "pct": comp["antisense_only_pct"].fillna(0).to_numpy(), "panel": PANELS[1]}))
comp_long = [d for d in rows if d.pct.sum() > 0 or d.type.iloc[0] == "mRNA"]   # drop types absent everywhere
type_levels = ["antisense-only", "type-ambiguous"] + [TYPE_LABEL[t] for t in reversed(TYPES)]   # mRNA at the bottom
comp_long = pd.concat(comp_long, ignore_index=True)
comp_long["type"] = pd.Categorical(comp_long["type"], categories=type_levels)
comp_long["sample"] = pd.Categorical(comp_long["sample"], categories=list(comp.index))
comp_long["panel"] = pd.Categorical(comp_long["panel"], categories=PANELS)
fills = {TYPE_LABEL[t]: TYPE_COLOR[t] for t in TYPES}
fills["type-ambiguous"] = "#c9c9c9"
fills["antisense-only"] = "#ececec"
p_comp = (ggplot(comp_long, aes("sample", "pct", fill="type"))
          + geom_col(width=0.75)
          + facet_wrap("panel", ncol=2)
          + scale_fill_manual(values=fills, breaks=[x for x in reversed(type_levels) if x in set(comp_long["type"])])
          + scale_y_continuous(expand=(0, 0)) + coord_cartesian(ylim=(0, 100))   # zoom, never drop bars
          + labs(x="", y="% of reads", fill="")
          + theme_bw(base_size=10)
          + theme(figure_size=(max(9, 1.1 * len(comp) + 4), 4), legend_position="right",
                  axis_text_x=element_text(rotation=45, ha="right"),
                  panel_grid_major_x=element_blank(), panel_grid_minor=element_blank()))
body += img_tag(p_comp, "Mapped reads by target type")
add("composition", "Read composition by target type", body)

# ----------------------------------------------------------------- genome check
if "genome_check_summary" in inp:
    gc = pd.read_csv(inp["genome_check_summary"], sep="\t").set_index("sample")
    body = ('<p class="note">Reads Bowtie2 could not align to the transcriptome, aligned to the genome with HISAT2 '
            '(splice-aware, end-to-end, up to 10 loci per read) and placed by their 5′ end. '
            '<b>transcriptome exon</b>: in an exon of a transcript the reads were aligned to (typically an '
            'unannotated splice junction or indel); <b>&lt;biotype&gt; exon</b>: exon of an annotated transcript '
            'left out of the transcriptome; <b>intron</b>, <b>antisense</b> (inside a gene on the other strand only) '
            'and <b>intergenic</b>. Category percentages are of the transcriptome-unaligned reads; BAMs are in '
            '<code>genome_check/</code>.</p>')
    pct_cols = [c for c in gc.columns if c.endswith("_pct")]
    show = gc[["total_reads", "transcriptome_unaligned", "pct_of_total_unaligned", "genome_aligned",
               "pct_unaligned_on_genome", "genome_aligned_unique"] + pct_cols].copy()
    show.columns = (["total reads", "not on transcriptome", "% of total", "on genome", "% on genome", "on genome, 1 locus"]
                    + [f"% {c[:-4]}" for c in pct_cols])
    body += table(show, intcols=("total reads", "not on transcriptome", "on genome", "on genome, 1 locus"))
    body += img_tag(inp["genome_check_plot"], "Transcriptome-unaligned reads by genome category")
    loci = pd.read_csv(inp["genome_check_loci"], sep="\t", keep_default_na=False)
    if len(loci):
        body += "<h3>Loci with the most transcriptome-unaligned reads (top 20)</h3>"
        cols = ["rank", "chrom", "strand", "start", "end", "reads", "pct_unique", "category", "gene_id", "gene_name"]
        body += table(loci[cols].head(20).set_index("rank"), floatfmt="{:,.1f}", intcols=("start", "end", "reads"))
        body += f'<p class="note">All {len(loci):,} top loci, with per-sample counts: <code>genome_check/top_loci.tsv</code>.</p>'
    add("genome-check", "Transcriptome-unaligned reads on the genome", body)

# ----------------------------------------------------------------- positive control
if "positive_control" in inp or "p4p6_all" in inp:
    body = ""
    if "positive_control" in inp:
        pc = pd.read_csv(inp["positive_control"], sep="\t").set_index("sample")
        body += ('<p class="note">Share of mapped alignments (primary and secondary) on the positive-control '
                 'sequence, per sample.</p>')
        body += table(pc, intcols=("total_mapped", "positive_control_reads"), floatfmt="{:,.4f}")
    steps = [name for key, name in [("3_prime_bias_correction", "3′ bias-corrected"), ("depth_matching", "depth-matched")]
             if P.settings.get(key)]
    if "p4p6_all" in inp or "p4p6_roc" in inp:
        body += ('<p class="note">Final reactivities' + (f" ({', '.join(steps)})" if steps else " (no 3′ bias "
                 "correction or depth matching)") + ', i.e. the ones written for downstream use.</p>')
    if "p4p6_all" in inp:
        body += img_tag(inp["p4p6_all"], "p4p6 reactivity on the P4-P6 secondary structure, every reactivity ID")
    if "p4p6_roc" in inp:
        body += ('<h3>Reactivity vs known structure</h3><p class="note">ROC over A/C residues: positive = unpaired, '
                 'negative = in a Watson–Crick–Franklin pair (A–U, C–G); A/C in other pairs are left out. '
                 'High reactivity predicts unpaired. AUC = chance an unpaired A/C is more reactive than a WCF-paired one.</p>')
        if "p4p6_auc" in inp:
            body += table(pd.read_csv(inp["p4p6_auc"]).set_index("label"), floatfmt="{:,.4f}")
        body += img_tag(inp["p4p6_roc"], "p4p6 ROC: unpaired vs WCF-paired A/C")
    add("positive-control", "Positive control", body)

# ----------------------------------------------------------------- specificity
rows = []
for path in inp["specificity_tables"]:
    t = pd.read_csv(path)
    name = os.path.basename(path)[len("specificity_"):-len(".csv")]
    frac = dict(zip(t.iloc[:, 0], t.iloc[:, 2]))
    rows.append({"sample": name, **{f"% {b}": 100 * frac.get(b, float("nan")) for b in "ACGT"}})
spec = pd.DataFrame(rows).set_index("sample")
spec["% A+C"] = spec["% A"] + spec["% C"]
cond = samples.drop_duplicates("sample").set_index("sample")["condition"]
spec.insert(0, "condition", cond.reindex(spec.index))
body = ('<p class="note">Base at each RT stop (per sample). DMS modifies A and C, so +DMS libraries should be '
        'enriched for A/C stops relative to −DMS.</p>')
body += table(spec)
body += img_tag(inp["specificity_plot"], "RT-stop base specificity per sample")
add("specificity", "RT-stop specificity", body)

# ----------------------------------------------------------------- coverage
rows = []
for path in inp["coverage"]:
    t = pd.read_csv(path)
    v = pd.to_numeric(t.iloc[:, 1], errors="coerce")
    rows.append({"reactivity ID": os.path.basename(os.path.dirname(path)), "transcripts": int(v.notna().sum()),
                 "coverage ≥ 1": int((v >= 1).sum()), "coverage ≥ 2": int((v >= 2).sum()),
                 "median coverage": float(v.median())})
body = ('<p class="note">+DMS RT-stop coverage per transcript (A/C RT stops per A/C position, as in '
        '<code>coverage.csv</code>).</p>')
body += table(pd.DataFrame(rows).set_index("reactivity ID"))
GENOTYPE_NAMES = list(getattr(P, "genotype_names", []) or [])


def imgs(key, cap):
    """img_tag for one figure, or one per genotype when the input is a list (genotype mode)."""
    if key not in inp:
        return ""
    paths = inp[key] if isinstance(inp[key], (list, tuple)) else [inp[key]]
    if len(paths) == 1 and not GENOTYPE_NAMES:
        return img_tag(paths[0], cap)
    out = ""
    for path in paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        g = next((g for g in GENOTYPE_NAMES if stem.endswith("_" + g)), "")
        out += img_tag(path, f"{cap} ({g})" if g else cap)
    return out


for key, cap in [("upset_replicates", "Transcripts with coverage ≥ 1, shared between replicates"),
                 ("upset_merged", "Transcripts with coverage ≥ 1, shared between pooled groups")]:
    body += imgs(key, cap)
add("coverage", "Transcript coverage", body)

# ----------------------------------------------------------------- replicate correlation
if "rt_corr" in inp or "rt_corr_minus" in inp:
    figs = imgs("rt_corr", "+DMS RT-stop correlation between replicates") + \
        imgs("rt_corr_minus", "−DMS RT-stop correlation between replicates")
    add("correlation", "Replicate correlation", f'<div class="grid">{figs}</div>')

# ----------------------------------------------------------------- genotype comparison
if "genotype_summary" in inp:
    lift = pd.read_csv(inp["genotype_liftover"], sep="\t", keep_default_na=False)
    ident = pd.to_numeric(lift["identity"], errors="coerce")
    ga, gb = GENOTYPE_NAMES[:2]
    body = (f'<p class="note">{esc(ga)} and {esc(gb)} each map to their own transcriptome. Transcripts are '
            f'matched {"through transcript_matching" if P.settings.get("transcript_matching", "identical IDs") != "identical IDs" else "by identical ID"}; '
            'equal-length pairs with ≥ 80% identity are mapped position by position, others are aligned '
            '(gap-free blocks only). Reactivity per genotype is the mean over its reactivity IDs, compared at '
            'A/C positions reactive in both.</p>')
    body += kv_table({"matched transcript pairs": f"{len(lift):,}",
                      "mapped 1:1 (ungapped)": f"{(lift['method'] == 'ungapped').sum():,}",
                      "aligned": f"{(lift['method'] == 'aligned').sum():,}",
                      "too long to align": f"{(lift['method'] == 'too_long').sum():,}",
                      "median sequence identity": f"{ident.median():.4f}",
                      "−DMS pooled across genotypes": str(P.settings.get("pool_minus_across_genotypes", False))})
    body += "<h3>Per-transcript reactivity agreement</h3>"
    body += table(pd.read_csv(inp["genotype_summary"], sep="\t").set_index("comparison"), floatfmt="{:,.3f}",
                  intcols=("transcripts",))
    body += ('<p class="note">Single-replicate rows are the like-for-like baseline: if r between genotypes '
             '(one replicate each) is close to r within a genotype, most transcripts fold alike in both.</p>')
    body += img_tag(inp["genotype_plot"], "Per-transcript Pearson r between and within genotypes")
    add("genotypes", "Genotype comparison", body)

# ----------------------------------------------------------------- metagene
body = ('<p class="note">Per transcript, ln(coverage + 1) divided by its transcript mean, averaged by relative '
        'position across transcripts; one line per sample.</p>')
body += img_tag(inp["metagene_depth"], "Metagene: read depth (samtools depth)")
body += img_tag(inp["metagene_rtsc"], "Metagene: RT stops")
add("metagene", "Metagene coverage", body)

# ----------------------------------------------------------------- 3' bias correction
if "bias_models" in inp:
    rows = []
    for path in inp["bias_models"]:
        m = pd.read_csv(path, sep="\t", index_col=0, dtype={"value": str})["value"]
        rows.append({"reactivity ID": path.split(os.sep)[-3], "model form": m.get("model_form", "linear"),
                     "R²": float(m.get("r2", "nan")),
                     "fit positions": int(float(m.get("n_fit_positions", "nan"))),
                     "training transcripts": int(float(m.get("n_training_transcripts", "nan")))})
    body = ('<p class="note">The model is fit on protein-coding training transcripts and applied to every '
            "transcript's +DMS RT-stops (mRNA, ncRNA and the positive control), so every reactivity, p4p6 "
            'plot and ROC in this report is built from corrected stops.</p>')
    body += table(pd.DataFrame(rows).set_index("reactivity ID"), floatfmt="{:,.3f}")
    for cov, rt in zip(inp["bias_coverage_plots"], inp["bias_rtsc_plots"]):
        rid = cov.split(os.sep)[-3]
        body += f"<h3>{esc(rid)}</h3>"
        body += f'<div class="grid">{img_tag(cov, "Coverage log-ratio, fit and residual")}{img_tag(rt, "+DMS RT-stop metagene before / after")}</div>'
    add("bias", "3′ bias correction", body)

# ----------------------------------------------------------------- depth matching
if "depth_matched_corr" in inp or "depth_matched_summary" in inp:
    body = ('<p class="note">Each library\'s +DMS RT-stops are thinned, per transcript, to the smallest '
            '+DMS total across reactivity IDs (and the −DMS likewise when depth_matching_match_minus is on); '
            'reactivity is the mean over the draws.'
            + (' Genotype mode: each matched pair is thinned to one total across the libraries of both '
               'genotypes and listed under its first genotype\'s transcript; unmatched transcripts are not '
               'depth-matched.' if "genotype_summary" in inp else "") + '</p>')
    if "depth_matched_summary" in inp:
        dm = pd.read_csv(inp["depth_matched_summary"], sep="\t")   # one row per transcript
        rows = []
        n_draws = dm.filter(like="_draws_defined").to_numpy().max()   # reported if >= half the draws defined
        for col in [c for c in dm.columns if c.endswith("_plus_total")]:
            rid = col[: -len("_plus_total")]
            total, defined = dm[col], dm[f"{rid}_draws_defined"]
            kept = (dm["target_plus"] / total).where(total > 0)
            rows.append({"reactivity ID": rid, "transcripts matched": len(dm),
                         "+DMS stops before": int(total.sum()), "+DMS stops after": int(dm["target_plus"].sum()),
                         "median fraction kept": kept.median(),
                         "transcripts never thinned": int((dm["target_plus"] >= total).sum()),
                         "transcripts with reactivity": int((defined >= 0.5 * n_draws).sum())})
        body += table(pd.DataFrame(rows).set_index("reactivity ID"), floatfmt="{:,.3f}",
                      intcols=("transcripts matched", "+DMS stops before", "+DMS stops after",
                               "transcripts never thinned", "transcripts with reactivity"))
    if "depth_matched_corr" in inp:
        body += img_tag(inp["depth_matched_corr"], "+DMS RT-stop correlation after depth matching (first draw)")
    add("depth-matching", "Depth-matched reactivity", body)

# ----------------------------------------------------------------- heat correction
if "heat_shared" in inp:
    with open(inp["heat_shared"]) as f:
        n = sum(1 for line in f if line.strip())
    add("heat", "Heat correction",
        f'<p>{n:,} transcripts have +DMS RT-stop coverage ≥ 1 at both temperatures and were used for the '
        f'heat correction.</p>')

# ----------------------------------------------------------------- links
if P.multiqc_link:
    add("links", "Other reports",
        f'<p><a href="{esc(P.multiqc_link)}">MultiQC report of fastp trimming</a> '
        f'(separate file, next to this report).</p>')

# ----------------------------------------------------------------- write
# Section order in the report (sections are built above in whatever order is convenient)
ORDER = ["run", "mapping", "composition", "genome-check", "specificity", "correlation", "positive-control", "genotypes", "coverage",
         "metagene", "bias", "depth-matching", "heat", "links"]
sections.sort(key=lambda sec: ORDER.index(sec[0]) if sec[0] in ORDER else len(ORDER))
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
logo_html = ""
if inp.get("logo") and os.path.exists(inp["logo"]):
    with open(inp["logo"], "rb") as f:
        logo_html = (f'<img class="logo" alt="structure-fold2-snakemake logo" '
                     f'src="data:image/png;base64,{base64.b64encode(f.read()).decode()}">')
navlinks = "".join(f'<a href="#{a}">{esc(t)}</a>' for a, t, _ in sections)
body = "".join(f'<section id="{a}"><h2>{esc(t)}</h2>{b}</section>' for a, t, b in sections)
doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{esc(P.title)}</title><meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400&family=Source+Sans+3:wght@400;600&display=swap">
<style>{CSS}</style></head>
<body>
<header>{logo_html}<div class="titles"><h1>{esc(P.title)}</h1><p>{esc(P.settings.get("output_dir", ""))} · generated {now}</p></div></header>
<div class="layout"><nav>{navlinks}</nav><main>{body}</main></div>
<footer>structure-fold2-snakemake · self-contained report (figures embedded)</footer>
</body></html>"""
with open(sm.output[0], "w") as f:
    f.write(doc)
print(f"wrote {sm.output[0]} ({len(doc) / 1e6:.1f} MB, {len(sections)} sections)", file=sys.stderr)
