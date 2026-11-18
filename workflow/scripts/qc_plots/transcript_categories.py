"""Transcript ID -> target type for the read-composition QC (transcript_categories rule, report.smk).

Sources, later ones override earlier ones:
  --base FASTA          base transcriptome: mRNA, except records whose header carries an Ensembl-style
                        transcript_biotype:<x> tag (Ensembl cDNA/ncRNA FASTA): rRNA biotypes -> rRNA,
                        other non-protein_coding biotypes -> ncRNA
  --base-annotation PREFIX ANNOTATION
                        (repeatable) biotypes for base records built from a GTF/GFF3 (gffread headers
                        carry none): transcript <PREFIX><id> -> mRNA (protein_coding), rRNA (rRNA
                        biotypes) or ncRNA (any other biotype); records without a biotype stay mRNA.
                        PREFIX is "-" for none (genotype mode: "<genotype>__").
  --set TYPE FASTA [ANNOTATION]
                        every record of FASTA -> TYPE (ncRNA, lncRNA, rRNA, positive_control). In an
                        ncRNA set, a transcript_biotype:<x> tag in the FASTA header (Ensembl ncRNA FASTA)
                        gives rRNA for rRNA biotypes and mRNA for protein_coding (e.g. a set of organelle
                        transcripts); transcripts marked as rRNA -> rRNA also, for a set built from an annotation
                        (GFF3/GTF), GTF transcript_biotype/gene_biotype, GFF3 biotype attributes, or a
                        GFF3 feature type of rRNA (RAP-DB style), matched by transcript ID.
Writes "<id>\t<type>" lines to --out. Windows line endings in the FASTAs are ignored.
"""
import argparse
import gzip
import re

RRNA_BIOTYPES = {"rrna", "rrna_pseudogene", "rrna_gene"}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def headers(path):
    with opener(path) as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].rstrip("\r\n")


def rrna_ids(path):
    """Transcript IDs an annotation marks as rRNA."""
    out = set()
    gtf_tid = re.compile(r'transcript_id "([^"]+)"')
    gtf_bt = re.compile(r'(?:transcript_biotype|gene_biotype|transcript_type|gene_type) "([^"]+)"')
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.rstrip("\r\n").split("\t")
            if len(col) < 9:
                continue
            ftype, attr = col[2], col[8]
            m = gtf_tid.search(attr)
            if m:                                          # GTF
                b = gtf_bt.findall(attr)
                if any(x.lower() in RRNA_BIOTYPES for x in b) or ftype.lower() == "rrna":
                    out.add(m.group(1))
                continue
            a = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)   # GFF3
            bts = [a.get(k, "") for k in ("biotype", "transcript_biotype", "gene_biotype", "type")]
            if ftype.lower() == "rrna" or any(x.lower() in RRNA_BIOTYPES for x in bts):
                for k in ("ID", "transcript_id"):
                    if k in a:
                        out.add(a[k])
                        out.add(a[k].split(":", 1)[-1])      # Ensembl GFF3 "transcript:<id>"
    return out


def biotypes(path):
    """Transcript ID -> biotype from a GTF (transcript_biotype / transcript_type, else gene_*) or GFF3."""
    out = {}
    gtf_tid = re.compile(r'transcript_id "([^"]+)"')
    pref = [re.compile(k + r' "([^"]+)"') for k in ("transcript_biotype", "transcript_type", "gene_biotype", "gene_type")]
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.rstrip("\r\n").split("\t")
            if len(col) < 9:
                continue
            ftype, attr = col[2], col[8]
            m = gtf_tid.search(attr)
            if m:                                          # GTF: transcript (or first exon) rows
                if m.group(1) in out:
                    continue
                for r in pref:
                    b = r.search(attr)
                    if b:
                        out[m.group(1)] = b.group(1); break
                continue
            a = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)   # GFF3
            tid = a.get("ID") or a.get("transcript_id")
            if not tid or "Parent" in a and ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                continue
            b = next((a[k] for k in ("transcript_biotype", "biotype", "transcript_type") if k in a), None)
            if b is None and ftype.lower() == "rrna":
                b = "rRNA"
            if b is None and ftype == "mRNA":
                b = "protein_coding"
            if b:
                out[tid] = b; out[tid.split(":", 1)[-1]] = b
    return out


def category(biotype):
    b = biotype.lower()
    return "mRNA" if b == "protein_coding" else ("rRNA" if b in RRNA_BIOTYPES else "ncRNA")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--base-annotation", nargs=2, action="append", default=[], metavar=("PREFIX", "ANNOTATION"))
    ap.add_argument("--set", nargs="+", action="append", default=[], metavar=("TYPE", "FASTA"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cat = {}
    for h in headers(args.base):
        tid, c = h.split()[0], "mRNA"
        m = re.search(r"transcript_biotype:(\S+)", h)
        if m and m.group(1) != "protein_coding":
            c = "rRNA" if m.group(1).lower() in RRNA_BIOTYPES else "ncRNA"
        cat[tid] = c
    for prefix, ann in args.base_annotation:
        prefix = "" if prefix == "-" else prefix
        bt = biotypes(ann)
        n = {}
        for tid, b in bt.items():
            if prefix + tid in cat:
                c = category(b); cat[prefix + tid] = c; n[c] = n.get(c, 0) + 1
        print(f"{ann}: base records labelled from biotypes: {n}")
    for s in args.set:
        typ, fasta, ann = s[0], s[1], (s[2] if len(s) > 2 else None)
        rr = rrna_ids(ann) if (ann and typ == "ncRNA") else set()
        n = n_rr = 0
        for h in headers(fasta):
            tid = h.split()[0]
            m = re.search(r"transcript_biotype:(\S+)", h) if typ == "ncRNA" else None
            c = "rRNA" if tid in rr else (category(m.group(1)) if m else typ)
            n += 1; n_rr += c != typ
            cat[tid] = c
        print(f"{fasta}: {n} records -> {typ}" + (f" ({n_rr} relabelled from biotype tags/annotation)" if n_rr else ""))
    with open(args.out, "w") as f:
        for tid, c in cat.items():
            f.write(f"{tid}\t{c}\n")
    counts = {}
    for c in cat.values():
        counts[c] = counts.get(c, 0) + 1
    print("types:", counts)


if __name__ == "__main__":
    main()
