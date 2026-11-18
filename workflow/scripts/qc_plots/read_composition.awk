# Reads per target type from a SAM body (samtools view, no header; unmapped reads included), grouped by
# read (bowtie2 -a writes all alignments of a read together). Only forward-strand alignments assign a
# type, as in the reactivity (sam_filter drops reverse-strand alignments, -F 20).
#   total_reads, unmapped_reads, mapped_reads (any alignment, either strand)
#   unique_reads / multimapped_reads: exactly one / more than one alignment (either strand)
#   antisense_only: mapped reads with no forward-strand alignment (no type)
#   unique_<type> : reads with exactly one alignment, forward strand, by type
#   typed_<type>  : reads whose forward alignments all hit transcripts of that one type
#   type_ambiguous: reads whose forward alignments hit more than one type
#   pie_<slice>   : every read in exactly one slice: unmapped; otherwise the first type in -v pie
#                   (priority order, e.g. "rRNA mRNA lncRNA") its forward alignments hit; otherwise other
#                   (other ncRNA, positive control, unannotated, antisense-only)
# -v s=<sample> -v types="<type list>" -v pie="<priority list>"; first file: "<id>\t<type>".
function flush(   i, slice) {
    if (q == "") return
    reads++
    if (unm) { unmapped++; pies["unmapped"]++; return }
    mapped++
    if (n == 1) unique++; else multi++
    if (nf == 0) { anti++; pies["other"]++; return }
    if (n == 1) uniq[first]++
    if (ncat == 1) typed[first]++; else ambiguous++
    slice = "other"
    for (i = 1; i <= kp; i++) if (pt[i] in seen) { slice = pt[i]; break }
    pies[slice]++
}
NR == FNR { cat[$1] = $2; next }
$1 != q { flush(); q = $1; n = 0; nf = 0; ncat = 0; unm = 0; delete seen }
{
    if (int($2 / 4) % 2) { unm = 1; next }              # flag 4: unmapped
    n++
    if (int($2 / 16) % 2) next                          # flag 16: reverse strand, no type
    nf++
    c = ($3 in cat) ? cat[$3] : "other"
    if (!(c in seen)) { seen[c] = 1; ncat++; first = c }
}
BEGIN {
    if (types == "") types = "mRNA ncRNA positive_control other"
    kp = split(pie, pt, " ")
}
END {
    flush()
    print "sample", "metric", "reads"
    print s, "total_reads", reads + 0
    print s, "unmapped_reads", unmapped + 0
    print s, "mapped_reads", mapped + 0
    print s, "unique_reads", unique + 0
    print s, "multimapped_reads", multi + 0
    print s, "antisense_only", anti + 0
    print s, "type_ambiguous", ambiguous + 0
    k = split(types, t, " ")
    for (i = 1; i <= k; i++) print s, "unique_" t[i], uniq[t[i]] + 0
    for (i = 1; i <= k; i++) print s, "typed_" t[i], typed[t[i]] + 0
    for (i = 1; i <= kp; i++) print s, "pie_" pt[i], pies[pt[i]] + 0
    print s, "pie_other", pies["other"] + 0
    print s, "pie_unmapped", pies["unmapped"] + 0
}
