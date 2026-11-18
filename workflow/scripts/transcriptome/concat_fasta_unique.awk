# Concatenate FASTAs, keeping only the first record of each ID (first word of the header).
# Later records with an ID already seen are dropped (e.g. an ncRNA/lncRNA annotation that overlaps
# the base annotation, or an lncRNA also in the ncRNA set), so Bowtie2 never sees duplicate
# references. Windows line endings (CR) are stripped. Drop counts per file go to stderr.
FNR == 1 { if (prev != "" && dropped[prev]) printf "%s: %d duplicate IDs dropped\n", prev, dropped[prev] > "/dev/stderr"; prev = FILENAME }
{ sub(/\r$/, "") }
/^>/ {
    id = substr($1, 2)
    keep = !(id in seen)
    if (keep) seen[id] = 1; else dropped[FILENAME]++
}
keep { print }
END { if (prev != "" && dropped[prev]) printf "%s: %d duplicate IDs dropped\n", prev, dropped[prev] > "/dev/stderr" }
