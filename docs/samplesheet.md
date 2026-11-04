## Samplesheet

Tab-separated, one row per (sample, sequencing run):

| column | required | meaning |
|---|---|---|
| `sample` | yes | sample identifier, matches your FASTQ naming |
| `condition` | yes | exactly `plus` (+DMS/treated) or `minus` (-DMS/untreated) |
| `r1` | yes | path to that (sample, run)'s raw FASTQ |
| `ID` | no | groups samples into biological replicates (e.g. `rep1`) -- a given ID needs both a `plus` and `minus` row somewhere in the sheet, since the reactivity calculation pairs one condition's ID against the other's. **If omitted**, defaults to each sample's row position *within its condition* (`rep1`, `rep2`, ... independently for `plus` and `minus`) -- i.e. the 1st `minus` row pairs with the 1st `plus` row as one replicate, the 2nd with the 2nd, and so on. This assumes equal `plus`/`minus` counts given in matching order; give an explicit `ID` column if that's not your layout. |
| `run` | no | distinguishes multiple sequencing runs of the same sample (e.g. resequenced on a second flow cell); same-sample rows get concatenated before trimming. **Defaults to a single run** if omitted. |
| `temperature` | no | probing temperature for that row (e.g. `30`, `55`, or labels like `low`/`high` -- any two distinct values). **If present**, enables heat correction (`react_heat_correct.py`) between the two temperature groups -- see Outputs below. Every row of a given `ID` must share the same temperature, and exactly two distinct values must appear across the whole sheet. Omit the column entirely to skip heat correction. |