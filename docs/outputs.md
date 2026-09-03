# Outputs

The result of a run is one tab-separated file of calls plus twelve QC figures. Everything
else written along the way is intermediate state.

- [The calls](#the-calls)
- [Mutation-string notation](#mutation-string-notation)
- [Zygosity](#zygosity)
- [Decision rules](#decision-rules)
- [Figures](#figures)
- [Intermediate files](#intermediate-files)

## The calls

`<results_dir>/<run_name>_mutation_statuses.txt` — tab-separated, one header line, one row per cell barcode, rows sorted by barcode:

```
barcode	maternal	paternal
ACACCCTGTCTCACCT	wt	wt
AGCTCTCTCCAAATGC	D20994977-20994987,D20995488-20995499	low_cov_uncalled
```

Each haplotype field is either a status word or a comma-separated mutation signature,
ordered by reference position.

| Value                | Meaning                                                                          |
| -------------------- | -------------------------------------------------------------------------------- |
| a mutation signature | The mutations called for that chromosome, comma-separated                        |
| `wt`                 | High coverage, but no mutation observed with signal above the calling threshold. |
| `low_cov_uncalled`   | Too little coverage to decide                                                    |
| `n/a`                | That haplotype was never observed for this cell                                  |

> **A barcode absent from this file is not wild type.** A row is written only if at least one of its haplotype fields is a call — rows where every field is `n/a`, or `low_cov_uncalled` are dropped. The file is therefore not a census of the cells you supplied. Compare against your input barcode list if you need one.

## Mutation-string notation

`scptseq.misc.parse_mutation` is definitional for this notation.

| Form | Example | Meaning |
|---|---|---|
| `{ref}{pos}{alt}` | `G20995489A` | Substitution: reference base `G` at 20995489 observed as `A` |
| `D{start}-{end}` | `D20995488-20995499` | Deletion of the reference span, **both ends inclusive** — this example is 12 bp |
| `I{pos}{bases}` | `I20995489ACGT` | Insertion of `ACGT` immediately after position 20995489 |
| `J{start}_{end}` | `J20994993_20995101` | A standard splice junction; the bounds are the first and last intronic bases |
| `J{s}_{e}>{Δs}_{Δe}` | `J20994993_20995101>-26_+0` | A junction offset from the standard one: add the deltas to get the observed intron, here `(20994967, 20995101)` |
| `J{s}_{e}>-` | `J20994993_20995101>-` | The standard junction is **absent** from this chromosome's reads |

The `>-` form is emitted only by the `call` stage. Similarly, `call` re-resolves
target-adjacent deletions against the junction catalogue, so the same event can be recorded
as a deletion in `*_stat_cntrs_before_and_after.yml` and as an offset junction in the final
calls. When cross-reading those two files, expect that translation.

## Zygosity

The zygosity pie assigns each called cell one category. Two of the names are vestigial:

- **`hom` is 2PC and does not mean the two chromosomes carry identical mutations.** It means neither chromosome was called `wt`. 
- **`het maternal` names the chromosome that carries the mutation**, i.e. maternal is
  mutated and paternal is `wt`. `het paternal` is the reverse.

The figure relabels these: `het` becomes `1PC` (one perturbed copy), `hom` becomes `2PC`, so the slices read `wt`, `1PC`, `2PC`, `uncalled and wt`, and `uncalled and PC`.

## Decision rules

All thresholds in source. Editing them makes results non-comparable with a run at the defaults.

| Decision                                            | Rule                                                                              |
| --------------------------------------------------- | --------------------------------------------------------------------------------- |
| Call `mut` (emit a signature)                       | `max_frac ≥ 0.3` **and** at least 3 reads with defining mutation                  |
| Call `wt`                                           | `max_frac ≤ 0.3` **and** `cov ≥ 10`                                               |
| Otherwise                                           | `low_cov_uncalled`                                                                |
| A mutation joins the signature                      | within 10 bp of a cut site **and** `frac > 0.5`                                   |
| A standard junction is marked absent (`>-`)         | `frac < 0.05` **and** `cov ≥ 3`                                                   |
| A junction enters the control catalogue             | ≥ 2 UMIs in at least one control cell **and** median per-cell UMI fraction > 0.05 |
| A reference gap is recorded as a junction candidate | span longer than 15, so 17 bases at minimum                                       |

Notes on reading this table:

- `max_frac` is the largest coverage fraction among non-splice mutations within 10 bp of
  a cut site.
- `cov` counts a read if either endpoint of the mutation falls inside the read's outer
  alignment bounds.
- When no target-adjacent mutation exists at all, `cov` falls back to the haplotype's total
  read count. That fallback is what permits a `wt` call in the absence of any candidate.

## Figures

Written to `<results_dir>/figures/`. Two come from `refsplice`, ten from `call`. The
per-target figures have one panel per CRISPR target.

| File | Plot | Axes |
|---|---|---|
| `*_umis_per_control_cell.pdf` | histogram | x: UMIs per control cell; y: cells (unlabelled) |
| `*_fraction_umis_cdf.pdf` | CDF, one curve per junction | x: fraction of UMIs carrying the junction; y: CDF |
| `*_max_frac_v_cov.pdf` | scatter, colour = number of chromosomes | x: `max_frac`; y: `cov`, log scale |
| `*_reads_per_cell_hist.pdf` | stepped histogram | x: reads per cell; y: cells |
| `*_reads_per_cell_haplotype.pdf` | stepped histogram per haplotype, before and after overlaid | x: reads per cell; y: cells |
| `*_most_common_muts.pdf` | horizontal bar, top 10 signatures plus a combined remainder | x: chromosomes |
| `*_mut_pie_chart.pdf` | pie, cell count in the title | zygosity categories |
| `*_target_pie.pdf` | pie, chromosome count in the title | which targets are mutated |
| `*_muts_by_type.pdf` | bar, per target | x: `sub`, `del`, `ins`, `splice`; y: chromosomes |
| `*_del_lens.pdf` | bar, per target | x: deletion length; y: chromosomes |
| `*_ins_lens.pdf` | bar, per target | x: insertion length; y: chromosomes |
| `*_comb_indel_lens.pdf` | bar | x: combined indel length per chromosome; y: chromosomes |

**`*_max_frac_v_cov.pdf` is the decision plot.** Every chromosome appears as one point, and
the dashed grey guides mark the rules in the table above: the vertical line at
`max_frac = 0.3`, the horizontal line at `cov = 10`, and the curve for the
evidence-count rule. Points to the right of the vertical line and above the curve are
called `mut`; points to its left and above the horizontal line are called `wt`; and the rest
are uncalled.

## Intermediate files

The YAML files written by `refsplice` and `count` are internal state, not a published
format. No schema is documented for them, and it may change.

If you need to inspect one: the junction files
(`*_splicing_junctions.yml`, `*_standard_splicing_junctions.yml`,
`*_splicing_junction_str.yml`) load with `yaml.FullLoader`, while
`*_stat_cntrs_before_and_after.yml` and `*_splicing_junction_frac_umis_per_cell.yml`
serialize `collections.Counter` objects and require the unsafe `yaml.Loader`:

```python
import yaml
stats = yaml.load(open('results/RUN_stat_cntrs_before_and_after.yml'), Loader=yaml.Loader)
```

Its top level is keyed by barcode, then haplotype, then `before` / `after` /
`defining mutation`. `before` counts all reads for that haplotype; `after` counts only
those carrying the defining mutation.
