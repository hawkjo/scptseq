# Tutorial

This walks from the bundled example to a run on your own data, in five parts. Read it in
order the first time.

- [Part 1 — run the bundled example](#part-1--run-the-bundled-example)
- [Part 2 — what just happened](#part-2--what-just-happened)
- [Part 3 — align your own reads](#part-3--align-your-own-reads)
- [Part 4 — write your own target file](#part-4--write-your-own-target-file)
- [Part 5 — read your results](#part-5--read-your-results)

## Part 1 — run the bundled example

Run from inside `examples/`; the paths in the scripts are relative to it.

```
cd examples
./splitfastqs_example.sh
./refsplice_example.sh
./count_example.sh
./call_example.sh
```

The first script splits the bundled FASTQ into one file per cell barcode. The remaining
three do not consume its output — the aligned perturbed BAMs are already checked in, so
the example skips the alignment step. That step is Part 3.

The three analysis stages take about two minutes together. They share
`--results-dir=results` and the run name `TX46_Prosalpha3`, and produce:

```
results/TX46_Prosalpha3_splicing_junctions.yml
results/TX46_Prosalpha3_standard_splicing_junctions.yml
results/TX46_Prosalpha3_splicing_junction_frac_umis_per_cell.yml
results/TX46_Prosalpha3_splicing_junction_str.yml
results/TX46_Prosalpha3_stat_cntrs_before_and_after.yml
results/TX46_Prosalpha3_mutation_statuses.txt
results/figures/                       (12 PDFs)
example_data/TX46_Prosalpha3_bams/mutations/   (30 per-cell YAMLs)
```

To check the run, look at the final file: `results/TX46_Prosalpha3_mutation_statuses.txt`. You should see 30 barcode rows with `wt`, `low_cov_uncalled`, and distinct mutation signatures. The result should match the checked-in `example_data/ref_TX46_Prosalpha3_mutation_statuses.txt`.

`count` writes a `mutations/` directory inside the perturbed BAM directory, so delete it after a test run.

### Running the per-cell stages in parallel

`refsplice` and `count` do one independent unit of work per cell, so both take a
`--threads` option that spreads that work over worker processes:

```
./refsplice_example.sh --threads=8
```

Output does not depend on the value: the per-cell results are reduced in file order, not
completion order, so every file is byte-identical at any thread count. The bundled example
is too small to gain from it — thirty cells is less work than starting the workers — but on
a real run of thousands of cells this is most of the wall-clock time. `count` benefits more
than `refsplice`. Each worker holds its own copy of the target contig, so memory grows with
the thread count.

## Part 2 — what just happened

| Stage | Read | Wrote | The number to look at |
|---|---|---|---|
| `refsplice` | 30 control BAMs | the junction catalogue, 2 figures | how many junctions were found, and how many became "standard" |
| `count` | 30 perturbed BAMs, genome FASTA, the catalogue | per-cell mutation lists, one aggregate statistics file | reads per cell per haplotype |
| `call` | the catalogue and the statistics file | the calls, 10 figures | the `wt` / uncalled / signature split |

On the bundled control data `refsplice` finds 5 junctions and promotes 2 of them to
standard. With `-vvv` it prints the table it built.

Three local definitions:

- **A "splice junction" here is any gap in the reference alignment** longer than the
  threshold — the code does not consult the CIGAR operation, so an `N` and a long `D` are
  indistinguishable at this stage. Whether a given gap ends up reported as a junction or
  as a deletion depends on whether it matches the catalogue built from the control sample.
- **"Standard" means prevalent in the control sample**, not canonical in the annotation:
  a junction supported by at least 2 UMIs in some control cell, whose median per-cell UMI
  fraction exceeds 5%.
- **"before" and "after"** in the statistics file are two views of the same haplotype's
  reads: *before* is all of them, *after* is only those carrying that haplotype's defining
  mutation — the target-adjacent mutation covering the largest fraction of its reads.

## Part 3 — align your own reads

Alignment sits between `splitfastqs` and `refsplice`, and this package does not do it. Run
it on both the control and the perturbed FASTQs.

The bundled example BAMs were produced with minimap2 2.22 and samtools 1.21, one cell at a
time:

```
minimap2 -w 4 -k 13 -ax splice -t 20 <genome.fa> <barcode>.fq > <barcode>.sam
samtools view -b <barcode>.sam | samtools sort - > <barcode>.sorted.bam
samtools index <barcode>.sorted.bam
```

The published run used a whole-genome index rather than the single-contig FASTA that
ships in `example_data/`.

Whatever aligner you use, its output must satisfy four properties. Each one has a distinct
failure mode:

| Requirement                                                                               | If it is missing                                                                                                      |
| ----------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| **Spliced alignment** (reference gaps represented as gaps, not clipped or forced through) | Splicing changes are misread or invisible                                                                             |
| **One cell per BAM, named `<barcode>.sorted.bam`**                                        | The barcode is taken from the filename, so calls are mislabelled or the cell is skipped                               |
| **Coordinate-sorted and indexed**                                                         | The stages fetch only the target region, which requires an index. Without one, pysam raises on open                   |
| **Read names preserved verbatim**                                                         | The barcode and UMI are parsed out of the read name. See [inputs.md](inputs.md) for the specification and what breaks |

Contig names must agree across three places: the `chrm` field of your
target file, the BAM header, and the record ID in your genome FASTA. The BAM and FASTA should match automatically from the alignment step.

## Part 4 — write your own target file

> **All coordinates in this file are 0-based**, matching the BAM internals that pysam
> exposes. Coordinates copied from IGV, Ensembl, or a VCF are 1-based and must be
> converted. This is the single most consequential field to get wrong: the pipeline will
> run to completion, write every figure, and produce a mutation-statuses file containing
> nothing but its header line.

A minimal entry needs eight values. Everything else in the bundled
`example_data/targets.yml` is provenance that the code never reads.

```yaml
MyGene:
  chrm: 2R                  # must match the BAM header and the FASTA record ID
  start: 20994900           # region to analyse, half-open [start, end)
  end: 20996118
  targets:
    t1:
      cutsite: 20995489     
    t2:
      cutsite: 20994988
  seg_info:
    seg_sites:              # positions that distinguish the two haplotypes
      - 20994903
      - 20995104
    seg_bases:              # [maternal, paternal] base at each site, same order
      - [G, C]
      - [T, G]
    skip_sites: []          # positions to ignore
```

The gene key is what you pass as `<gene_name>` on the command line. The file is read with
`yaml.safe_load`, so write plain YAML with no type tags.

`seg_sites` and `seg_bases` are paired by position in the list, so their order must
correspond. See [inputs.md](inputs.md) for what the maternal/paternal assignment means and
how haplotypes are decided.

## Part 5 — read your results

`*_mutation_statuses.txt` is a tab-separated file with one row per cell barcode and one
column per haplotype:

```
barcode	maternal	paternal
ACACCCTGTCTCACCT	wt	wt
CGTAGCGAGCCACGCT	D20994977-20994987,D20995488-20995499	D20994985-20994987,D20995489-20995523
CCATTCGAGAATAGGG	J20994993_20995101>-26_+0,D20995489-20995497	D20994983-20994987,D20995489-20995491
```

A field is either a status word or a comma-separated **mutation signature**. In this gene
the two cut sites are at 20995489 (`t1`) and 20994988 (`t2`).

A couple examples:

**Row 2, maternal — `D20994977-20994987,D20995488-20995499`.** Two deletions, one per
target. `D20994977-20994987` is an 11 bp deletion ending one base before the `t2` cut
site; `D20995488-20995499` is a 12 bp deletion spanning `t1`. Both chromosomes of this
cell are edited at both targets.

**Row 3, maternal — `J20994993_20995101>-26_+0,D20995489-20995497`.** The `J` form names a
junction *relative to a standard one*: the standard junction is `(20994993, 20995101)`, and
the offsets shift its start by −26 and its end by 0, so the intron actually observed runs
from 20994967 to 20995101 — starting 26 bases upstream of the normal donor, and spanning
the `t2` cut site at 20994988. Alongside it, a 9 bp deletion at `t1`.

Two status words appear in this example:

- **`wt`**
- **`low_cov_uncalled`** — too few reads covering the relevant position to decide.

**A barcode missing from this file is not wild type.** Cells whose every haplotype field
would be uncalled are omitted from the output entirely, so the file is not a census of the
cells you put in.

[outputs.md](outputs.md) has the complete notation table, the decision rules and their
thresholds, and what each of the 12 figures plots.
