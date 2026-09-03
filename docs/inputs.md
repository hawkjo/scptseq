# Inputs

Reference for everything you supply: the target file, the BAMs, the read names, and the
genome FASTA.

- [Coordinates are 0-based](#coordinates-are-0-based)
- [The target file](#the-target-file)
- [Segregating sites and haplotypes](#segregating-sites-and-haplotypes)
- [BAM files](#bam-files)
- [Read names](#read-names)
- [Genome FASTA](#genome-fasta)
- [splitfastqs input](#splitfastqs-input)
- [Running more than one gene](#running-more-than-one-gene)

## Coordinates are 0-based

> **Every coordinate in the target file is 0-based**, matching the BAM internals that pysam
> exposes. `start` and `end` describe a half-open interval `[start, end)`; `cutsite`, `seg_sites`
> and `skip_sites` are single positions compared directly against pysam reference positions.
>
> Coordinates taken from IGV, Ensembl, or a VCF are 1-based and must be converted.
>
> **The failure is silent.** With coordinates off by one, the segregating-site bases match neither
> the maternal nor the paternal allele, so no read can be haplotyped. Every stage exits
> successfully, every figure is written, and `*_mutation_statuses.txt` contains its header line and
> nothing else. If you see that, check this first.

## The target file

A YAML mapping keyed by gene name. The key is what you pass as `<gene_name>`. These values
per gene are read by the code:

```yaml
MyGene:
  chrm: 2R
  start: 20994900
  end: 20996118
  targets:
    t1:
      cutsite: 20995489
    t2:
      cutsite: 20994988
  seg_info:
    seg_sites:
      - 20994903
      - 20995104
    seg_bases:
      - [G, C]
      - [T, G]
    skip_sites: []
```

| Field                                                                                                             | Read by the code | Meaning                                                                               |
| ----------------------------------------------------------------------------------------------------------------- | ---------------- | ------------------------------------------------------------------------------------- |
| `chrm`                                                                                                            | yes              | Contig name. Must match the BAM header and the FASTA record ID                        |
| `start`, `end`                                                                                                    | yes              | Region analysed, half-open. Reads are fetched for this interval only                  |
| `targets.t<n>.cutsite`                                                                                            | yes              | Cut site of each target. `t1` is required, and is also the reference point for ordering segregating sites |
| `seg_info.seg_sites`                                                                                              | yes              | Positions that distinguish the haplotypes                                             |
| `seg_info.seg_bases`                                                                                              | yes              | `[maternal, paternal]` base at each site, paired with `seg_sites` by list position    |
| `seg_info.skip_sites`                                                                                             | yes              | Individual reference positions to ignore when collecting mutations                    |
| `primer`, `strand`, `targets.t*.start` / `end` / `seq`, `seg_info.seg_signatures`, `seg_info.seg_sites_and_bases` | **no**           | Present in the bundled example as provenance. Never consulted, omit them freely       |

Three notes on writing your own:

- **Plain YAML only.** The file is read with `yaml.safe_load`. `examples/example_data/targets.yml` is a working example.
- **Targets must be named `t1`, `t2`, ... `tN`, and `t1` must be present.** A gene may
  define any number of them. Unlike the rest of the file, `targets` accepts no other keys.
  Anything else under it is an error, not an ignored field. Keep cut sites more than 20 bp apart,
  and inside `[start, end)`. 
- **`skip_sites` is a flat list of individual positions**, tested by membership against one
  reference position at a time. An interval will not work — a two-element list is read as
  two positions, not as a range. Use `[]` for none.

## Segregating sites and haplotypes

A **segregating site** is a position where your two chromosomes carry different bases, and
you know which base belongs to which. If you do not know which base is maternal versus paternal, an arbitrary choice suffices but this must be remembered for downstream interpretation.

How a read is assigned:

1. The segregating sites are ordered by distance from the `t1` cut site.
2. The nearest site that the read covers is consulted.
3. If the observed base equals the first base in that site's `seg_bases` pair, the read is
   maternal, if it equals the second, paternal.
4. If it matches neither, the next-nearest covered site is tried, and so on.
5. A read matching no site at any position is excluded from the analysis.

## BAM files

- **One cell per file**, in a directory per sample. `refsplice` takes the control
  directory; `count` takes the perturbed one.
- **Named `<barcode>.sorted.bam`.** In `count` the barcode is taken from the filename —
  everything before the first `.` — so the filename determines how the cell is labelled in
  the output.
- **Coordinate-sorted and indexed.** Only the target region is fetched, which requires the
  index alongside the BAM.
- Files are matched by a `*bam` glob, so put nothing else ending in `bam` in the directory.
- No SAM tag is read anywhere. A single pooled BAM with `CB`/`UB` tags will not work; the
  reads must be split into per-cell files with the barcode and UMI in the read name.

## Read names

```
{barcode}_{umi}#{original_read_name}
```

Everything before the first `#` is split on `_`, giving exactly two fields. A real example
from the bundled data:

```
ACACCCTGTCTCACCT_CCGTCGGGTA#m64368e_250322_211917/1444384/ccs
```

The original read name may itself contain `_` and `#`; only the first `#` matters. 

If your read names do not follow this requirement, the two stages fail differently:
`refsplice` unpacks the split into two variables and raises a `ValueError`, while `count`
takes the second field by index and will silently use the wrong string as the UMI. A
barcode or UMI containing `_` breaks both.

## Genome FASTA

Required by `count` only.

- Uncompressed. It is opened directly, so a `.gz` will not be read.
- Loaded fully into memory.
- The record ID must match `chrm` in your target file and the contig name in the BAM header.
- It need not be a whole genome — a single-contig FASTA is fine. The bundled example ships
  chromosome 2R alone against BAMs whose headers list 1870 contigs.

## splitfastqs input

Takes the barcoded, UMI-tagged FASTQs produced upstream by [freediv10Xcellbcs](https://github.com/hawkjo/freediv10Xcellbcs), and writes one `<barcode>.fq` per barcode.

- `<fastq_files>` is **one argument**: a comma-separated list, not several arguments.
- Inputs ending in `.gz` are decompressed transparently.
- The barcode is the text before the first `_` in the read ID.
- **It exits if the output directory already exists**, rather than overwriting. Remove it
  or choose another path.
- `--output-dir` is not `~`-expanded; use an absolute or relative path.

Note that the output is `<barcode>.fq` while the analysis stages expect
`<barcode>.sorted.bam` — the alignment step in between produces that name.

## Running more than one gene

A target file can define many genes, and the bundled example defines two. Each run handles
one gene, and outputs are keyed only by `<run_name>` within a single flat results
directory. So for a second gene, change `<run_name>`, `--results-dir`, or both.

With `count`, its per-cell files are named `<barcode>.mutations.yml` with no gene or run name, and they are written into the perturbed BAM directory. Two genes processed over the same BAM directory will overwrite each other's per-cell files, and the stage re-reads every file in that directory when it aggregates. Use a separate BAM directory per gene, or delete `mutations/` between runs.
