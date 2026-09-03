# scPT-seq: Single-Cell Perturbation and Transcriptome Sequencing

A suite of tools for processing scPT-seq data, as described in
### Allele-resolved CRISPR editing and splicing outcomes with their transcriptional consequences in vivo

**John A. Hawkins, Siamak Redhai, Svenja Leible, Mireia Osuna Lopez, Hilal Ozgur, Tianyu Wang, Michaela Holzem, Lars M. Steinmetz, Michael Boutros, Oliver Stegle**

`scptseq` calls CRISPR-induced mutations — deletions, insertions, substitutions and
splicing changes — per cell and per haplotype near targeted cut sites. It takes aligned
reads and produces per-cell, per-haplotype mutation-status calls and QC figures.

Read alignment is not part of this package. Align with the mapper of your choice,
between the first and second stages of the pipeline.

## Data requirements

1. **Long reads that span a segregating site and a cut site.** A read contributes to a call only if
   it covers both. Haplotypes are assigned by reading the base at a segregating site whose maternal and paternal alleles you supply. A homozygous sample cannot be haplotyped.
2. **A matched control sample of the same gene.** The control is what distinguishes
   perturbation-induced splicing changes from the background of low-frequency
   non-canonical junctions.

The input format requirements — one coordinate-sorted, indexed BAM per cell, barcode and UMI
encoded in the read name, spliced alignment — are described in [docs/inputs.md](docs/inputs.md).

## Installation

For the easiest installation, use pip:

```
python -m pip install git+https://github.com/hawkjo/scptseq
```

To install from a clone, optionally into a virtual environment:

```
git clone https://github.com/hawkjo/scptseq.git
cd scptseq
python -m venv envscptseq
. envscptseq/bin/activate
python -m pip install .
```

For a reproducible environment, use [uv](https://docs.astral.sh/uv/) instead — a
`uv.lock` is checked in, and uv is the only installer that reads it:

```
uv sync
```

Dependencies are listed in `pyproject.toml`. Linux and macOS; on Windows, use WSL, since
pysam publishes no Windows wheels. Install time ~1.5 minutes.

### Verifying the installation

The bundled end-to-end test runs the last three stages against the example data and
compares the result against a checked-in reference output:

```
uv run pytest -m slow
```

It takes about two minutes. 

## The pipeline

Five steps, in order. 

1. **`scptseq splitfastqs`** — split reads into one FASTQ per cell barcode. Run it separately on
   the control and the perturbed reads.
2. **Align, externally.** Spliced alignment, one BAM per cell, coordinate-sorted and
   indexed. See [docs/tutorial.md](docs/tutorial.md) for the command used to produce the
   bundled example data.
3. **`scptseq refsplice`** — catalogue the real splice junctions, canonical and
   non-canonical, from the **control** BAMs.
4. **`scptseq count`** — collate per-cell, per-haplotype, per-mutation read and UMI
   statistics from the **perturbed** BAMs.
5. **`scptseq call`** — threshold those counts into mutation-status calls and QC figures.

| Command | Reads | Writes |
|---|---|---|
| `splitfastqs` | one or more FASTQs | one FASTQ per barcode |
| `refsplice` | control BAMs, target info | 4 junction YAMLs, 2 figures |
| `count` | perturbed BAMs, genome FASTA, target info, `refsplice` output | per-cell mutation YAMLs, 1 aggregate statistics YAML |
| `call` | `refsplice` and `count` output, target info | `*_mutation_statuses.txt`, 10 figures |

Two things to know before your first run:

- **`refsplice`, `count` and `call` must share both `--results-dir` and `<run_name>`.**
  Each stage reads what the previous one wrote, and both values are part of how it finds
  those files.
- **The tool is silent by default.** Pass `-v`, `-vv` or `-vvv` to see progress.

For the full argument lists, run `scptseq --help`.

## Documentation

- **[docs/tutorial.md](docs/tutorial.md)** — start here. Run the bundled example, then
  adapt it to your own reads, targets and genome.
- **[docs/inputs.md](docs/inputs.md)** — the input format reference: target information,
  BAM layout, read names, genome FASTA.
- **[docs/outputs.md](docs/outputs.md)** — how to read the calls, the mutation-string
  notation, the decision rules, and the QC figures.
- **[docs/troubleshooting.md](docs/troubleshooting.md)** — symptoms and their causes.

Upstream barcode and UMI extraction is handled by the separate
[freediv10Xcellbcs](https://github.com/hawkjo/freediv10Xcellbcs) tool.

## Examples

Runnable scripts for each stage, using the bundled example data for 30 cells, are in
[examples/](examples/). Runtime ~2 minutes.
