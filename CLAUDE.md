# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project: scptseq

**What it is:** A command-line suite for processing **scPT-seq** (single-cell
Perturbation and Transcriptome sequencing) data. It calls CRISPR-induced
mutations (indels, substitutions, and splicing changes) per cell, per
haplotype, near targeted cut sites — starting from aligned reads and ending in
per-cell mutation-status calls plus QC figures. Accompanies the Hawkins et al.
scPT-seq paper.

**Scope:** This package handles only the post-alignment steps: barcode
demultiplexing of FASTQs, splice-junction reference building, per-cell mutation
counting, and mutation calling. **Read alignment is done externally** (by the
user's mapper of choice) between `splitfastqs` and the rest of the pipeline.
Upstream barcode/UMI extraction is done by the separate
[freediv10Xcellbcs](https://github.com/hawkjo/freediv10Xcellbcs) tool.

## Commands

- Install (editable): `pip install -e .` — or `uv sync` (a `uv.lock` is checked in; dev deps `pytest`, `ruff`).
- Lint + format: `ruff check . --fix && ruff format .`
- Run tests: `pytest` (no `tests/` directory exists yet — add under `tests/` per the global preferences).
- Run the CLI: `scptseq <command> ...` (entry point is `scptseq.main:main`). See `scptseq/main.py` for the docopt usage string.
- Sanity-check by hand: run the scripts in `examples/` in pipeline order from inside the `examples/` directory — they use the bundled `example_data/` and run in ~2 min total.

There is **no Cython** in this repo despite the templated mentions of it — all
pure Python. Don't go looking for `.pyx` files.

## The pipeline (commands run in this order)

The four subcommands are stages of one pipeline and **must share the same
`--results-dir`** — later stages read YAML files written by earlier ones (and
raise a clear error if `--results-dir` doesn't match). Each is dispatched in
`main.py` to a single top-level function:

1. **`splitfastqs`** (`split_fastq_by_bc.py`) — splits one or more FASTQs into one file per cell barcode. Run on both control and perturbed reads, then align externally.
2. **`refsplice`** (`preprocessing.py:find_all_ref_splice_junctions`) — from **control** BAMs, discovers all real splice junctions (canonical + non-canonical), so that perturbation-induced splicing changes can be distinguished from background. Writes `*_splicing_junctions.yml`, `*_standard_splicing_junctions.yml`, and `*_splicing_junction_str.yml`.
3. **`count`** (`preprocessing.py:haplotyped_mutation_preprocessing`) — from **perturbed** BAMs, collates per-cell/per-haplotype/per-mutation read & UMI statistics into `*_stat_cntrs_before_and_after.yml`. Needs the genome FASTA and the `refsplice` outputs.
4. **`call`** (`call_mutations.py:call_mutations`) — thresholds the counts into final per-cell mutation calls (`*_mutation_statuses.txt`) and writes the QC figures into `<results-dir>/figures/`.

`config.py:CommandLineArguments` wraps the raw docopt dict; access args via its
properties (e.g. `arguments.run_name`, `arguments.results_dir`), never the dict
directly.

## Key data formats and conventions

- **BAM input**: one sorted+indexed BAM per cell, named `<barcode>.sorted.bam`, in a per-sample directory. Code globs `*bam` and reads only the target gene's region via pysam `.fetch(chrm, start, end)`.
- **Read name format**: `{barcode}_{umi}#{original_read_name}`. Parsed by `umi_tools.bc_and_umi_given_read_name` — barcode and UMI are split on `_`, after stripping anything from `#` on.
- **Target info YAML** (`examples/example_data/targets.yml`): keyed by gene name. Per gene it holds `chrm/start/end`, two CRISPR targets `t1`/`t2` (each with a `cutsite`), and `seg_info` (segregating sites + their maternal/paternal bases). Loaded via `misc.TargetInfo`.
- **Haplotyping**: each read is assigned `maternal` or `paternal` (see `constants.haplotypes`) by reading the segregating site nearest the t1 cut site — `TargetInfo.maternal_or_paternal`.
- **Mutation string notation** (parsed by `misc.parse_mutation`, the single source of truth — read it before touching anything mutation-related):
  - `{ref}{pos}{alt}` substitution (e.g. `G20995489A`)
  - `D{start}-{end}` deletion
  - `I{pos}{bases}` insertion
  - `J{start}_{end}` standard splice junction; `J{s}_{e}>{Δstart}_{Δend}` a near-standard junction offset from a standard one; trailing `>-` means a standard junction is missing.
- **UMI correction** (`umi_tools.py`): UMIs within edit distance 1 are merged via connected components (scipy sparse graph), collapsing to the highest-count UMI in each component. Used both for splice-junction support counts and for fractional read/UMI tallies.
- **"before"/"after"** in the count/call structures: "before" = all reads for the haplotype; "after" = reads filtered to those carrying the haplotype-defining mutation (the target-adjacent mutation with the highest fraction of coverage).

## Gotchas

- `scipy` is imported in `umi_tools.py` but is **not** listed in `setup.py`'s `install_requires` (it comes in transitively via `statsmodels`). Add it explicitly if you touch dependency declarations.
- The mutation/splice-junction matching logic in `preprocessing.py` and `call_mutations.py` (transitive matching to standard junctions, the `before`/`after` counters, distance-to-cutsite filtering) is the load-bearing, easy-to-get-silently-wrong part. Any tests should target these with small hand-built inputs where the answer is known by construction.
