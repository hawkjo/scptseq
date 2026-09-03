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

## User documentation lives in `docs/`

Behavior, formats, and interpretation are documented for users — read these
before re-deriving anything from source:

- `docs/tutorial.md` — the pipeline end to end, including the external
  alignment step and the command that produced the example BAMs.
- `docs/inputs.md` — target-file schema (and which fields the code actually
  reads), BAM layout, read-name contract, genome FASTA requirements.
- `docs/outputs.md` — the calls, the mutation-string notation, every threshold
  in the decision rules, and the figure inventory.
- `docs/troubleshooting.md` — error paths and silent-failure symptoms.

## Commands

- Install (editable): `pip install -e .` — or `uv sync` (a `uv.lock` is checked in; dev deps `pytest`, `ruff`).
- Lint + format: `ruff check . --fix && ruff format .`
- Run tests: `uv run pytest -m slow` — the suite is currently one end-to-end
  regression test, so `pytest -m "not slow"` selects nothing and exits 5.
- Run the CLI: `scptseq <command> ...` (entry point is `scptseq.main:main`). The
  docopt usage string in `scptseq/main.py` is the single source of truth for
  arguments; the README deliberately does not duplicate it.
- Sanity-check by hand: run the scripts in `examples/` in pipeline order from inside the `examples/` directory — they use the bundled `example_data/` and run in ~2 min total.

There is **no Cython** in this repo despite the templated mentions of it — all
pure Python. Don't go looking for `.pyx` files.

## Pipeline structure

Four subcommands, dispatched in `main.py` to a single top-level function each.
`refsplice`, `count` and `call` must share both `--results-dir` and
`<run_name>`; each reads YAML written by the previous stage and raises a clear
error otherwise.

1. **`splitfastqs`** (`split_fastq_by_bc.py`) — one FASTQ per cell barcode.
2. **`refsplice`** (`preprocessing.py:find_all_ref_splice_junctions`) — from
   **control** BAMs, builds the real-junction catalogue. Writes four YAMLs:
   `*_splicing_junctions.yml`, `*_standard_splicing_junctions.yml`,
   `*_splicing_junction_frac_umis_per_cell.yml`, `*_splicing_junction_str.yml`,
   plus two figures.
3. **`count`** (`preprocessing.py:haplotyped_mutation_preprocessing`) — from
   **perturbed** BAMs, writes per-cell `mutations/*.mutations.yml` plus
   `*_stat_cntrs_before_and_after.yml`.
4. **`call`** (`call_mutations.py:call_mutations`) — thresholds the counts into
   `*_mutation_statuses.txt` and ten figures.

`config.py:CommandLineArguments` wraps the raw docopt dict; access args via its
properties (e.g. `arguments.run_name`), never the dict directly.

## Gotchas for development

- The mutation/splice-junction matching logic in `preprocessing.py` and
  `call_mutations.py` (transitive matching to standard junctions, the
  `before`/`after` counters, distance-to-cutsite filtering) is the load-bearing,
  easy-to-get-silently-wrong part. Any tests should target these with small
  hand-built inputs where the answer is known by construction.
- `misc.parse_mutation` is the single source of truth for mutation-string
  notation. Read it before touching anything mutation-related.
- In `call_mutations.py`, `yy_thresh` is assigned inside the
  `max_frac_v_cov` figure block and consumed by `get_cell_haplotype_status`
  further down. The calling logic depends on that figure code having run.
- **Output behavior is frozen while the paper is under revision.** Verify any
  change by diffing a fresh example run against
  `examples/example_data/ref_TX46_Prosalpha3_mutation_statuses.txt` (the
  reference is unsorted; `sort` both sides).
