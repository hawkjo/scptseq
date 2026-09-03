# Troubleshooting

Indexed by symptom. Section A is for runs that stop, section B is for runs that finish and
produce something wrong.

## A. It raised an exception

| Message                                                                                   | Cause                                                                                                                   | Fix                                                                            |
| ----------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `<dir> does not exist. results-dir must match count`                                      | `call` cannot find the results directory                                                                                | Run `count` first, and pass the same `--results-dir`                           |
| `<file>_splicing_junction_str.yml does not exist. results-dir must match refsplice`       | `count` or `call` cannot find the `refsplice` output                                                                    | Run `refsplice` first, with the same `--results-dir` and the same `<run_name>` |
| `<file>_standard_splicing_junctions.yml does not exist. results-dir must match refsplice` | as above                                                                                                                | as above                                                                       |
| `FileNotFoundError: ..._stat_cntrs_before_and_after.yml`                                  | `call` ran before `count`                                                                                               | Run `count` first                                                              |
| `Output directory already exists`                                                         | `splitfastqs` refuses to write into an existing directory                                                               | Remove it, or choose another `--output-dir`                                    |
| `ValueError: too many values to unpack` while reading control BAMs                        | A read name does not match `{barcode}_{umi}#{...}`                                                                      | See [read names](inputs.md#read-names)                                         |
| `ValueError: max() arg is an empty sequence` during `call`                                | No cell reached the calling stage                                                                                       | Check `count`'s `-vvv` output for how many reads were haplotyped               |
| `OSError: [Errno 24] Too many open files` during `splitfastqs`                            | One file handle is held open per barcode for the whole run, which exceeds the default limit at roughly a thousand cells | Raise the limit for the shell that runs it, e.g. `ulimit -n 4096`              |
| `<file> has no record named <chrm>` during `count`                                        | The genome FASTA has no record matching the target file's `chrm`                                                        | Check `chrm` against the FASTA record ids and the BAM header                   |
| `BrokenProcessPool` during `refsplice` or `count`                                         | A `--threads` worker died, most often killed for memory                                                                 | Lower `--threads`, or re-run with `--threads=1` to get the underlying error    |

**A note on `--threads` and error messages.** A worker process inherits none of the
parent's logging configuration, so a warning raised inside one prints without the
timestamp prefix that every other line carries. Re-running with `--threads=1` runs the
same code in the parent process, which gives both the normal formatting and an ordinary
traceback instead of a `BrokenProcessPool`.

**A note on the three "results-dir must match" messages.** The filenames they look for
embed `<run_name>`, so all three also fire when `--results-dir` is correct but
`<run_name>` was mistyped between stages. Check both before assuming the earlier stage
failed.

## B. It finished, but the output looks wrong

| Symptom                                                          | Possible cause                                                                                                                                                                                                       |
| ---------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `*_mutation_statuses.txt` has only its header line               | Coordinates in the target file are 1-based. Nothing can be haplotyped, so nothing is called, and no stage errors. See [inputs.md](inputs.md#coordinates-are-0-based)                                                 |
| Same symptom, coordinates confirmed 0-based                      | `chrm` does not match the BAM contig name or the FASTA record ID; or `start`/`end` do not overlap the reads; or `seg_bases` pairs are in the wrong order relative to `seg_sites`                                     |
| Every chromosome is `wt`                                         | The perturbed and control directories were swapped — the control sample has no edits at the cut sites                                                                                                                |
| The junction catalogue contains short "introns" at the cut sites | The perturbed and control directories were swapped: the perturbed sample's cut-site deletions were catalogued as if they were real junctions, so they are then treated as background                                 |
| The figures are dominated by substitutions                       | The genome FASTA is a different assembly or build from the one the reads were aligned against.                                                                                                                       |
| Cells appear that are not in your BAM directory                  | Stale per-cell files in `mutations/` from an earlier run. The stage re-reads every file in that directory. Delete it and re-run                                                                                      |
| The maternal and paternal columns look transposed                | `seg_bases` lists `[maternal, paternal]` per site. If yours are `[paternal, maternal]`, every call is on the wrong chromosome, and nothing detects it                                                                |
| It appears to hang                                               | The default log level is `ERROR`, so nothing prints while it works. Re-run with `-vv` or `-vvv`                                                                                                                      |
| A deletion appears where you expect a splicing change            | A genuine intron that the control sample did not support strongly enough is reported as a large deletion. Check whether the junction reached "standard" status in `refsplice`                                        |

## Checking a run against the bundled example

If your own data misbehaves and you want to confirm the installation itself is sound:

```
uv run pytest -m slow
```

This runs `refsplice`, `count` and `call` against the bundled example data and compares the
calls to a checked-in reference. It takes about two minutes. If it passes, the pipeline
works and the problem is in the inputs.
