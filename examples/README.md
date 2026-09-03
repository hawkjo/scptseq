# Examples

Run these from inside this directory — the paths in each script are relative to it.

```
cd examples
./splitfastqs_example.sh
./refsplice_example.sh
./count_example.sh
./call_example.sh
```

`splitfastqs_example.sh` splits the bundled sample FASTQ into one file per barcode. The
other three scripts do not consume its output. The aligned perturbed BAMs are already
checked in under `example_data/TX46_Prosalpha3_bams/`, so the analysis runs without
repeating the alignment step. All three share `--results-dir=results` and the run name
`TX46_Prosalpha3`, and write into `results/`.

Two notes:

- `count_example.sh` writes a `mutations/` directory inside
  `example_data/TX46_Prosalpha3_bams/`, which is part of the tree. Delete it after
  a run.
- The final `results/TX46_Prosalpha3_mutation_statuses.txt` should match the checked-in
  `example_data/ref_TX46_Prosalpha3_mutation_statuses.txt`. The reference file's rows are in a different order, so compare with `sort` on both sides.

For a walkthrough of what each stage does and how to adapt this to your own data, see
[../docs/tutorial.md](../docs/tutorial.md).
