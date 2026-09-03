"""The per-cell stages produce identical output at every thread count.

``refsplice`` and ``count`` distribute their per-cell work over ``--threads``
worker processes. Their reductions depend on the order of the per-cell results,
so the parallel path must return results in input order rather than completion
order. This runs both stages serially and in parallel over the bundled example
and compares every file they write, byte for byte.

The figures are excluded: matplotlib stamps a wall-clock ``/CreationDate`` into
every PDF, so two consecutive serial runs already differ. ``call`` is excluded
because it is unaffected by this option -- it reads the YAML ``count`` writes,
which is written with sorted keys.

This is a slow integration test (~2 min).
"""

import shutil
import sys
from pathlib import Path

import pytest

from scptseq.main import main

EXAMPLE_DATA = Path(__file__).resolve().parents[1] / "examples" / "example_data"

RUN_NAME = "TX46_Prosalpha3"
GENE_NAME = "Prosalpha3"
TARGETS = EXAMPLE_DATA / "targets.yml"
GENOME = EXAMPLE_DATA / "Drosophila_melanogaster.BDGP6.28.dna.toplevel.chrm_2R.fa"
CONTROL_BAMS = EXAMPLE_DATA / "TX46_Prosalpha3_ctrl_bams"
PERTURBED_BAMS = EXAMPLE_DATA / "TX46_Prosalpha3_bams"


def _run_stages(tmp_path, monkeypatch, threads):
    """Run refsplice and count at the given thread count.

    Returns:
        `(results_dir, mutations_dir)` for that run. The perturbed BAMs are
        copied because `count` writes its per-cell files into that directory.
    """
    results_dir = tmp_path / f"results_{threads}"
    perturbed_bams = tmp_path / f"perturbed_bams_{threads}"
    shutil.copytree(PERTURBED_BAMS, perturbed_bams)

    results_arg = f"--results-dir={results_dir}"
    threads_arg = f"--threads={threads}"
    for argv in (
        ["refsplice", str(CONTROL_BAMS), RUN_NAME, GENE_NAME, str(TARGETS),
         results_arg, threads_arg],
        ["count", str(GENOME), str(perturbed_bams), RUN_NAME, GENE_NAME, str(TARGETS),
         results_arg, threads_arg],
    ):
        monkeypatch.setattr(sys, "argv", ["scptseq", *argv])
        main()

    return results_dir, perturbed_bams / "mutations"


def _yml_bytes(directory):
    """Map each ``*.yml`` file in `directory` to its contents, keyed by name."""
    return {path.name: path.read_bytes() for path in sorted(directory.glob("*.yml"))}


@pytest.mark.slow
def test_threads_do_not_change_output(tmp_path, monkeypatch):
    serial_results, serial_mutations = _run_stages(tmp_path, monkeypatch, 1)
    parallel_results, parallel_mutations = _run_stages(tmp_path, monkeypatch, 4)

    serial_yml = _yml_bytes(serial_results)
    parallel_yml = _yml_bytes(parallel_results)
    # Name the expected files, so that neither run silently writing nothing passes.
    assert set(serial_yml) == {
        f"{RUN_NAME}_splicing_junctions.yml",
        f"{RUN_NAME}_standard_splicing_junctions.yml",
        f"{RUN_NAME}_splicing_junction_frac_umis_per_cell.yml",
        f"{RUN_NAME}_splicing_junction_str.yml",
        f"{RUN_NAME}_stat_cntrs_before_and_after.yml",
    }, sorted(serial_yml)
    assert serial_yml.keys() == parallel_yml.keys()
    for name in serial_yml:
        assert serial_yml[name] == parallel_yml[name], name

    serial_cells = _yml_bytes(serial_mutations)
    parallel_cells = _yml_bytes(parallel_mutations)
    assert len(serial_cells) == 30, len(serial_cells)
    assert serial_cells.keys() == parallel_cells.keys()
    for name in serial_cells:
        assert serial_cells[name] == parallel_cells[name], name
