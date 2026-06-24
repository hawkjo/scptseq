"""End-to-end test: the bundled example pipeline reproduces its reference output.

Runs the three analysis stages (``refsplice`` -> ``count`` -> ``call``) on the
data in ``examples/example_data`` exactly as the ``examples/*.sh`` scripts do,
and asserts that the produced per-cell mutation calls match the checked-in,
author-validated reference file ``ref_TX46_Prosalpha3_mutation_statuses.txt``.

This is a slow integration test (~2 min). 
"""

import shutil
import sys
from pathlib import Path

import pytest

from scptseq.main import main

# examples/example_data lives two levels up from this test file.
EXAMPLE_DATA = Path(__file__).resolve().parents[1] / "examples" / "example_data"

RUN_NAME = "TX46_Prosalpha3"
GENE_NAME = "Prosalpha3"
TARGETS = EXAMPLE_DATA / "targets.yml"
GENOME = EXAMPLE_DATA / "Drosophila_melanogaster.BDGP6.28.dna.toplevel.chrm_2R.fa"
CONTROL_BAMS = EXAMPLE_DATA / "TX46_Prosalpha3_ctrl_bams"
PERTURBED_BAMS = EXAMPLE_DATA / "TX46_Prosalpha3_bams"
REFERENCE = EXAMPLE_DATA / "ref_TX46_Prosalpha3_mutation_statuses.txt"


def _parse_mutation_statuses(fpath):
    """Parse a *_mutation_statuses.txt TSV into {barcode: (maternal, paternal)}.

    The file has a ``barcode\tmaternal\tpaternal`` header followed by one row
    per called cell. Returning a dict makes equality checks order-independent
    and gives a readable diff on failure.
    """
    statuses = {}
    with open(fpath) as f:
        header = next(f).rstrip("\n").split("\t")
        assert header == ["barcode", "maternal", "paternal"], header
        for line in f:
            barcode, maternal, paternal = line.rstrip("\n").split("\t")
            statuses[barcode] = (maternal, paternal)
    return statuses


def _run_stage(monkeypatch, argv):
    """Invoke the scptseq CLI in-process with the given argv (sans program name)."""
    monkeypatch.setattr(sys, "argv", ["scptseq", *argv])
    main()


@pytest.mark.slow
def test_example_pipeline_reproduces_reference(tmp_path, monkeypatch):
    results_dir = tmp_path / "results"

    # `count` writes a mutations/ subdir into the perturbed BAM dir, so work on a
    # copy to keep the checked-in example data pristine. The genome and control
    # BAMs are only read, so they are used in place.
    perturbed_bams = tmp_path / "perturbed_bams"
    shutil.copytree(PERTURBED_BAMS, perturbed_bams)

    results_arg = f"--results-dir={results_dir}"
    _run_stage(
        monkeypatch,
        [
            "refsplice",
            str(CONTROL_BAMS),
            RUN_NAME,
            GENE_NAME,
            str(TARGETS),
            results_arg,
        ],
    )
    _run_stage(
        monkeypatch,
        [
            "count",
            str(GENOME),
            str(perturbed_bams),
            RUN_NAME,
            GENE_NAME,
            str(TARGETS),
            results_arg,
        ],
    )
    _run_stage(monkeypatch, ["call", RUN_NAME, GENE_NAME, str(TARGETS), results_arg])

    # Key intermediate outputs should exist and be non-empty (structural smoke
    # check; we do not compare their contents).
    for name in (
        f"{RUN_NAME}_standard_splicing_junctions.yml",
        f"{RUN_NAME}_splicing_junction_str.yml",
        f"{RUN_NAME}_stat_cntrs_before_and_after.yml",
    ):
        intermediate = results_dir / name
        assert intermediate.is_file() and intermediate.stat().st_size > 0, name

    produced = results_dir / f"{RUN_NAME}_mutation_statuses.txt"
    assert produced.is_file(), produced

    assert _parse_mutation_statuses(produced) == _parse_mutation_statuses(REFERENCE)
