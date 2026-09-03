"""Unit tests for ``TargetInfo``: target ordering, validation, and cut-site distance.

These cover the parts of the target-file contract whose failures would be silent rather
than loud — a mis-ordered target list mislabels every per-target figure, and a
cut-site distance measured against the wrong target flips the ``max_dist`` filter that
decides which mutations join a signature.
"""

import pytest
import yaml

from scptseq.misc import TargetInfo

GENE = "TestGene"


def write_targets(tmp_path, cutsite_given_target, seg_sites=(), seg_bases=()):
    """Write a one-gene target file and return its path.

    Written with ``yaml.dump``, which sorts keys — so the targets land in the file in the
    same alphabetical order that real, machine-written target files have.

    Args:
        tmp_path: Directory to write into.
        cutsite_given_target: Maps target name to cut-site position.
        seg_sites: Segregating-site positions.
        seg_bases: ``[maternal, paternal]`` base pair per segregating site.

    Returns:
        Path to the written file.
    """
    target_info = {
        GENE: {
            "chrm": "2R",
            "start": 0,
            "end": 100000,
            "targets": {
                tname: {"cutsite": cutsite} for tname, cutsite in cutsite_given_target.items()
            },
            "seg_info": {
                "seg_sites": list(seg_sites),
                "seg_bases": [list(pair) for pair in seg_bases],
                "skip_sites": [],
            },
        }
    }
    fpath = tmp_path / "targets.yml"
    with open(fpath, "w") as out:
        yaml.dump(target_info, out)
    return fpath


def test_target_names_ordered_numerically_not_alphabetically(tmp_path):
    # The cut-site value encodes the target number, so the expectation below does not
    # depend on how the implementation derives the order.
    fpath = write_targets(tmp_path, {"t1": 1000, "t2": 2000, "t10": 10000})

    # Guard the guard: if the file stopped carrying the hazardous order, this test would
    # silently stop testing anything.
    assert list(yaml.safe_load(open(fpath))[GENE]["targets"]) == ["t1", "t10", "t2"]

    target_info = TargetInfo(fpath, GENE)
    assert target_info.target_names == ["t1", "t2", "t10"]
    assert target_info.cutsites == [1000, 2000, 10000]


def test_seg_site_order_uses_t1_not_the_first_target_in_the_file(tmp_path):
    # Written as raw text so t2 genuinely precedes t1 in the file.
    fpath = tmp_path / "targets.yml"
    fpath.write_text(
        f"""
{GENE}:
  chrm: 2R
  start: 0
  end: 100000
  targets:
    t2:
      cutsite: 100
    t1:
      cutsite: 900
  seg_info:
    seg_sites: [200, 800]
    seg_bases: [[G, C], [T, A]]
    skip_sites: []
"""
    )
    target_info = TargetInfo(fpath, GENE)

    assert target_info.t1_cutsite == 900
    # Distances from t1 are 700 and 100. Ordering against the first target in the file
    # (t2, at 100) would instead give [200, 800].
    assert target_info.sorted_seg_sites == [800, 200]


def test_single_target_gene_is_accepted(tmp_path):
    target_info = TargetInfo(write_targets(tmp_path, {"t1": 500}), GENE)

    assert target_info.target_names == ["t1"]
    assert target_info.cutsites == [500]
    assert target_info.t1_cutsite == 500


@pytest.mark.parametrize(
    "cutsite_given_target, expected_in_message",
    [
        ({"t2": 100, "t3": 200}, "t2"),  # t1 absent
        ({"t1": 100, "t2a": 200}, "t2a"),  # not t followed by digits
        ({"t0": 100, "t1": 200}, "t0"),  # t0 would sort ahead of the reference target
        ({"t1": 100, "t01": 200}, "t01"),  # ties with t1 on the numeric sort key
        ({"t1": 100, "notes": "free text"}, "notes"),  # targets takes no other keys
    ],
)
def test_malformed_target_names_raise_value_error(
    tmp_path, cutsite_given_target, expected_in_message
):
    fpath = write_targets(tmp_path, cutsite_given_target)

    with pytest.raises(ValueError) as excinfo:
        TargetInfo(fpath, GENE)

    # The gene name is the actionable part when a file holds many genes.
    assert GENE in str(excinfo.value)
    assert expected_in_message in str(excinfo.value)


def test_mut_min_dist_to_cutsite_uses_nearest_of_all_targets(tmp_path):
    target_info = TargetInfo(write_targets(tmp_path, {"t1": 1000, "t2": 2000, "t3": 3000}), GENE)

    # A two-cut-site implementation returns 996 for A2996T, which silently flips the
    # `max_dist <= 10` filter that decides whether a mutation is target-adjacent.
    assert target_info.mut_min_dist_to_cutsite("A2996T") == 4
    assert target_info.mut_min_dist_to_cutsite("I2005ACGT") == 5
    assert target_info.mut_min_dist_to_cutsite("D1010-1020") == 10
    assert target_info.mut_min_dist_to_cutsite("D990-1010") == 0  # span crosses the cut site
    assert target_info.mut_min_dist_to_cutsite("J1500_1600") == 400
    # Offset junctions (J{start}_{end}>{d}_{d}) are deliberately untested here: the
    # offset adjustment in mut_min_dist_to_cutsite is guarded on the deltas being an int
    # while parse_mutation returns a tuple, so it never fires. Correcting that changes
    # calls and is scheduled separately in ROADMAP.md.


def test_bundled_example_target_file_is_unchanged(example_targets_file):
    """Guard the frozen two-target values the paper's results were produced with."""
    prosalpha3 = TargetInfo(example_targets_file, "Prosalpha3")
    assert prosalpha3.target_names == ["t1", "t2"]
    assert prosalpha3.cutsites == [20995489, 20994988]
    assert prosalpha3.sorted_seg_sites == [20995104, 20994903]

    # In both genes t1's cut site is the larger of the two, so this also catches an
    # accidental sort by position rather than by target number.
    lamc = TargetInfo(example_targets_file, "LamC")
    assert lamc.target_names == ["t1", "t2"]
    assert lamc.cutsites == [14575540, 14572821]
