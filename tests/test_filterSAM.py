"""
Tests for the filterSAM dispatcher: argument validation and the
single-vs-parallel routing logic.
"""

import pytest

from filtersam import filtersam as fs
from filtersam.filtersam import filterSAM

from conftest import IDENTITY_KEPT, read_segment_names


def test_invalid_filter_by_raises(sam_path, tmp_path):
    with pytest.raises(ValueError):
        filterSAM(sam_path, tmp_path / "out.sam", filter_by="nonsense", cutoff=95.0)


@pytest.mark.parametrize("cutoff", [-1.0, 100.1, 1000.0])
def test_out_of_range_cutoff_raises(sam_path, tmp_path, cutoff):
    with pytest.raises(ValueError):
        filterSAM(sam_path, tmp_path / "out.sam", filter_by="identity", cutoff=cutoff)


def test_none_processes_uses_single_path(bam_path, tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(fs, "parallelizeBAMoperation",
                        lambda *a, **k: calls.append((a, k)))
    out = tmp_path / "out.bam"
    filterSAM(bam_path, out, filter_by="identity", cutoff=95.0, n_processes=None)
    assert calls == []
    assert read_segment_names(out) == IDENTITY_KEPT[95.0]


def test_single_process_does_not_split(bam_path, tmp_path, monkeypatch):
    # Regression guard for #3 / PR #6: `-p 1` must take the direct path and
    # never invoke the (expensive) parallel splitting machinery.
    calls = []
    monkeypatch.setattr(fs, "parallelizeBAMoperation",
                        lambda *a, **k: calls.append((a, k)))
    out = tmp_path / "out.bam"
    filterSAM(bam_path, out, filter_by="identity", cutoff=95.0, n_processes=1)
    assert calls == [], "n_processes=1 should not call parallelizeBAMoperation"
    assert read_segment_names(out) == IDENTITY_KEPT[95.0]


def test_multiple_processes_use_parallel_path(bam_path, tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(fs, "parallelizeBAMoperation",
                        lambda *a, **k: calls.append((a, k)))
    out = tmp_path / "out.bam"
    filterSAM(bam_path, out, filter_by="identity", cutoff=95.0, n_processes=2)
    assert len(calls) == 1, "n_processes>1 should call parallelizeBAMoperation"


def test_parallel_run_matches_serial(bam_path, tmp_path):
    # End-to-end parallel run (uses samtools-backed splitting from parallelbam).
    out = tmp_path / "out_parallel.bam"
    filterSAM(bam_path, out, filter_by="identity", cutoff=70.0, n_processes=2)
    assert read_segment_names(out) == IDENTITY_KEPT[70.0]
