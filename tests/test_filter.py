"""
Tests for the single-file filtering functions:
filterSAMbyIdentity and filterSAMbyPercentMatched.
"""

from pathlib import Path

import pytest

from filtersam.filtersam import filterSAMbyIdentity, filterSAMbyPercentMatched

from conftest import IDENTITY_KEPT, MATCHED_KEPT, detect_format, read_segment_names


@pytest.mark.parametrize("cutoff", sorted(IDENTITY_KEPT))
def test_filter_by_identity_keeps_expected(sam_path, tmp_path, cutoff):
    out = tmp_path / "out.sam"
    filterSAMbyIdentity(sam_path, out, identity_cutoff=cutoff)
    assert read_segment_names(out) == IDENTITY_KEPT[cutoff]


@pytest.mark.parametrize("cutoff", sorted(MATCHED_KEPT))
def test_filter_by_matched_keeps_expected(sam_path, tmp_path, cutoff):
    out = tmp_path / "out.sam"
    filterSAMbyPercentMatched(sam_path, out, matched_cutoff=cutoff)
    assert read_segment_names(out) == MATCHED_KEPT[cutoff]


def test_segments_without_md_tag_are_always_dropped(sam_path, tmp_path):
    # Cutoff 0 keeps everything that has an MD tag, but never the MD-less read.
    out = tmp_path / "out.sam"
    filterSAMbyIdentity(sam_path, out, identity_cutoff=0.0)
    assert "read_nomd" not in read_segment_names(out)


def test_works_on_bam_input(bam_path, tmp_path):
    out = tmp_path / "out.bam"
    filterSAMbyIdentity(bam_path, out, identity_cutoff=95.0)
    assert read_segment_names(out) == IDENTITY_KEPT[95.0]


# --- Output format selection (regression guard for the single-process BAM fix) ---

def test_output_sam_is_text(sam_path, tmp_path):
    out = tmp_path / "out.sam"
    filterSAMbyIdentity(sam_path, out, identity_cutoff=95.0)
    assert detect_format(out) == "sam"


def test_output_bam_is_binary(sam_path, tmp_path):
    out = tmp_path / "out.bam"
    filterSAMbyIdentity(sam_path, out, identity_cutoff=95.0)
    assert detect_format(out) == "bam"


def test_output_format_follows_output_extension_not_input(bam_path, tmp_path):
    # BAM input but a .sam output request must yield text SAM.
    out = tmp_path / "out.sam"
    filterSAMbyIdentity(bam_path, out, identity_cutoff=95.0)
    assert detect_format(out) == "sam"
    assert read_segment_names(out) == IDENTITY_KEPT[95.0]


# --- Default output path naming (regression guard for the suffix fix) ---

def test_default_output_path_naming(tmp_path):
    # A filename containing 'sam' before the extension used to break the old
    # regex-based extension detection; the suffix-based logic handles it.
    src = tmp_path / "mysample.bam"
    # Build a tiny BAM from the SAM fixture text.
    from conftest import SAM_TEXT
    import pysam
    sam = tmp_path / "_seed.sam"
    sam.write_text(SAM_TEXT)
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(sam), "r") as s, \
            pysam.AlignmentFile(str(src), "wb", template=s) as d:
        for seg in s:
            d.write(seg)
    pysam.set_verbosity(save)

    filterSAMbyIdentity(src, identity_cutoff=95.0)

    expected = tmp_path / "mysample.identity_filtered_at_95.0.bam"
    assert expected.is_file()
    assert detect_format(expected) == "bam"
    assert read_segment_names(expected) == IDENTITY_KEPT[95.0]
