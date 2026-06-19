"""
Shared fixtures for the filtersam test suite.

A small, hand-crafted SAM file with segments whose percent-identity and
percent-matched values are known exactly is used throughout. The values were
verified against pysam's own ``get_aligned_pairs``/CIGAR parsing:

    name            %identity   %matched   MD tag
    read_perfect       100.0      100.0     MD:Z:10
    read_90id           90.0       90.0     MD:Z:5A4
    read_softclip      100.0       80.0     MD:Z:8     (CIGAR 2S8M)
    read_70id           70.0       70.0     MD:Z:1A2C2A2
    read_nomd           ----       ----     (no MD tag, always dropped)
"""

import pysam
import pytest

SAM_TEXT = """@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:ref\tLN:100
read_perfect\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tMD:Z:10
read_90id\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTAGGTAC\t*\tMD:Z:5A4
read_softclip\t0\tref\t1\t60\t2S8M\t*\t0\t0\tTTACGTACGT\t*\tMD:Z:8
read_70id\t0\tref\t1\t60\t10M\t*\t0\t0\tAGGTGGTGAC\t*\tMD:Z:1A2C2A2
read_nomd\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0
"""

# Expected segments retained for a given filter and cutoff.
IDENTITY_KEPT = {
    95.0: {"read_perfect", "read_softclip"},
    90.0: {"read_perfect", "read_90id", "read_softclip"},
    70.0: {"read_perfect", "read_90id", "read_softclip", "read_70id"},
}
MATCHED_KEPT = {
    100.0: {"read_perfect"},
    85.0: {"read_perfect", "read_90id"},
    50.0: {"read_perfect", "read_90id", "read_softclip", "read_70id"},
}


def _write_sam(path):
    path.write_text(SAM_TEXT)
    return path


def read_segment_names(path):
    """Return the set of query names in a SAM/BAM file (format auto-detected)."""
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(path), "r") as handle:
        names = {seg.query_name for seg in handle}
    pysam.set_verbosity(save)
    return names


def detect_format(path):
    """Return 'bam' if the file is BGZF/BAM-compressed, 'sam' if plain text."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    return "bam" if magic == b"\x1f\x8b" else "sam"


@pytest.fixture
def sam_path(tmp_path):
    """A plain-text SAM fixture file."""
    return _write_sam(tmp_path / "sample.sam")


@pytest.fixture
def bam_path(tmp_path):
    """A BAM fixture file built from the same records as ``sam_path``."""
    sam = _write_sam(tmp_path / "_src.sam")
    bam = tmp_path / "sample.bam"
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(sam), "r") as src:
        with pysam.AlignmentFile(str(bam), "wb", template=src) as dst:
            for seg in src:
                dst.write(seg)
    pysam.set_verbosity(save)
    return bam


@pytest.fixture
def segments(sam_path):
    """The parsed AlignedSegment objects, keyed by query name."""
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(str(sam_path), "r") as handle:
        segs = {seg.query_name: seg for seg in handle}
    pysam.set_verbosity(save)
    return segs
