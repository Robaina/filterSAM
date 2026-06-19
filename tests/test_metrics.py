"""
Unit tests for the per-segment metric helpers in filtersam.filtersam.
"""

import pytest

from filtersam import filtersam as fs


def test_has_md_tag(segments):
    assert fs.has_MD_tag(segments["read_perfect"])
    assert not fs.has_MD_tag(segments["read_nomd"])


def test_sum_matches_and_mismatches(segments):
    # Sum of CIGAR M lengths.
    assert fs.sumMatchesAndMismatches(segments["read_perfect"]) == 10
    assert fs.sumMatchesAndMismatches(segments["read_90id"]) == 10
    # Soft-clipped bases (2S) do not count towards M.
    assert fs.sumMatchesAndMismatches(segments["read_softclip"]) == 8


def test_get_number_of_matches(segments):
    assert fs.getNumberOfMatches(segments["read_perfect"]) == 10
    assert fs.getNumberOfMatches(segments["read_90id"]) == 9
    assert fs.getNumberOfMatches(segments["read_softclip"]) == 8
    assert fs.getNumberOfMatches(segments["read_70id"]) == 7


def test_get_query_length(segments):
    # M + I + S + = + X consume query; soft clip is included.
    assert fs.getQueryLength(segments["read_perfect"]) == 10
    assert fs.getQueryLength(segments["read_softclip"]) == 10


@pytest.mark.parametrize(
    "name,expected",
    [
        ("read_perfect", 100.0),
        ("read_90id", 90.0),
        ("read_softclip", 100.0),
        ("read_70id", 70.0),
    ],
)
def test_percent_identity(segments, name, expected):
    assert fs.percent_identity(segments[name]) == pytest.approx(expected)


@pytest.mark.parametrize(
    "name,expected",
    [
        ("read_perfect", 100.0),
        ("read_90id", 90.0),
        ("read_softclip", 80.0),
        ("read_70id", 70.0),
    ],
)
def test_percent_matched(segments, name, expected):
    assert fs.percent_matched(segments[name]) == pytest.approx(expected)
