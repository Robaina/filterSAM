#!/usr/bin/env python3
# coding: utf-8

"""
Tools to filter SAM/BAM files by percent identity and percent of matched sequence
"""

import re
from pathlib import Path

import pysam
from parallelbam.parallelbam import parallelizeBAMoperation

from filtersam.utils import terminalExecute, sam2bam


def extractSegmentsWithMDtag(sam_dir: Path, output_dir: Path=None,
                             suppress_output: bool = False) -> None:
    """
    Use samtools to filter out segments that do not have an MD tag
    """
    if output_dir is None:
        output_dir = Path(sam_dir.parent) / f"{sam_dir.stem}_only_md.{sam_dir.suffix}"
    samtools_command = f'samtools view -h -d MD {sam_dir} > {output_dir}'
    terminalExecute(samtools_command, suppress_output=suppress_output)
    
def sumMatchesAndMismatches(segment):
    """
    Get total matches/mismatches from CIGAR string (M field)
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    return sum(
        [value for (code, value) in segment.cigartuples if code == 0]
    )

def getNumberOfMatches(segment):
    """
    Get numnber of matches from alignment
    Do not consider insertion/deletion as mismatches
    """
    parsed_MD = segment.get_aligned_pairs(with_seq=True)
    return len([
        base for (read_pos, ref_pos, base) in parsed_MD 
        if ((base is not None and base.isupper()) and read_pos is not None)
    ])

def getQueryLength(segment):
    """
    Compute query length from CIGAR field corresponding to
    query sequence. 
    Following: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
    
    Cigar fields which 'consume sequence': M, I, S, =, X
    
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    codes = [0, 1, 4, 7, 8]
    return sum(
        [value for (code, value) in segment.cigartuples if code in codes]
    )

def percent_matched(segment):
    """
    Compute percentage of sequence that has been matched to reference
    """
    seq_length = getQueryLength(segment)
    n_matches = getNumberOfMatches(segment)
    return 100 * (n_matches / seq_length)
    
def percent_identity(segment):
    """
    Compute percent identity from MD tag of aligned segment.
    segment: pysam AlignedSegment object.
    """
    return 100 * (getNumberOfMatches(segment) / sumMatchesAndMismatches(segment))

def has_MD_tag(segment):
    return 'MD' in [tag for (tag, _) in segment.get_tags()]

def filterSAMbyIdentity(input_path: Path, output_path: Path = None,
                        identity_cutoff: float = 95.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent identity
    equal or above identity_cutoff value.
    """
    input_path = Path(input_path)
    file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = (f'{input_path.as_posix().split(file_ext)[0]}'
                       f'.identity_filtered_at_{identity_cutoff}{file_ext}')
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)                                      
    for segment in samfile:
        if (has_MD_tag(segment) and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)           
    filtered_sam.close()
    samfile.close()

def filterSAMbyPercentMatched(input_path: Path, output_path: Path = None,
                              matched_cutoff: float = 50.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent of matched
    based equal or higher than matched_cutoff. 
    
    Percent of matched bases is computed as the fraction of matches in
    the total query length.
    """
    file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = (f'{input_path.as_posix().split(file_ext)[0]}'
                       f'_matched_filtered_at_{matched_cutoff}{file_ext}') 
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_matched(segment) >= matched_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def filterSAM(input_path: Path, output_path: Path = None,
              filter_by: str = 'identity', cutoff: float = 95.0,
              n_processes: int = None) -> None:
    """
    Filter aligned segments in BAM or SAM file by percent identity or percent
    of matched sequence.
    """
    if filter_by == 'identity':
        filter_method = filterSAMbyIdentity
    elif filter_by == 'matched':
        filter_method = filterSAMbyPercentMatched
    else:
        raise ValueError('Invalid filter, available filters are: "identity" and "matched"')

    if not (cutoff >=0 and cutoff <= 100):
        raise ValueError('Cutoff value must be between 0 and 100.')
    
    if n_processes is None:
        filter_method(input_path, output_path, cutoff)
    else:
        if '.sam' in input_path.name:
            print('Converting sam file to bam for processing')
            input_path = Path(input_path.as_posix().replace('.sam', '.bam'))
            sam2bam(input_path)
        parallelizeBAMoperation(path_to_bam=input_path.as_posix(), callback=filter_method,
                                callback_additional_args=[cutoff],
                                n_processes=n_processes, output_path=output_path.as_posix())
