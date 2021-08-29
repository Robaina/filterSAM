#!/usr/bin/env python3
# coding: utf-8

"""
Tools to filter SAM/BAM files by percent identity and percent of matched sequence
"""

import pysam
import os
import sys
import re
import time
import numpy as np
from .utils import terminalExecute


def extractSegmentsWithMDtag(sam_dir: str, output_dir: str=None,
                             suppress_output: bool = False) -> None:
    """
    Use samtools to filter out segments that do not have an MD tag
    """
    if output_dir is None:
        basename, ext = os.path.basename(sam_dir).split('.')
        output_dir = os.path.join(os.path.dirname(sam_dir), f'{basename}_only_md.{ext}')
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

def filterSAMbyIdentity(input_path: str, output_path: str = None,
                        identity_cutoff: float = 95.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent identity
    equal or above identity_cutoff value.
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'.identity_filtered_at_{identity_cutoff}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def filterSAMbyPercentMatched(input_path: str, output_path: str = None,
                              matched_cutoff: float = 50.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent of matched
    based equal or higher than matched_cutoff. 
    
    Percent of matched bases is computed as the fraction of matches in
    the total query length.
    """
    file_ext = re.search('.(s|b)am', input_path).group()
    if output_path is None:
        output_path = (f'{input_path.split(file_ext)[0]}'
                       f'_matched_filtered_at_{matched_cutoff}{file_ext}') 
    
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_matched(segment) >= matched_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def filterSAM(input_path: str, output_path: str = None,
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
        pass



# Run script: python3 filter_by_identity.py input.bam identity_cutoff [output_path]
if __name__ == "__main__":
    """
    Filter records in SAM/BAM file by given percent identity.
    Usage: 
    python3 filter_by_identity.py input.{bam|sam} identity_cutoff [output_path]
    """
    
    input_file = sys.argv[1]
    if (len(sys.argv) > 2):
        identity_cutoff = int(sys.argv[2])
    else:
        identity_cutoff = 95
    if (len(sys.argv) > 3):
        output_path = int(sys.argv[3])
    else:
        output_path = None
    
    start = time.time()
    filterSAMbyIdentity(input_file, identity_cutoff, output_path)
    # out = computeSAMstatistics(input_file, identity_cutoff)
    end = time.time()
    
    print(f'Execution time: {end - start}')