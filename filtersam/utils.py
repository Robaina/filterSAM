"""
Functions for general purposes
"""

import os
from pathlib import Path


def terminalExecute(command_str: str, suppress_output=False) -> None:
    """
    Execute given command in terminal through Python
    """
    if suppress_output:
        suppress_code = '>/dev/null 2>&1'
        command_str = f'{command_str} {suppress_code}'
    os.system(command_str)

def sam2bam(sam_file: Path, output_dir: Path = None) -> None:
    """
    Convert sam file to bam
    """
    if output_dir is None:
        output_dir = Path(sam_file.as_posix().replace('.sam', '.bam'))
    is_sam = '.sam' in sam_file.name
    if is_sam:
        terminalExecute(f'samtools view {sam_file} -S -b -o {output_dir}')