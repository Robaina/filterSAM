"""
Functions for general purposes
"""

import os
import subprocess
import numpy as np


def terminalExecute(command_str: str, suppress_output=False) -> None:
    """
    Execute given command in terminal through Python
    """
    if suppress_output:
        suppress_code = '>/dev/null 2>&1'
        command_str = f'{command_str} {suppress_code}'
    os.system(command_str)

def sam2bam(sam_file: str, output_dir: str = None) -> None:
    """
    Convert sam file to bam
    """
    if output_dir is None:
        output_dir = sam_file.replace('.sam', '.bam')
    is_sam = '.sam' in sam_file
    if is_sam:
        terminalExecute(f'samtools view {sam_file} -S -b -o {output_dir}')
    
# def deleteTemporaryFiles(dir_path: str) -> None:
#     """
#     Remove temporary SAM files from directory
#     """
#     for fname in os.listdir(dir_path):
#         os.remove(os.path.join(dir_path, fname))

# def getFastqPairedFiles(data_dir: str, pattern: tuple=('_1.', '_2.')) -> dict:
#     """
#     Group paired-end fastq files by condition
#     """
#     fnames = os.listdir(data_dir)
#     conditions = np.unique([fname.split(pattern[0])[0].split(pattern[1])[0] 
#                             for fname in fnames]).tolist()
#     return {
#         condition: np.sort(
#             [fname for fname in fnames if condition in fname]
#         ).tolist() 
#         for condition in conditions
#     }
    
# def sortSAMbyName(sam_file: str, output_dir=None, suppress_output=False) -> None:
#     """
#     Sort SAM entries by name. Required by htseq-count to process paired-end data.
#     """
#     if output_dir is None:
#         output_dir = f'{sam_file.split(".sam")[0]}_sorted.sam'
#     samtools_command = f'samtools sort -n -O sam {sam_file} > {output_dir}'
#     terminalExecute(samtools_command, suppress_output=suppress_output)  
