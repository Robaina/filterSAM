"""
Command line interface to filtersam
"""

import os
import sys
import argparse
from filtersam.filtersam import filterSAM

  
parser = argparse.ArgumentParser(prog ='filtersam',
                                    description ='Tools to filter SAM/BAM files by percent identity or percent of matched sequence',
                                    epilog = 'Developed by Semidán Robaina Estévez (srobaina@ull.edu.es)')

parser.add_argument('bam', type = str,
                    help = 'path to bam / sam file')
parser.add_argument('-c', '--cutoff', metavar = '', type = float,
                    default = 95.0, dest = 'cutoff',
                    help ='percent cutoff value')
parser.add_argument('-f', '--filter', metavar = None, type = str,
                    default = 'identity', dest ='filter',
                    choices = ['identity', 'matched'],
                    help = 'filter to be appplied')
parser.add_argument('-p', '--processes', metavar = '', type = int,
                    dest = 'processes',
                    help = 'number of processes for parallelization')

args = parser.parse_args()
if not os.path.isdir(args.bam):
    print('The path specified does not exist')
    sys.exit()

filterSAM(input_path=os.path.abspath(args.bam), output_path=args.out, filter_by=args.filter,
          cutoff=args.cutoff, n_processes=args.processes)
    