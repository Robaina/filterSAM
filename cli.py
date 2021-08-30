"""
Command line interface to filtersam
"""

import os
import sys
import argparse
from filtersam.filtersam import filterSAM


def main():
  
    parser = argparse.ArgumentParser(
        prog ='filtersam',
        description ='Tools to filter SAM/BAM files by percent identity or percent of matched sequence',
        epilog = 'Developed by Semidán Robaina Estévez (srobaina@ull.edu.es)'
    )

    parser.add_argument('bam', type = str,
                        help = 'path to bam / sam file')
    parser.add_argument('-i', '--identity', metavar = '', type = float,
                        dest ='identity',
                        help = 'filter by given percent identity')
    parser.add_argument('-m', '--matched', metavar = '', type = float,
                        dest ='matched',
                        help = 'filter by given percent of matched sequence')
    parser.add_argument('-p', '--processes', metavar = '', type = int,
                        dest = 'processes',
                        help = 'number of processes for parallelization')
    parser.add_argument('-o', '--output', metavar = '', type = str,
                        default = None, dest = 'out',
                        help = 'path to output file')

    args = parser.parse_args()
    bam = os.path.abspath(args.bam)

    if not os.path.isfile(bam):
        print('Specified bam file does not exist')
        sys.exit()

    if args.identity is not None:
        print(f'Filtering by percent identity at {args.identity}%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='identity',
                  cutoff=args.identity, n_processes=args.processes)
    if args.matched is not None:
        print(f'Filtering by percent of matched sequence at {args.matched}%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='matched',
                  cutoff=args.matched, n_processes=args.processes)
    if args.identity is None and args.matched is None:
        print(f'Defaulting to filtering by percent identity at 95%')
        filterSAM(input_path=bam, output_path=args.out, filter_by='identity',
                  cutoff=95.0, n_processes=args.processes)
        
        
        
main()   

    
    
