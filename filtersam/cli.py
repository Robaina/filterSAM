"""
Command line interface to filtersam
"""

import sys
import argparse
from pathlib import Path

from filtersam.filtersam import filterSAM


def main():

    def printLogo():
        print(
            (
            """
   __ _ _ _                                
  / _(_) | |                               
 | |_ _| | |_ ___ _ __ ___  __ _ _ __ ___  
 |  _| | | __/ _ \ '__/ __|/ _` | '_ ` _ \ 
 | | | | | ||  __/ |  \__ \ (_| | | | | | |
 |_| |_|_|\__\___|_|  |___/\__,_|_| |_| |_|
                                           
                                           
"""
            "Tools to filter SAM or BAM files\n"
            "Semidán Robaina Estévez (srobaina@ull.edu.es), 2022\n"
            " \n"
            )
            )
            
    printLogo()
  
    parser = argparse.ArgumentParser(
        prog='filtersam',
        description='Tools to filter SAM/BAM files by percent identity or percent of matched sequence',
        epilog="  \n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('bam', type=Path,
                        help = 'path to bam / sam file')
    parser.add_argument('-i', '--identity', metavar='', type=float,
                        dest='identity',
                        help='filter by given percent identity')
    parser.add_argument('-m', '--matched', metavar='', type=float,
                        dest='matched',
                        help='filter by given percent of matched sequence')
    parser.add_argument('-p', '--processes', metavar='', type=int,
                        dest='processes',
                        help='number of processes for parallelization')
    parser.add_argument('-o', '--output', metavar='', type=Path,
                        default=None, dest='out',
                        help='path to output file')
    
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    bam = Path(args.bam).absolute()

    if not bam.is_file():
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

    
    
