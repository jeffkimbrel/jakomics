import gzip
import argparse
import os
from natsort import natsorted
from multiprocessing import Pool
from jakomics import utilities, colors

from Bio.SeqIO.QualityIO import FastqGeneralIterator

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='get basic stats on a directory (--in_dir) or list (-f) of fastq.gz files')

parser.add_argument('-f', '--files',
                    help="fastq.gz files",
                    nargs='*',
                    required=False,
                    default="")

parser.add_argument('--in_dir',
                    help="Directory with fastq.gz files",
                    required=False,
                    default="")

parser.add_argument('-t', '--threads',
                    help="Threads for multiprocessing",
                    required=False,
                    default=8,
                    type=int)

args = parser.parse_args()

# FUNCTIONS ###################################################################

def main(file):
    count = 0
    bp = 0

    with gzip.open(file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            count += 1
            bp += len(seq)

    print(f'{file}\t{count}\t{bp}')


# MAIN ########################################################################

if __name__ == "__main__":

    file_list = utilities.get_file_list(args.files, ending=['fastq.gz'])
    file_list = utilities.get_directory_file_list(
        args.in_dir, ending=['fastq.gz'], file_list=file_list)

    pool = Pool(processes=args.threads)
    pool.map(main, natsorted(file_list))
    pool.close()
