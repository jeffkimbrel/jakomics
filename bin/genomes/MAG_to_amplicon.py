import argparse
import os
import sys
from natsort import natsorted
from multiprocessing import Manager, Pool

from jakomics import blast, colors, utilities

print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('--in_dir', help="Directory with MAGs", required=False, default="")
parser.add_argument('-f', '--files',
                    help="Paths to individual genome files",
                    nargs='*',
                    required=False,
                    default=None)
parser.add_argument('-db', '--database', help="Amplicon FASTA file", required=True)
parser.add_argument('-s', '--suffix', help="File suffixes to include",
                    required=False, nargs="*", default=["fa", "fna", "fasta"])

args = parser.parse_args()

manager = Manager()
shared_list = manager.list()

# FUNCTIONS ###################################################################


def blast_call(file):
    global shared_list
    global bar

    results = blast.run_blast(type="nucl",
                              q=file,
                              db=args.database,
                              e=1e-50)

    passed = []

    for result in sorted(results.keys()):
        for hit in results[result]:
            if hit.id >= 99:
                passed.append(hit)

    shared_list.append(passed)
    print(f'> Processed {len(shared_list)} of {len(file_list)} genomes...',
          end="\r", file=sys.stderr)

### MAIN ######################################################################


if __name__ == "__main__":

    try:
        blast.make_blast_db("nucl", args.database)
    except:
        sys.exit(
            f"{colors.bcolors.RED}Error: Blast database file not found, or not a fasta file{colors.bcolors.END}")

    file_list = utilities.get_file_list(args.files, args.suffix)
    file_list = utilities.get_directory_file_list(args.in_dir, args.suffix, file_list)

    if len(file_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid genomes were found ({str(args.suffix)}){colors.bcolors.END}")

    pool = Pool(processes=8)
    pool.map(blast_call, natsorted(file_list))
    pool.close()

    for query in shared_list:
        for subject in query:
            subject.print_full_result()
