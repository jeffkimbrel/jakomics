import argparse
import os
import sys
from natsort import natsorted
from multiprocessing import Manager, Pool

from jakomics.utilities import blast, colors

print(f'{colors.bcolors.GREEN}THIS IS STILL A WORK IN PROGRESS!{colors.bcolors.END}')

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('--in_dir', help="Directory with MAGs", required=True)
parser.add_argument('-db', '--database', help="Amplicon FASTA file", required=True)

args = parser.parse_args()
args.in_dir = os.path.abspath(args.in_dir) + '/'

manager = Manager()
shared_list = manager.list()

# FUNCTIONS ###################################################################


def get_files():
    files = []
    dirs = os.listdir(args.in_dir)
    for fileName in dirs:
        if fileName.endswith('fa'):
            files.append(args.in_dir + fileName)

    return natsorted(files)


def blast_call(file):
    global shared_list

    results = blast.do_blast(type="nucl",
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

    blast.make_blast_db("nucl", args.database)
    file_list = get_files()

    pool = Pool(processes=8)
    pool.map(blast_call, natsorted(file_list))
    pool.close()

    for query in shared_list:
        for subject in query:
            subject.print_full_result()
