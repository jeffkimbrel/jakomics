import argparse
import os
from natsort import natsorted

from utilities import blast

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('--in_dir', help="Directory with MAGs", required=True)
parser.add_argument('-db', '--database', help="Amplicon FASTA file", required=True)

args = parser.parse_args()
args.in_dir = os.path.abspath(args.in_dir) + '/'

# MAIN ########################################################################

blast.make_blast_db("nucl", args.database)

for file in natsorted(os.listdir(args.in_dir)):
    if file.endswith('.fa'):
        # print(file)
        results = blast.do_blast(type="nucl",
                                 q=args.in_dir+file,
                                 db=args.database,
                                 e=1e-50)

        for result in sorted(results.keys()):
            for hit in results[result]:
                if hit.id >= 99:
                    hit.print_full_result()
