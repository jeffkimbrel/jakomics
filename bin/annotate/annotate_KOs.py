import sys
import os
import argparse
import time
from jakomics import utilities, colors, kegg
from multiprocessing import Pool, Manager
import uuid
import pandas as pd

print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

# utilities.test()
# colors.test()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

parser.add_argument('--profile',
                    help="kofamscan profile",
                    default='prokaryote.hal',
                    required=False)

parser.add_argument('--threads',
                    help="Threads",
                    default=8,
                    type=int,
                    required=False)


args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

manager = Manager()
counter = manager.list()

# FUNCTIONS ###################################################################


def main(file):
    output_path = args.out_dir
    profile = args.profile

    global counter
    print(
        f'Finished {len(counter)} of {len(file_list)}',
        end="\r",
        file=sys.stderr)

    if profile == 'prokaryote.hal':
        output_file = os.path.join(output_path, file.name + '.kofam94.txt')
    else:
        output_file = os.path.join(output_path, file.name + '.' +
                                   os.path.basename(profile) + '.txt')

    if os.path.exists(output_file):
        print(f'{colors.bcolors.RED}Skipping {output_file} because it already exists{colors.bcolors.END}')
    else:
        hits = kegg.run_kofam(file.file_path, args.profile)
        df = kegg.kofam_to_df(hits)
        df.to_csv(output_file, sep="\t", index=False)

    counter.append(file.short_name)
    print(
        f'Finished {len(counter)} of {len(file_list)}',
        end="\r",
        file=sys.stderr)
    # time.sleep(1)  # sleep so control c will cancel easier


## MAIN LOOP ###################################################################


if __name__ == "__main__":

    file_list = utilities.get_files(args.files, args.in_dir, ["faa"])

    if len(file_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid .faa files were found!{colors.bcolors.END}")

    pool = Pool(processes=args.threads)
    pool.map(main, file_list)
    pool.close()

    print()
