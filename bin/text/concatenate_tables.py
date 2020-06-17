import os
import sys
import argparse
import pandas as pd
from jakomics import utilities, colors

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with table files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual table files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('-o', '--out',
                    help="Output file",
                    required=True)

parser.add_argument('-s', '--sep',
                    help="Column separator",
                    default="\t",
                    required=False)

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == '__main__':

    file_list = utilities.get_files(args.files, args.in_dir, ["txt", "csv", "tsv"])

    merged = pd.DataFrame()

    for file in file_list:
        df = pd.read_csv(file.file_path, sep=args.sep)
        df.insert(loc=0, column="FILE", value=file.name, allow_duplicates=True)
        merged = pd.concat([merged, df], axis=0, sort=False)

    merged.to_csv(args.out, sep=args.sep, index=False)
