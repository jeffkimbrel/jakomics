import argparse
import sys
import pandas as pd
from jakomics import utilities, colors

# FUNCTIONS ###################################################################


def read_file(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.set_index("ID")
    return df


def main(file, min_abundance, min_columns, out_file):
    print(f'Reading in {file}')

    df = read_file(file)
    row_drop = []
    keeper_count = 0
    drop_count = 0

    for index, row in df.iterrows():
        # print(row)
        threshold_count = row[row >= min_abundance].count()
        if threshold_count < min_columns:
            row_drop.append(index)
            drop_count += 1
        else:
            keeper_count += 1
        print(f'Filtering results: {colors.bcolors.GREEN}Good rows = {keeper_count}{colors.bcolors.END}, {colors.bcolors.RED}Bad rows = {drop_count}{colors.bcolors.END}',
              end="\r", file=sys.stderr)
    print()

    df_passed = df.drop(row_drop)

    kept_data = 100 * df_passed.sum().sum() / df.sum().sum()
    kept_rows = 100 * len(df_passed.index) / len(df.index)

    print(
        f'Output has {kept_rows:.2f}% of the original rows, but {kept_data:.2f}% of the original data.')

    print(f'Writing table to {out_file}')
    df_passed.to_csv(out_file, index=True, sep="\t")


if __name__ == "__main__":

    # OPTIONS #################################################################

    # print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

    parser = argparse.ArgumentParser(description="",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-f', '--file',
                        help="Table or matrix file",
                        required=True)

    parser.add_argument('-o', '--out',
                        help="Out file",
                        required=True)

    parser.add_argument('--min_abundance',
                        help="Minimum count per column",
                        required=False,
                        type=float,
                        default=3.0)

    parser.add_argument('--min_columns',
                        help="Minimum columns above abundance",
                        required=False,
                        type=int,
                        default=3)

    args = parser.parse_args()

    main(args.file, args.min_abundance, args.min_columns, args.out)
