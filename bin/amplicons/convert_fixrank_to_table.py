import pandas as pd
import argparse

from jakomics.taxa import RDP

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-f', '--file',
                    help="RDP fixrank file",
                    required=True)

parser.add_argument('-t', '--threshold',
                    help="RDP confidence threshold",
                    type=float,
                    required=True)

args = parser.parse_args()

if args.threshold > 1:
    args.threshold /= 100

# MAIN ########################################################################


if __name__ == "__main__":

    column_names = ['ASV', 'STRAND', 'D', 'D_p', 'P', 'P_p',
                    'C', 'C_p', 'O', 'O_p', 'F', 'F_p', 'G', 'G_p']

    df = pd.read_csv(args.file,
                     skiprows=7,
                     sep=";",
                     header=None,
                     names=column_names,
                     index_col=None)

    print("", "domain", "phylum", "class", "order", "family", "genus", sep="\t")

    for id, row in df.iterrows():
        asv = RDP(row, args.threshold)
        asv.view()
