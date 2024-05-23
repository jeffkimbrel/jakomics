import os
import argparse
import pandas as pd
from natsort import natsorted

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Subset a table on whether a column value is found in a column from a subset table',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-t', '--table', help="Table with full data", required=True)
parser.add_argument('-s', '--subset', help="Table with subset data", required=True)
parser.add_argument('-o', '--out', help="Path to write subsetted table", required=True)

parser.add_argument('--tc', help="Table Column with IDs (0-based)",
                    required=False, default=0, type=int)
parser.add_argument('--sc', help="Subset Column with IDs (0-based)",
                    required=False, default=0, type=int)
parser.add_argument('--strip', help="Strip this text off either ID", required=False, default="")

args = parser.parse_args()

# Get IDs from Subset #########################################################

subset_array = []

s = pd.read_csv(args.subset, sep="\t", engine='python')

for j in natsorted(s[s.columns[args.sc]]):
    j = j.replace(args.strip, "")
    subset_array.append(j)

# Now Subset the Main Table ###################################################

t = pd.read_csv(args.table, sep=None, engine='python')

row_drop = []

for index, row in t.iterrows():
    if row[args.tc].replace(args.strip, "") not in subset_array:
        row_drop.append(index)

subsetted_t = t.drop(row_drop)

subsetted_t.to_csv(args.out, index=False, sep="\t")
