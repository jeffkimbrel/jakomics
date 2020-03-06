import pandas as pd
import argparse
import os
from Bio import SeqIO

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-s', '--stats', help="Metawrap Stats File", required=True)
parser.add_argument('-n', '--name', help="Name for MAGs", required=True)

parser.add_argument('--in_dir', help="Directory with original name bins", required=True)
parser.add_argument(
    '--out_dir', help="Directory to write renamed bins and stats file", required=True)

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'
args.out_dir = os.path.abspath(args.out_dir) + '/'

print("\nMetaWrap Stats File:\n\t", args.stats)
print("Original Bins Directory:\n\t", args.in_dir)
print("Renamed Bins Directory:\n\t", args.out_dir)

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Read and Sort Stats File ####################################################
stats = pd.read_csv(args.stats, sep="\t")
stats = stats.sort_values(['completeness', 'contamination'], ascending=[False, True])
stats = stats.reset_index(drop=True)
stats["MAG"] = stats.index + 1
stats["MAG"] = args.name + stats["MAG"].astype(str)
stats['OLD_PATH'] = args.in_dir + stats['bin'] + '.fa'
stats['NEW_PATH'] = args.out_dir + stats['MAG'] + '.fa'
stats = stats.set_index('MAG')
stats["contigs"] = 0

# Interate over bins ##########################################################

for MAG, row in stats.iterrows():
    counter = 1
    renamed = []
    for record in SeqIO.parse(row['OLD_PATH'], "fasta"):
        record.id = MAG + '_' + str(counter)
        renamed.append(record)
        counter += 1
    stats.loc[MAG, 'contigs'] = len(renamed)
    SeqIO.write(renamed, row['NEW_PATH'], "fasta")

# Print Stats File ############################################################
stats = stats[['bin', 'completeness', 'contamination', 'GC',
               'lineage', 'contigs', 'N50', 'size', 'OLD_PATH', 'NEW_PATH']]

stats = stats.rename(columns={'lineage': 'checkM_lineage'})

out_file = args.out_dir + args.name + "_rename_stats.txt"

stats.to_csv(out_file, sep="\t")

print("\nNew stats file:\n\t", out_file)
print("\nFinished processing " + str(len(stats.index)) + " MAGs!!")
