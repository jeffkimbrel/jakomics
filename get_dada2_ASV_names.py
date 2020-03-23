import argparse
import pandas as pd
import sys
from natsort import natsorted

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-i', '--id', help="ASV ID", required=True)

parser.add_argument('-f1', '--fasta1', help="Original Fasta File", required=True)
parser.add_argument('-t1', '--table1', help="ASV table with seq names", required=True)

parser.add_argument('-f2', '--fasta2', help="New Fasta File", required=True)
parser.add_argument('-t2', '--table2', help="New ASV table with merged names", required=True)

args = parser.parse_args()

# Classes #####################################################################


class ORIGINAL:
    def __init__(self, record):
        self.record = record
        self.id = int(record.id.replace(args.id, ""))
        self.seq = str(record.seq)
        self.found_in_table = False

# Functions ###################################################################


def get_last_id(originals):
    last = 0
    for original in originals:
        if original.id > last:
            last = original.id
    return(last)


def missing_fasta_seqs(originals):
    missing = []
    for original in originals:
        if original.found_in_table == False:
            missing.append(original.record)

    if len(missing) > 0:
        print("Uh-oh, there FASTA seqs not found in the ASV table. Please re-run on the full FASTA dataset.")
        for record in missing:
            print(record.id, record.seq, sep="\t")
        sys.exit()


def write_fasta(t):

    # natural sort the df on ID
    t['ID'] = pd.Categorical(t['ID'], ordered=True, categories=natsorted(t['ID'].unique()))
    t = t.sort_values('ID')

    seqrecords = []
    for SEQ, row in t.iterrows():
        record = SeqRecord(
            Seq(
                SEQ
            ),
            id=row['ID'],
            description=''
        )
        seqrecords.append(record)

    SeqIO.write(seqrecords, args.fasta2, "fasta")


# Main ########################################################################


# # get original info from fasta into individual classes
originals = []

for record in SeqIO.parse(args.fasta1, "fasta"):
    originals.append(ORIGINAL(record))

# # get ASV table

asv_table = pd.read_csv(args.table1, sep="\t")
asv_table['ID'] = None
asv_table = asv_table.set_index('SEQ')

# # add IDs to existing ASVs in the table

for original in originals:
    if original.seq in asv_table.index:
        asv_table.loc[original.seq, 'ID'] = original.id
        original.found_in_table = True

# # Find any fasta seqs that are missing from the ASV table
missing_fasta_seqs(originals)
last_id = get_last_id(originals) + 1

print(f'Will start numbering new ASVs at {args.id}{last_id}')
# print(asv_table)

# number new ASVs
for SEQ, row in asv_table.iterrows():
    if row['ID'] == None:
        asv_table.loc[SEQ, 'ID'] = args.id + str(last_id)
        last_id += 1
    else:
        asv_table.loc[SEQ, 'ID'] = args.id + str(row['ID'])

write_fasta(asv_table)

asv_table = asv_table.set_index('ID')
asv_table = asv_table.reindex(index=natsorted(asv_table.index))
asv_table.to_csv(args.table2, sep="\t")
