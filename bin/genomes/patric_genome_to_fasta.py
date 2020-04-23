import os
import sys
import argparse
import pandas as pd
import json
from multiprocessing import Manager, Pool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with top-level patric folders",
                    required=True)

parser.add_argument('--aa_out',
                    help="Directory to write protein sequences",
                    required=False,
                    default=None)

parser.add_argument('--nt_out',
                    help="Directory to write gene nt sequences",
                    required=False,
                    default=None)

args = parser.parse_args()

if args.aa_out is not None:
    if not os.path.exists(args.aa_out):
        print("Creating directory " + os.path.abspath(args.aa_out))
        os.makedirs(args.aa_out)

if args.nt_out is not None:
    if not os.path.exists(args.nt_out):
        print("Creating directory " + os.path.abspath(args.nt_out))
        os.makedirs(args.nt_out)

manager = Manager()
counter = manager.dict()

# CLASSES #####################################################################


class Patric():

    def __init__(self, full_name):
        self.full_name = full_name
        self.features_path = os.path.join(args.in_dir, full_name, full_name + ".txt")
        self.metadata_path = os.path.join(args.in_dir, full_name, "load_files", "genome.json")
        self.apply_metadata()

    def apply_metadata(self):
        metadata = json.loads(open(self.metadata_path, "r").read())
        self.genome_name = metadata[0]['genome_name']
        self.genome_id = metadata[0]['genome_id']

# FUNCTIONS ###################################################################


def get_files():
    genome_list = []
    folder_list = [d for d in os.listdir(
        args.in_dir) if os.path.isdir(os.path.join(args.in_dir, d))]

    for folder in folder_list:
        genome_list.append(Patric(folder))

    return(genome_list)


def make_seq_record(locus_tag, sequence, function):
    record = SeqRecord(
        Seq(
            str(sequence)
        ),
        id=locus_tag,
        description=locus_tag + " " + function
    )
    return record


def main(file):
    global counter

    print(f'> Finished {len(counter)} of {len(files)} genomes...',
          end="\r", file=sys.stderr)

    df = pd.read_csv(file.features_path, sep="\t")
    aa_records = []
    nt_records = []

    for id, feature in df.iterrows():
        locus_tag = feature['feature_id'].replace("fig|" + file.genome_id, file.genome_name)
        aa = feature['aa_sequence']
        nt = feature['nucleotide_sequence']

        # only make seq objects if sequence is present
        if pd.isna(aa) is False:
            aa_records.append(make_seq_record(locus_tag, aa, feature['function']))
        if pd.isna(nt) is False:
            nt_records.append(make_seq_record(locus_tag, nt, feature['function']))

    if args.aa_out is not None:
        SeqIO.write(aa_records, os.path.join(args.aa_out, file.genome_name + ".faa"), "fasta")
    if args.nt_out is not None:
        SeqIO.write(nt_records, os.path.join(args.nt_out, file.genome_name + ".ffn"), "fasta")

    counter[file.genome_name] = file.genome_name

    print(f'> Finished {len(counter)} of {len(files)} genomes...',
          end="\r", file=sys.stderr)


if __name__ == "__main__":

    files = get_files()
    pool = Pool()
    pool.map(main, files)
    print()  # print an extra line to overcome sys.stderr stuff
    pool.close()
