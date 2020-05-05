import os
import sys
import argparse
import json
from multiprocessing import Manager, Pool
from jakomics import colors
from natsort import natsorted

from Bio import SeqIO


# RAST output by default has an incorrect header which leads to a BioPython warning. This will suppress all warnings.
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with top-level patric folders",
                    required=True)

args = parser.parse_args()

# manager = Manager()
# counter = manager.dict()

# CLASSES #####################################################################

class Patric():

    def __init__(self, full_name):
        self.full_name = full_name
        self.gb_path = os.path.join(args.in_dir, full_name, full_name + ".gb")
        self.metadata_path = os.path.join(args.in_dir, full_name, "load_files", "genome.json")
        self.apply_metadata()
        self.gb_out_path = os.path.join(args.in_dir, full_name, self.genome_name + ".gbk")

    def apply_metadata(self):
        metadata = json.loads(open(self.metadata_path, "r").read())
        self.genome_name = metadata[0]['genome_name']
        self.genome_id = metadata[0]['genome_id']

# FUNCTIONS ###################################################################

def get_files():
    genome_list = []
    folder_list = [d for d in os.listdir(
        args.in_dir) if os.path.isdir(os.path.join(args.in_dir, d))]

    for folder in natsorted(folder_list):
        genome_list.append(Patric(folder))

    return(genome_list)

def add_locus_tags(file):

    if os.path.exists(file.gb_out_path):
        print(f'{colors.bcolors.RED}Overwriting {file.gb_out_path}{colors.bcolors.END}')
        os.remove(file.gb_out_path)

    print(f'Adding locus tags to {colors.bcolors.GREEN}{file.genome_name}{colors.bcolors.END} on PID {os.getpid()}')


    for seq_record in SeqIO.parse(file.gb_path, "genbank"):
        seq_record.id = seq_record.description

        new_features = []
        for feature in seq_record.features:
            for db_xref in feature.qualifiers['db_xref']:
                if db_xref.startswith('RAST2'):
                    feature.qualifiers['locus_tag'] = db_xref.replace('RAST2:fig|' + file.genome_id, file.genome_name)
                    # print(db_xref, feature.qualifiers['locus_tag'])
            new_features.append(feature)
        seq_record.features = new_features


        output_handle = open(file.gb_out_path, "a")
        SeqIO.write(seq_record, output_handle, "genbank")
        output_handle.close()

# MAIN ########################################################################

if __name__ == "__main__":

    files = get_files()

    pool = Pool()
    pool.map(add_locus_tags, files)
    pool.close()
