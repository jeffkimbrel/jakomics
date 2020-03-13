import pandas as pd
import argparse
import os
from os import path
import sys

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Submit a folder of MAGs to Patric. Must log-in to the Patric shell first.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--in_dir', help="Directory with original name bins", required=True)
parser.add_argument('-t', '--taxa', help="File with taxa information", required=True)
parser.add_argument('-p', '--patric', help="Patric workspace path", required=True)

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'


# Classes #####################################################################

class MAG:

    def __init__(self, mag_name, tax_name, tax_id):
        self.mag_name = mag_name
        self.tax_name = tax_name.replace(" ", "_")
        self.tax_id = str(tax_id)
        self.path = args.in_dir + mag_name + ".fa"

    def validate_path(self):
        try:
            f = open(self.path)
            f.close()
        except IOError:
            sys.exit("File not found: " + self.path)

    def patric(self):
        command = 'p3-submit-genome-annotation --contigs-file ' + self.path + ' -n "' + self.mag_name + '" -t ' + self.tax_id + \
            ' -d Bacteria /jkimbrel@patricbrc.org/' + args.patric + \
            ' "' + self.mag_name + '_' + self.tax_name + '"'
        print(command)
        os.system(command)


# Main ########################################################################

taxids = pd.read_csv(args.taxa, sep=None, engine='python')
mags = {}


for index, row in taxids.iterrows():
    mags[index] = MAG(row['GENOME'], row['NAME'], row['NCBI_ID'])
    mags[index].validate_path()
    mags[index].patric()
