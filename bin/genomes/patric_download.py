from natsort import natsorted
import argparse
import os
from os import path
import sys
import subprocess

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Download genomes in a Patrick workspace. Must log-in to the Patric shell first.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-p', '--patric', help="Patric workspace path", required=True)
parser.add_argument('--out_dir', help="Directory to download the Patric data to", required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Functions ###################################################################


def get_genome_list():
    genome_list = []

    command = 'p3-ls -l --type /jkimbrel@patricbrc.org/' + args.patric
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    ls_results = stdout.decode().split("\n")

    for result in ls_results:
        split = result.split()
        if len(split) > 0:
            if split[6] == 'job_result':
                genome_list.append('\ '.join(split[7:]))
    return genome_list


def download_patric_folder(genome):

    command_main = 'p3-cp -R ws:/jkimbrel@patricbrc.org/' + \
        args.patric + '/.' + genome + ' ' + args.out_dir + genome
    print(command_main)
    os.system(command_main)

# Main ########################################################################


genome_list = get_genome_list()
for genome in natsorted(genome_list):
    download_patric_folder(genome)
