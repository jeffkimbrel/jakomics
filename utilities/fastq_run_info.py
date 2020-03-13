import gzip
import argparse
import subprocess
import os
import sys
from natsort import natsorted

from Bio.SeqIO.QualityIO import FastqGeneralIterator

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Check the Illumina run information on a directory (--in_dir) or list (-f) of fastq.gz files')

parser.add_argument('-f', '--files',
                    help="fastq.gz files",
                    nargs='*',
                    required=False,
                    default=None)

parser.add_argument('--in_dir', help="Directory with fastq.gz files", required=False, default=None)

args = parser.parse_args()

# FUNCTIONS ###################################################################


def get_file_list():
    # combine args.files and args.in_dir to a master file list

    file_list = []

    if args.in_dir is not None:
        for file in os.listdir(args.in_dir):
            if file.endswith('fastq.gz'):
                abs_path = os.path.abspath(file)
                file_list.append(abs_path)

    if args.files is not None:
        for file in args.files:
            if file.endswith('fastq.gz'):
                abs_path = os.path.abspath(file)
                file_list.append(abs_path)

    file_list = list(set(file_list))

    if len(file_list) == 0:
        sys.exit("Error: No valid fastq.gz files given")

    return file_list


def get_md5(file):
    call = "md5 " + file
    p1 = subprocess.Popen(call, shell=True,
                          stdin=None,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    out, err = p1.communicate()
    out = out.decode()
    junk, md5 = out.split(" = ")
    return(md5.strip())


def count_headers(file, md5):
    headers = {}

    with gzip.open(file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            split = title.split(":")
            merge = split[0] + ":" + split[1] + ":" + split[2] + ":" + split[3] + ":" + md5

            if merge in headers:
                headers[merge] += 1
            else:
                headers[merge] = 1

    return(headers)


def format_out(info):
    print("FILE", "MD5", "INSTRUMENT", "RUN", "FLOWCELL", "LANE", "READS", sep=" | ")
    print(":---", ":---", ":---", ":---", ":---", ":---", ":---", sep=" | ")

    for file in info:
        for header in info[file]:
            instrument, run, flowcell, lane, md5 = header.split(":")
            print(file, md5, instrument, run, flowcell, lane, info[file][header], sep="|")


# MAIN ########################################################################


info = {}

file_list = get_file_list()

for file in natsorted(file_list):

    md5 = get_md5(file)
    headers = count_headers(file, md5)

    info[file] = headers

format_out(info)
