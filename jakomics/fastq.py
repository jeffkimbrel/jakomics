import os
import uuid
import pandas as pd

from jakomics.file import FILE
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

'''
The fastq class is meant to primarily keep track of pairing information for
files. The files are stored in the self.files slot, and these are of the jakomics
FILE class - position 0 for F, I and S, and position 1 for R if paired end.

Methods generally run on the FILEs in the self.file slot.
'''


class FASTQ():

    def __str__(self):
        return "<JAKomics FASTQ class>"

    def __init__(self, sample, meta):
        self.sample = sample
        self.meta = meta
        self.infer_pair_type()

    def infer_pair_type(self):
        if pd.notnull(self.meta['I']):
            self.type = "Interleaved"
            self.files = [FILE(self.meta['I'])]
            self.files[0].read = 'FR'
        elif pd.notnull(self.meta['F']) and pd.notnull(self.meta['R']):
            self.type = "Paired"
            self.files = [FILE(self.meta['F']), FILE(self.meta['R'])]
            self.files[0].read = 'F'
            self.files[1].read = 'R'
        elif pd.notnull(self.meta['F']):
            self.type = "Single"
            self.files = [FILE(self.meta['F'])]
            self.files[0].read = 'F'
        else:
            self.type = "Unknown"
            self.files = []

    def verify_read_pairs(self):
        if self.type == "Paired":
            call = 'reformat.sh in1=' + self.files[0].file_path + \
                ' in2=' + self.files[1].file_path + ' verifypaired=t'
        elif self.type == "Interleaved":
            call = 'reformat.sh in=' + self.files[0].file_path + ' verifypaired=t'

        if self.type == "Paired" | self.type == "Interleaved":
            print(call)
            lines = system_call(call)

            if "Names appear to be correctly paired." in lines:
                self.ordered = True
                print("Reads are in the same order in both files.")
            else:
                print(
                    f'\n***  ERROR  ***\nThe read pairs are not in the same order in both files for {self.sample} \n*** EXITING ***\n')
                sys.exit()
        else:
            print(f'FastQ file for {self.sample} is not paired')


def run_info(file):
    headers = {}

    with gzip.open(file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            split = title.split(":")
            merge = split[0] + ":" + split[1] + ":" + split[2] + ":" + split[3]

            if merge in headers:
                headers[merge] += 1
            else:
                headers[merge] = 1

    return headers
