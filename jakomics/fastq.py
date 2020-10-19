import os
import uuid
import pandas as pd

from jakomics.file import FILE
from jakomics.utilities import system_call
from jakomics import colors

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

        self.processed_fastq = self.files
        self.processed_sample_name = self.sample

    def view(self):
        print(self.sample, self.files)

    def verify_read_pairs(self):
        if self.type == "Paired":
            call = 'reformat.sh in1=' + self.processed_fastq[0].file_path + \
                ' in2=' + self.processed_fastq[1].file_path + ' verifypaired=t'
            lines = system_call(call)
        elif self.type == "Interleaved":
            call = 'reformat.sh in=' + self.processed_fastq[0].file_path + ' verifypaired=t'
            lines = system_call(call)
        else:
            print(f'{colors.bcolors.YELLOW}{self.sample} file is a single direction only{colors.bcolors.END}')
            call = None

        if call is not None:

            if "Names appear to be correctly paired." in lines:
                self.ordered = True
                print(
                    f"{colors.bcolors.GREEN}{self.sample} reads appear to be correctly paired{colors.bcolors.END}")
            else:
                self.ordered = False
                print(
                    f"{colors.bcolors.RED}{self.sample} reads are not correctly paired... exiting{colors.bcolors.END}")
                sys.exit()

    def adapter_trimming(self, db, echo=False, run=True):
        self.processed_sample_name = self.processed_sample_name + ".rt"
        if self.type == "Paired":
            in1 = self.processed_fastq[0].file_path
            out1 = self.processed_fastq[0].dir + "/" + self.processed_sample_name + ".R1.fastq.gz"

            in2 = self.processed_fastq[1].file_path
            out2 = self.processed_fastq[1].dir + "/" + self.processed_sample_name + ".R2.fastq.gz"

            call = 'bbduk.sh in1=' + in1 + ' in2=' + in2 + ' out1=' + out1 + \
                ' out2=' + out2 + ' stats=' + self.processed_sample_name + '_stats.rt.txt' + ' ref=' + db + \
                ' t=8 ftl=5 ktrim=r k=23 mink=11 hdist=1 tpe tbo minlen=50 -Xmx8g'

            system_call(call, echo=echo, run=run)


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
