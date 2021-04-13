import os
import sys
import uuid
import pandas as pd

from jakomics.file import FILE
from jakomics.utilities import system_call
from jakomics import colors, bbtools

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

        self.processed_fastq = []
        for file in self.files:
            self.processed_fastq.append(os.path.join(file.dir, file.name))

        self.processed_sample_name = self.sample

    def view(self):
        print(self.sample, self.files)

    def verify_read_pairs(self, echo=False, run=True):
        if self.type == "Paired":
            call = 'reformat.sh in1=' + self.processed_fastq[0] + \
                ' in2=' + self.processed_fastq[1] + ' verifypaired=t'
            lines = system_call(call, echo=echo, run=run)
        elif self.type == "Interleaved":
            call = 'reformat.sh in=' + self.processed_fastq[0] + ' verifypaired=t'
            lines = system_call(call, echo=echo, run=run)
        else:
            print(f'{colors.bcolors.YELLOW}{self.sample} file is a single direction only{colors.bcolors.END}')
            call = None

        if call is not None:

            if "Names appear to be correctly paired." in lines:
                self.ordered = True
                # print(
                #     f"{colors.bcolors.GREEN}{self.sample} reads appear to be correctly paired{colors.bcolors.END}")
            else:
                self.ordered = False
                print(
                    f"{colors.bcolors.RED}{self.sample} reads are not correctly paired... exiting{colors.bcolors.END}")
                sys.exit()

    def adapter_trimming(self, db, echo=False, run=True, mem="Xmx8g", threads=8):

        self.processed_sample_name = self.processed_sample_name + ".rt"
        self.rt = []
        if self.type == "Paired":
            in1 = self.processed_fastq[0]
            self.rt.append(os.path.dirname(self.processed_fastq[0]) + "/" +
                           self.processed_sample_name + ".R1.fastq.gz")

            in2 = self.processed_fastq[1]
            self.rt.append(os.path.dirname(self.processed_fastq[1]) + "/" +
                           self.processed_sample_name + ".R2.fastq.gz")

            call = 'bbduk.sh in1=' + in1 + ' in2=' + in2 + ' out1=' + self.rt[0] + \
                ' out2=' + self.rt[1] + ' stats=' + self.processed_sample_name + '_stats.txt ref=' + db + \
                ' t=' + str(threads) + \
                ' ftl=5 ktrim=r k=23 mink=11 hdist=1 tpe tbo minlen=50 -' + mem

        elif self.type == "Interleaved":
            in1 = self.processed_fastq[0]
            self.rt.append(os.path.dirname(self.processed_fastq[0]) + "/" +
                           self.processed_sample_name + ".fastq.gz")

            call = 'bbduk.sh in=' + in1 + ' out=' + self.rt[0] + \
                ' stats=' + self.processed_sample_name + '_stats.txt ref=' + db + \
                ' t=' + str(threads) + \
                ' ftl=5 ktrim=r k=23 mink=11 hdist=1 tpe tbo minlen=50 -' + mem

        bb = bbtools.extract_stats(system_call(call, echo=echo, run=run))
        self.processed_fastq = self.rt
        return(bb)

    def contaminant_filtering(self, db, echo=False, run=True, mem="Xmx8g", threads=8):

        self.processed_sample_name = self.processed_sample_name + ".cf"
        self.cf = []
        if self.type == "Paired":
            in1 = os.path.join(self.files[0].dir, os.path.basename(self.processed_fastq[0]))
            self.cf.append(os.path.join(self.files[0].dir,
                                        self.processed_sample_name + ".R1.fastq.gz"))

            in2 = os.path.join(self.files[1].dir, os.path.basename(self.processed_fastq[1]))
            self.cf.append(os.path.join(self.files[1].dir,
                                        self.processed_sample_name + ".R2.fastq.gz"))

            for file in self.files:
                print(file.__dict__)

            call = 'bbduk.sh in1=' + in1 + ' in2=' + in2 + ' out1=' + \
                self.cf[0] + ' out2=' + self.cf[1] + ' stats=' + \
                self.processed_sample_name + '_stats.txt t=' + str(threads) + \
                ' ref=' + db + ' k=31 hdist=1 minlen=50 -' + mem

        elif self.type == "Interleaved":
            in1 = os.path.join(self.files[0].dir, os.path.basename(self.processed_fastq[0]))
            self.cf.append(os.path.join(self.files[0].dir,
                                        self.processed_sample_name + ".fastq.gz"))

            call = 'bbduk.sh in=' + in1 + ' out=' + \
                self.cf[0] + ' stats=' + \
                self.processed_sample_name + '_stats.txt t=' + str(threads) + \
                ' ref=' + db + ' k=31 hdist=1 minlen=50 -' + mem

        bb = bbtools.extract_stats(system_call(call, echo=echo, run=run))
        self.processed_fastq = self.cf
        return(bb)

    def quality_filtering(self, echo=False, run=True, mem="Xmx8g", threads=8):
        self.processed_sample_name = self.processed_sample_name + ".qf"
        self.qf = []
        if self.type == "Paired":
            in1 = self.processed_fastq[0]
            self.qf.append(os.path.dirname(self.processed_fastq[0]) + "/" +
                           self.processed_sample_name + ".R1.fastq.gz")

            in2 = self.processed_fastq[1]
            self.qf.append(os.path.dirname(self.processed_fastq[1]) + "/" +
                           self.processed_sample_name + ".R2.fastq.gz")

            call = 'bbduk.sh in1=' + in1 + ' in2=' + in2 + ' out1=' + \
                self.qf[0] + ' out2=' + self.qf[1] + \
                ' stats=' + \
                self.processed_sample_name + '_stats.txt t=' + \
                str(threads) + ' qtrim=r trimq=10 minlen=50 -' + mem

        elif self.type == "Interleaved":
            in1 = self.processed_fastq[0]
            self.qf.append(os.path.dirname(self.processed_fastq[0]) + "/" +
                           self.processed_sample_name + ".fastq.gz")

            call = 'bbduk.sh in=' + in1 + ' out=' + \
                self.qf[0] + ' stats=' + \
                self.processed_sample_name + '_stats.txt t=' + \
                str(threads) + ' qtrim=r trimq=10 minlen=50 -' + mem

        bb = bbtools.extract_stats(system_call(call, echo=echo, run=run))
        self.processed_fastq = self.qf
        return(bb)


def run_info(file):
    headers = {}

    if file.endswith(".gz"):
        with gzip.open(file, "rt") as handle:
            for title, seq, qual in FastqGeneralIterator(handle):
                split = title.split(":")
                merge = split[0] + ":" + split[1] + ":" + split[2] + ":" + split[3]

                if merge in headers:
                    headers[merge] += 1
                else:
                    headers[merge] = 1
    else:
        with open(file, "rU") as handle:
            for title, seq, qual in FastqGeneralIterator(handle):
                split = title.split(":")
                merge = split[0] + ":" + split[1] + ":" + split[2] + ":" + split[3]

                if merge in headers:
                    headers[merge] += 1
                else:
                    headers[merge] = 1

    return headers
