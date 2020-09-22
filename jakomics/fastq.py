import os
import uuid
import pandas as pd


class FASTQ():

    def __str__(self):
        return "<JAKomics FASTQ class>"

    def __init__(self, sample, row):
        self.sample = sample
        self.F = row['F']
        self.R = row['R']
        self.I = row['I']

        if pd.notnull(self.I):
            self.type = "Interleaved"
        elif pd.notnull(self.F) and pd.notnull(self.R):
            self.type = "Paired"
        elif pd.notnull(self.F):
            self.type = "Single"
        else:
            self.type = "Unknown"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='test')

    parser.add_argument('-s', '--samples',
                        help="excel file with samples in S, F, R, I columns",
                        required=True)

    args = parser.parse_args()

    files = pd.read_excel(args.samples, index_col=0)
    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        print(d.sample, d.type, sep="\t")
