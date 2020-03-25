import pandas as pd
from natsort import natsorted
import argparse

# OPTIONS #####################################################################
parser = argparse.ArgumentParser(description='TPM transform columns in a table')

parser.add_argument('-l', '--lengths',
                    help="Excel file with 'MAG' and 'length' columns",
                    required=True)
parser.add_argument('-t', '--table', help="Table with raw counts", required=True)

parser.add_argument('-o', '--out',
                    help="Filename to write TPM table to",
                    required=True)

args = parser.parse_args()

# FUNCTIONS ###################################################################


def get_lengths():
    lengths = pd.read_excel(args.lengths)
    lengths = lengths.set_index('MAG')
    lengths = lengths['length']
    return lengths


def prepare_table():
    table = pd.read_csv(args.table, sep="\t")
    table = table.set_index('MAG')
    table = table.fillna(0)
    return table


def calc_rpk(df):
    rpk = df.div(df.length, axis=0)
    rpk = rpk.drop(columns=['length'])
    return rpk


def calc_tpm(rpk, per_million):
    tpm = rpk.div(per_million, axis=1)
    tpm.index = natsorted(tpm.index)
    return tpm

# MAIN ########################################################################


if __name__ == "__main__":
    lengths = get_lengths()
    table = prepare_table()
    df = pd.concat([table, lengths], axis=1, sort=True)
    rpk = calc_rpk(df)
    per_million = rpk.sum(axis=0) / 1000000
    tpm = calc_tpm(rpk, per_million)
    tpm.to_csv(args.out, sep="\t")
