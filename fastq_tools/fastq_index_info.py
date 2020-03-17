import argparse
import pandas as pd
from Bio import SeqIO
import sys

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='X')
parser.add_argument('-i1', '--index1', help="Index 1 fastq", required=True)
parser.add_argument('-i2', '--index2', help="Index 2 fastq", required=True)
parser.add_argument('-m',  '--mapping', help="Mapping File", required=True)

args = parser.parse_args()

index2_junk = 'TCTTTCCC'

# get mapping file

map = pd.read_csv(args.mapping, sep="\t")
map = map.set_index("BarcodeSequence")
map = map.rename(columns={'#SampleID': 'Sample'})

# Save Barcode Counts
barcodes = {}
counter = 0
index2_counter = 0

for index1, index2 in zip(SeqIO.parse(args.index1,  "fastq"),
                          SeqIO.parse(args.index2,  "fastq")):
    counter += 1

    if counter % 100000 == 0:
        print(f'Processed {counter} read pairs', end="\r", file=sys.stderr)

    F = str(index1.seq)
    R = str(index2.seq)

    if R == index2_junk:
        index2_counter += 1
    else:
        if F in barcodes:
            if R in barcodes[F]:
                barcodes[F][R] += 1
            else:
                barcodes[F][R] = 1
        else:
            barcodes[F] = {R: 1}

#

unmapped_F = {}
unmapped_R = {}
unmapped_total = 0

# final print
print(print(f'Processed {counter} read pairs', file=sys.stderr))

for barcode in sorted(barcodes.keys()):
    for reverse in barcodes[barcode]:
        merged = barcode + reverse
        sample = "-"

        if merged in map.index:
            sample = map.loc[merged, "Sample"]
        else:
            unmapped_F[barcode] = unmapped_F.get(barcode, 0) + barcodes[barcode][reverse]
            unmapped_R[reverse] = unmapped_R.get(reverse, 0) + barcodes[barcode][reverse]
            unmapped_total += barcodes[barcode][reverse]

        if barcodes[barcode][reverse] > 1000 or sample != '-':
            print(barcode + reverse, barcodes[barcode][reverse], sample, sep="\t")

print('Index2_TCTTTCCC', index2_counter, "-", sep="\t")
print('Total_Unmapped', unmapped_total, sep="\t")

print("-----")
for w in sorted(unmapped_F, key=unmapped_F.get, reverse=True):
    if unmapped_F[w] > 1000:
        print("FORWARD", w, unmapped_F[w], sep="\t")

print("-----")
for w in sorted(unmapped_R, key=unmapped_R.get, reverse=True):
    if unmapped_R[w] > 1000:
        print("REVERSE", w, unmapped_R[w], sep="\t")
