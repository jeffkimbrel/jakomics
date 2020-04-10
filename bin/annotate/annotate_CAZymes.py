import sys
import os
import argparse
import uuid
import re
import operator  # for sorting by class attribute

from jakomics import hmm

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='''

Find CAZYmes in an amino acid file

QC_CODE refers to different scoring rules set forth by DBCAN6:
1A = E-value <= 1e-18 AND HMM coverage >= 35%
1B = E-value <= 1e-15 AND HMM coverage >= 35%
2 = E-value <= 1e-5  AND Alignment Length >= 80
3 = E-value <= 1e-3  AND HMM coverage >= 30%
0 = Low-quality hit that did not pass any metrics

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-a', '--amino',
                    help="Amino acid file (fasta)",
                    required=True)

parser.add_argument('-n', '--name',
                    help="Name of files")

parser.add_argument('-b', '--berlemont',
                    action="store_true",
                    help="Save berlemont GH summary file")

parser.add_argument('-c', '--cazyclass',
                    action="store_true",
                    help="Save CAZy classsummary file")

args = parser.parse_args()

if args.name == None:
    args.name = os.path.basename(args.amino)

# MISC ########################################################################

extrasPath = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(sys.argv[0])))), "db")

hmm_output = open(args.name + ".dbcan8.txt", 'w')

berlemont_file = [line.strip() for line in open(extrasPath + '/berlemont_types.txt')]
berlemont_lookup = {}
for line in berlemont_file:
    berlemont_lookup[line.split("\t")[0]] = line.split("\t")[1]

# CLASSES #####################################################################

genes = {}


class GENE:

    def __init__(self, id):
        self.id = id
        self.cazymes = []

    def add_cazyme(self, cazyme):
        self.cazymes.append(cazyme)

    def class_count(self, class_summary):
        unique = list(set(self.cazymes))
        for cazyme in unique:
            cazy_class = re.split(r'(\d+)', cazyme)[0]
            class_summary[cazy_class] = class_summary.get(cazy_class, 0) + 1

        return(class_summary)

    def berlemont_count(self, berlemont_summary):
        stripped = []
        for cazyme in self.cazymes:
            stripped.append(cazyme.split("_")[0])

        unique = list(set(stripped))
        for cazyme in unique:
            if cazyme in berlemont_lookup:
                berlemont = berlemont_lookup[cazyme]
                berlemont_summary[berlemont] = berlemont_summary.get(berlemont, 0) + 1

        return(berlemont_summary)

# FUNCTIONS ###################################################################


def systemCall(command):
    # print(">$ " + command)
    os.system(command)


def getID():
    return(uuid.uuid4().hex)


def runHMMER(file):
    command = 'hmmsearch -o ' + run_id + '.log --domT 10 --domtblout ' + run_id + '.raw.txt ' + \
        extrasPath + '/dbCAN-HMMdb-V8.txt ' + file
    print(command)
    systemCall(command)


def removeTemp():
    systemCall('rm ' + run_id + '.log')
    systemCall('rm ' + run_id + '.raw.txt')


def processHMMER():
    '''
    read the raw HMMER file and produce HMM classes
    '''

    rawResult = [line.strip() for line in open(run_id + '.raw.txt')]

    formattedResult = {}

    for line in rawResult:
        if not line.startswith("#"):

            formattedResult[uuid.uuid4().hex] = hmm.HMM(line)

    return(formattedResult)

## MAIN LOOP ###################################################################


run_id = getID()

runHMMER(args.amino)

HMMs = processHMMER()

hmm_output.write(
    "LOCUS\tHMM\tEVAL\tSCORE\tC-EVAL\tI-EVAL\tSEQ_COORDS\tHMM_COORDS\tALIGN_LENGTH\tHMM_COVERAGE\tDESCRIPTION\tQC_CODE\n")

# for result in HMMs:
for result in (sorted(HMMs.values(), key=operator.attrgetter('gene'))):
    result.passCAZY()
    if result.passQC in ['1A', '1B', '2', '3']:
        hmm_output.write(result.write())
        if result.gene not in genes:
            genes[result.gene] = GENE(result.gene)

        genes[result.gene].add_cazyme(result.model)
hmm_output.close()

print("Wrote HMM data to " + args.name + ".dbcan8.txt")

removeTemp()

# EXTRAS ######################################################################

class_summary = {}
berlemont_summary = {}

for gene in genes:
    class_summary = genes[gene].class_count(class_summary)
    berlemont_summary = genes[gene].berlemont_count(berlemont_summary)

if args.cazyclass is True:
    cazy_output = open(args.name + ".class_summary.txt", 'w')
    cazy_output.write("CLASS\tCOUNT\n")
    for cazy_class in class_summary:
        cazy_output.write(cazy_class + "\t" + str(class_summary[cazy_class]) + "\n")
    cazy_output.close()
    print("Wrote CAZy class summary to " + args.name + ".class_summary.txt")


if args.berlemont is True:
    berlemont_output = open(args.name + ".berlemont_summary.txt", 'w')
    berlemont_output.write("SUBSTRATE\tCOUNT\n")
    for berlemont in berlemont_summary:
        berlemont_output.write(berlemont + "\t" + str(berlemont_summary[berlemont]) + "\n")
    berlemont_output.close()
    print("Wrote Berlemont substrate summary to " + args.name + ".berlemont_summary.txt")
