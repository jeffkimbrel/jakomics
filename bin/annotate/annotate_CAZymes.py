import sys
import os
import argparse
import uuid
from multiprocessing import Manager, Pool

from jakomics import hmm, utilities, colors
from bin.text import count_and_merge_tables

# print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

# hmm.test()
# utilities.test()
# colors.test()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='''

            Find CAZYmes in an amino acid file

            QC_CODE refers to different scoring rules set forth by DBCAN6:
                1A = E-value <= 1e-18 AND HMM coverage >= 35%
                1B = E-value <= 1e-15 AND HMM coverage >= 35%
                2 = E-value <= 1e-5  AND Alignment Length >= 80
                3 = E-value <= 1e-3  AND HMM coverage >= 30%
                0 = Low-quality hit that did not pass any metrics''',

                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('-b', '--berlemont',
                    action="store_true",
                    help="Save berlemont GH summary file")

parser.add_argument('--qc',
                    help="QC codes to include",
                    required=False,
                    nargs='*',
                    default=['1A', '1B'])

parser.add_argument('-o', '--out',
                    help="Write merged file",
                    required=False,
                    default=None)

parser.add_argument('-s', '--substrate',
                    help="Write merged substrates file",
                    required=False,
                    default=None)

args = parser.parse_args()

manager = Manager()
counter = manager.dict()

# MISC ########################################################################

extrasPath = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(sys.argv[0])))), "db")

dbcan_db_path = extrasPath + '/dbCAN-HMMdb-V8.txt'

# CLASSES #####################################################################

genes = {}


class File:

    def __init__(self, file_path):
        self.file_path = file_path
        self.run_id = utilities.get_unique_ID()
        self.file_name = os.path.basename(file_path)
        self.results_file = self.file_name + '.dbcan8.txt'
        self.results_class_file = self.file_name + '.class_summary.txt'
        self.results_berlemont_file = self.file_name + '.berlemont_summary.txt'
        self.temp_log = self.run_id + '.log'
        self.temp_output = self.run_id + '.temp.txt'

    def remove_temp(self):
        os.system('rm ' + self.temp_log)
        os.system('rm ' + self.temp_output)

    def process_dbcan_results(self):
        '''
        read the raw HMMER file and produce HMM classes
        '''

        rawResult = [line.strip() for line in open(self.temp_output)]

        self.cazymes = {}

        for line in rawResult:
            if not line.startswith("#"):

                cazyme = hmm.CAZYME(line)
                cazyme.pass_cazy()
                if cazyme.gene in self.cazymes:
                    self.cazymes[cazyme.gene].append(cazyme)
                else:
                    self.cazymes[cazyme.gene] = [cazyme]

    def write_results(self):
        hmm_output = open(self.results_file, 'w')
        hmm_output.write(
            "LOCUS\tHMM\tEVAL\tSCORE\tC-EVAL\tI-EVAL\tSEQ_COORDS\tHMM_COORDS\tALIGN_LENGTH\tHMM_COVERAGE\tQC_CODE\tCLASS\tSUBSTRATE\n")

        for gene in sorted(self.cazymes.keys()):
            for hit in self.cazymes[gene]:
                if hit.pass_qc in args.qc:
                    hit.assign_cazy_class()
                    hit.assign_substrate()
                    hmm_output.write(hit.write())

        hmm_output.close()


# FUNCTIONS ###################################################################

def processHMMER(file):
    '''
    read the raw HMMER file and produce HMM classes
    '''

    rawResult = [line.strip() for line in open(file.temp_output)]

    formattedResult = {}

    for line in rawResult:
        if not line.startswith("#"):

            formattedResult[uuid.uuid4().hex] = hmm.CAZYME(line)

    return(formattedResult)


def main(file_path):

    global counter

    print(f'> Finished {len(counter)} of {len(file_list)} files...',
          end="\r", file=sys.stderr)

    file = File(file_path)

    hmm.run_hmmsearch(file.file_path, file.temp_log, file.temp_output,
                      dbcan_db_path)

    file.process_dbcan_results()
    file.write_results()

    # cleanup
    file.remove_temp()

    # summary stuff
    counter[file.results_file] = file.file_name
    print(f'> Finished {len(counter)} of {len(file_list)} files...',
          end="\r", file=sys.stderr)


## MAIN LOOP ###################################################################


if __name__ == "__main__":

    file_list = utilities.get_file_list(args.files, [''])
    file_list = utilities.get_directory_file_list(args.in_dir, [''], file_list)

    if len(file_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid .faa files were found!{colors.bcolors.END}")

    pool = Pool()
    pool.map(main, file_list)
    pool.close()

    print("\n\n---\n")

    for file in sorted(counter.keys()):
        print(
            f'{file} results written to {colors.bcolors.GREEN}{counter[file]}{colors.bcolors.END}')

    print("\nFinished searching for cazymes!")

    # write merged results
    if args.out is not None:
        print("\n---\n")
        count_and_merge_tables.main(counter, args.out, 'HMM')

    if args.substrate is not None:
        print("\n---\n")
        count_and_merge_tables.main(counter, args.substrate, 'SUBSTRATE')
