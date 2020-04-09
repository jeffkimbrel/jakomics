import argparse
import os
import subprocess

from utilities import bbtools

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('--in_dir',
                    help="Directory with unfiltered fastq.gz files",
                    required=True)
parser.add_argument('--out_dir',
                    help="Directory to write filtered fastq.gz files",
                    default="contaminant_filtered",
                    required=False)
parser.add_argument('-c', '--contaminants',
                    help="Fasta file with contaminant sequences",
                    default=os.path.dirname(os.path.abspath(sys.argv[0])) + "/db/contam_seqs.fa",
                    required=False)

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'
args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Classes #####################################################################


class SAMPLE:

    def __init__(self, F_file):
        self.sample = F_file.replace("_R1.fastq.gz", "")
        self.F_file = F_file
        self.R_file = self.sample + "_R2.fastq.gz"
        self.F_in_path = args.in_dir + self.F_file
        self.R_in_path = args.in_dir + self.R_file
        self.F_out_path = args.out_dir + self.F_file
        self.R_out_path = args.out_dir + self.R_file
        self.stats_file = args.out_dir + self.sample + "_stats.txt"

    def p(self):
        print(f'---\nSAMPLE: {self.sample}\nF_out : {self.F_out_path}\nR_out: {self.R_out_path}')

    def filter(self, c):
        command = 'bbduk.sh -Xmx8g in1=' + self.F_in_path + ' in2=' + self.R_in_path + ' out1=' + self.F_out_path + \
            ' out2=' + self.R_out_path + ' ref=' + c + ' k=31 hdist=1 stats=' + self.stats_file

        # print("\n>$ "+command)
        subprocess.run('source activate ~/bin/bbmap/ && ' +
                       command + ' && conda deactivate', shell=True, capture_output=True)

# Main ########################################################################


samples = []

for file in os.listdir(args.in_dir):
    if file.endswith('R1.fastq.gz'):
        samples.append(SAMPLE(file))


for sample in samples:
    sample.filter(args.contaminants)
    total, matched = bbtools.bbduk_stats_parser(sample.stats_file)

    print(f'{matched} of {total} reads were considered contaminants for {sample.sample}')
