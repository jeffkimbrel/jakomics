import uuid
import os
import subprocess
import re

from jakomics import colors


class KOFAM:

    def __init__(self, line):

        parsed = re.split('\t', line)

        self.gene = parsed[1]
        self.KO = parsed[2]
        self.threshold = parsed[3]
        self.score = float(parsed[4])
        self.evalue = float(parsed[5])
        self.description = parsed[6]

        if len(self.threshold) == 0:
            self.threshold = 0
            print(f'{colors.bcolors.YELLOW}WARNING: {self.KO} does not have a KO threshold. All hits >0 will be included.{colors.bcolors.END}')

        if self.score >= float(self.threshold):
            self.threshold = float(self.threshold)
            self.passed = True
        else:
            self.passed = False

    def view(self):
        return [self.gene, self.KO, self.threshold, self.score, self.evalue]

    def result(self):
        return {'gene': self.gene,
                'annotation': self.KO,
                'score': self.score,
                'evalue': self.evalue}

    def __str__(self):
        return "<JAKomics KOFAM class>"


def run_kofam(faa_path, hal_path, cpus=1):
    temp_dir = 'KO_' + uuid.uuid4().hex

    command = 'exec_annotation --no-report-unannotated --tmp-dir ' + \
        temp_dir + ' ' + faa_path + ' --cpu ' + str(int(cpus))
    command = command + ' --profile ' + hal_path + ' -f detail-tsv ; rm -fR ' + temp_dir

    # print(command)
    kofam_results = subprocess.Popen(command, shell=True,
                                     stdin=None,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

    out, err = kofam_results.communicate()
    hits = []
    for line in out.decode().split("\n"):
        if len(line) > 0 and not line.lstrip().startswith('#'):
            hits.append(KOFAM(line))

    return hits


def parse_kofam_hits(run_kofam_out):
    '''
    Returns a dictionary of passed results with KO as key and list of kofam classes as value
    '''
    parsed = {}
    for hit in run_kofam_out:
        if hit.passed:
            # print(db['DB_NAME'], genome.short_name, hit.view(), sep="\t")
            if hit.KO in parsed:
                parsed[hit.KO].append(hit)
            else:
                parsed[hit.KO] = [hit]
    return parsed
