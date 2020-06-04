import uuid
import os
import subprocess
import re

class KOFAM:

    def __init__(self, line):

        parsed = re.split(' +', line)

        self.gene = parsed[1]
        self.KO = parsed[2]
        self.treshold = float(parsed[3])
        self.score = float(parsed[4])
        self.evalue = float(parsed[5])
        self.description = ' '.join(parsed[6:])

        if self.score >= self.treshold:
            self.passed = True
        else:
            self.passed = False

    def view(self):
        return [self.gene, self.KO, self.treshold, self.score, self.evalue]


def run_kofam(faa_path, hal_path):
    temp_dir = 'KO_' + uuid.uuid4().hex

    command = '/Users/kimbrel1/Science/kofam_93/kofamscan-1.1.0/exec_annotation --no-report-unannotated --tmp-dir ' + temp_dir + ' ' + faa_path + ' --cpu 1'
    command = command + ' --profile ' + hal_path + '; rm -fR ' + temp_dir

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
