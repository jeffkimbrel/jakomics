import uuid
import os
import subprocess
import re
import pandas as pd

from jakomics import colors
from jakomics.utilities import system_call


class KOFAM:

    def __init__(self, line, t_scale=1.0, score_as_ratio = False):

        parsed = re.split('\t', line)
        self.parsed = parsed
        self.gene = parsed[1]
        self.KO = parsed[2]
        self.threshold = float(parsed[3])
        self.evalue = float(parsed[5])
        self.description = parsed[6]
        self.score_as_ratio = score_as_ratio

        if self.score_as_ratio:
            if len(self.threshold) == 0:
                self.score = 1
            else:
                self.score = float(parsed[4]) / self.threshold
        else:
            self.score = float(parsed[4])

        if len(self.threshold) == 0:
            self.threshold = 0
            self.warning = f"WARNING: {self.KO} does not have a KO threshold. All hits >0 will be included."

        if self.score >= float(self.threshold) * float(t_scale):
            self.threshold = float(self.threshold)
            self.passed = True
        else:
            self.passed = False

    def view(self):
        return [self.gene, self.KO, self.threshold, self.score, self.evalue, self.description]

    def result(self):
        return {'gene': self.gene,
                'annotation': self.KO,
                'score': self.score,
                'evalue': self.evalue}

    def __str__(self):
        return "<JAKomics KOFAM class>"


def run_kofam(faa_path, hal_path, temp_dir, ko_list, cpus=1, t_scale=1, score_as_ratio = False, echo=False, run=True):

    #temp_dir = 'KO_' + uuid.uuid4().hex

    command = f'exec_annotation --no-report-unannotated -k {ko_list} --tmp-dir {temp_dir} "{faa_path}" -T {t_scale} --cpu {int(cpus)} --profile {hal_path} -f detail-tsv'
    kofam_out = system_call(command, return_type="out", echo=echo, run=run)

    hits = []
    for line in kofam_out:
        if len(line) > 0 and not line.lstrip().startswith('#'):
            hits.append(KOFAM(line, t_scale, score_as_ratio = score_as_ratio))

    return hits


def parse_kofam_hits(run_kofam_out):
    '''
    Returns a dictionary of passed results with KO as key and list of kofam classes as value
    '''
    parsed = {}
    for hit in run_kofam_out:
        if hit.score_as_ratio:
            # print(db['DB_NAME'], genome.short_name, hit.view(), sep="\t")
            if hit.KO in parsed:
                parsed[hit.KO].append(hit)
            else:
                parsed[hit.KO] = [hit]
        elif hit.passed:
            # print(db['DB_NAME'], genome.short_name, hit.view(), sep="\t")
            if hit.KO in parsed:
                parsed[hit.KO].append(hit)
            else:
                parsed[hit.KO] = [hit]
    return parsed


def kofam_to_df(run_kofam_out):

    results = pd.DataFrame(columns=['LOCUS_TAG', 'KO', 'SCORE', 'THRESHOLD', 'EVALUE',
                                    'DESCRIPTION'])
    for hit in run_kofam_out:
        if hit.passed:
            results = results.append(
                pd.Series(data={'LOCUS_TAG': hit.gene,
                                'KO': hit.KO,
                                'SCORE': hit.score,
                                'THRESHOLD': hit.threshold,
                                'EVALUE': hit.evalue,
                                'DESCRIPTION': hit.description
                                }
                          ),
                ignore_index=True)
    return results
