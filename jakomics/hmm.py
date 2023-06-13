import os
import re
import pandas as pd
from jakomics.utilities import system_call, check_executable


def test():
    print("hmm module loaded correctly")


class HMM:

    '''
    class for a line of hmmer3 hmmsearch output
    '''

    def __init__(self, line):
        split = line.split()

        self.parsed = split
        self.gene = split[0]
        self.model = split[3].replace(".hmm", "")
        self.evalue = float(split[6])
        self.c_evalue = float(split[11])
        self.i_evalue = float(split[12])
        self.score = float(split[7])
        self.gene_length = int(split[2])
        self.model_length = int(split[5])
        self.gene_coordinates = (int(split[17]), int(split[18]))
        self.model_coordinates = (int(split[15]), int(split[16]))
        self.description = split[22]

    def align_length(self):
        return(self.model_coordinates[1] - self.model_coordinates[0] + 1)

    def hmm_coverage(self):
        return((self.model_coordinates[1] - self.model_coordinates[0] + 1) / self.model_length)

    def gene_coordinates_str(self):
        return str(self.gene_coordinates[0]) + "-" + str(self.gene_coordinates[1])

    def model_coordinates_str(self):
        return str(self.model_coordinates[0]) + "-" + str(self.model_coordinates[1])

    def model_coverage_str(self):
        return str(round(100*self.hmm_coverage(), 2)) + '%'

    def view(self):
        return [self.gene, self.model, self.evalue, self.score, self.c_evalue, self.i_evalue, self.gene_coordinates_str(), self.model_coordinates_str(), self.align_length(), self.model_coverage_str(), self.description]

    def write(self):
        return(self.gene + "\t" + self.model + "\t" + str(self.evalue) + "\t" + str(self.score) + "\t" + str(self.c_evalue) + "\t" + str(self.i_evalue) + "\t" + self.gene_coordinates_str() + "\t" + self.model_coordinates_str() + "\t" + str(self.align_length()) + "\t" + self.model_coverage_str() + "\t" + self.description + "\t" + "\n")

    def result(self):
        return {'gene': self.gene,
                'annotation': self.model,
                'score': self.score,
                'evalue': self.evalue}

    def __str__(self):
        return "<JAKomics HMM class>"


class CAZYME(HMM):

    def __str__(self):
        return "<JAKomics HMM/CAZYME class>"

    def pass_cazy(self):
        if self.evalue <= 1e-18 and self.hmm_coverage() >= 0.35:
            self.pass_qc = '1A'
        elif self.evalue <= 1e-15 and self.hmm_coverage() >= 0.35:
            self.pass_qc = '1B'
        elif self.align_length() >= 80 and self.evalue <= 1e-05:
            self.pass_qc = '2'
        elif self.evalue <= 1e-03 and self.hmm_coverage() >= 0.3:
            self.pass_qc = '3'
        else:
            self.pass_qc = '0'

        return self.pass_qc

    def assign_cazy_class(self):
        self.cazy_class = re.split(r'(\d+)', self.model)[0]

    def assign_substrate(self):
        family = self.model.split("_")[0]
        self.substrate = '-'
        if self.cazy_class == 'GH':
            if family in ['GH1', 'GH2', 'GH3', 'GH4', 'GH20', 'GH31']:
                self.substrate = 'Oligosaccharides'
            elif family in ['GH13', 'GH14', 'GH15', 'GH57', 'GH65']:
                self.substrate = 'Starch'
            elif family in ['GH5', 'GH6', 'GH8', 'GH9', 'GH12', 'GH44', 'GH45', 'GH48']:
                self.substrate = "Cellulose"
            elif family in ['GH10', 'GH11', 'GH30']:
                self.substrate = 'Xylan'
            elif family in ['GH66', 'GH70']:
                self.substrate = 'Dextran'
            elif family in ['GH32', 'GH68']:
                self.substrate = 'Fructan'
            elif family in ['GH16', 'GH26', 'GH28', 'GH39', 'GH43', 'GH53', 'GH67', 'GH78']:
                self.substrate = 'OtherPlantPolysaccharides'
            elif family in ['GH18', 'GH19', 'GH85']:
                self.substrate = 'Chitin'
            elif family in ['GH38', 'GH88', 'GH92', 'GH101']:
                self.substrate = 'OtherAnimalPolysaccharides'
            elif family in ['GH29', 'GH35', 'GH42', 'GH46', 'GH49', 'GH59', 'GH71', 'GH75', 'GH76', 'GH97', 'GH100', 'GH108']:
                self.substrate = 'Mixed'

    def series(self):
        '''
        Return the results of a single hmm hit as a pandas dataseries
        '''
        s = pd.Series(data={
            'LOCUS': self.gene,
            'HMM': self.model,
            'EVAL': self.evalue,
            'SCORE': self.score,
            'C-EVAL': self.c_evalue,
            'I-EVAL': self.i_evalue,
            'SEQ_COORDS': self.gene_coordinates_str(),
            'HMM_COORDS': self.model_coordinates_str(),
            'ALIGN_LENGTH': self.align_length(),
            'HMM_COVERAGE': self.model_coverage_str(),
            'QC_CODE': self.pass_qc,
            'CLASS': self.cazy_class,
            'SUBSTRATE': self.substrate
        })

        return s

    def view(self):
        print(self.gene, self.model, self.evalue, self.score, self.c_evalue, self.i_evalue, self.gene_coordinates_str(),
              self.model_coordinates_str(), self.align_length(), self.model_coverage_str(), self.pass_qc, self.cazy_class, self.substrate, sep="\t")

    def write(self):
        return(self.gene + "\t" + self.model + "\t" + str(self.evalue) + "\t" + str(self.score) + "\t" + str(self.c_evalue) + "\t" + str(self.i_evalue) + "\t" + self.gene_coordinates_str() + "\t" + self.model_coordinates_str() + "\t" + str(self.align_length()) + "\t" + self.model_coverage_str() + "\t" + str(self.pass_qc) + "\t" + self.cazy_class + "\t" + self.substrate + "\n")


def run_hmmsearch(path, log, raw, db, eval=0.001, score=10, cut_tc=False, echo=False, run=True):
    
    try:
        check_executable("hmmsearch")

        if cut_tc == True:
            command = f'hmmsearch -o {log} --cut_tc --domtblout {raw}'
        else:
            command = f'hmmsearch -o {log} -E {str(eval)} --domtblout {raw}'
        command += f' "{db}" "{path}"'

        return(system_call(command, echo=echo, run=run, return_type='err'))
    except:
        print(f"something went wrong with hmmsearch")

    
    


def parse_hmm_hits(file_path):
    lines = [line.strip() for line in open(file_path)]
    parsed = {}
    for line in lines:
        if not line.startswith("#"):
            hit = HMM(line)
            if hit.model in parsed:
                parsed[hit.model].append(hit)
            else:
                parsed[hit.model] = [hit]

    return parsed


def cazymes_to_df(raw_results, qc=['1A', '1B']):
    '''
    pass in a path to a raw hmm results file (output from run_hmmsearch), and parse for cazymes. Returns a
    nicely formatted dataframe
    '''

    rawResult = [line.strip() for line in open(raw_results)]

    results = pd.DataFrame(columns=['LOCUS', 'HMM', 'EVAL', 'SCORE', 'C-EVAL', 'I-EVAL', 'SEQ_COORDS',
                                    'HMM_COORDS', 'ALIGN_LENGTH', 'HMM_COVERAGE', 'QC_CODE', 'CLASS', 'SUBSTRATE'])

    for line in rawResult:
        if not line.startswith("#"):

            cazyme = CAZYME(line)
            if cazyme.pass_cazy() in qc:
                cazyme.assign_cazy_class()
                cazyme.assign_substrate()
                results = results.append(
                    cazyme.series(),
                    ignore_index=True)

    return results
