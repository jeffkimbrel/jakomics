from jakomics import hmm
import re
from Bio import SwissProt
import gzip

class Annotation:
    def __init__(self, gene, annotations_template):
        self.gene = gene

        for key in annotations_template:
            setattr(self, key, annotations_template[key])

    def add_kofam(self, kofam_results):

        if hasattr(self, 'kofam_hits') == False:
            self.kofam_hits = []

        if hasattr(self, 'KO'):
            for kofam in kofam_results:
                if kofam.passed:
                    if kofam.KO in self.KO:
                        self.kofam_hits.append(kofam)

    def add_blast(self, blast_results):

        if hasattr(self, 'blast_hits') == False:
            self.blast_hits = []

        if hasattr(self, 'BLASTP'):
            for gene in blast_results:
                for hit in blast_results[gene]:
                    if hit.subject in self.BLASTP:
                        if hit.filter(e=1e-15):
                            self.blast_hits.append(hit)

    def add_hmm(self, hmm_results):

        if hasattr(self, 'hmm_hits') == False:
            self.hmm_hits = []

        if hasattr(self, 'HMM'):
            for line in hmm_results:
                if not line.startswith("#"):
                    hit = hmm.HMM(line)
                    if hit.model in self.HMM:
                        self.hmm_hits.append(hit)

    def return_genes(self):
        blast_genes = []
        if hasattr(self, 'blast_hits'):
            for hit in self.blast_hits:
                blast_genes.append(hit.query)

        kofam_genes = []
        if hasattr(self, 'kofam_hits'):
            for kofam in self.kofam_hits:
                kofam_genes.append(kofam.gene)

        hmm_genes = []
        if hasattr(self, 'hmm_hits'):
            for hmm in self.hmm_hits:
                hmm_genes.append(hmm.gene)

        return list(set(kofam_genes + blast_genes + hmm_genes))


def find_ec(x):
    '''
    regex to find EC numbers in a given string. Needs to start with EC, and then have 
    either ":" or "=". The last field can be a dash, or contain an "n".

    returns a list of hits, empty if none
    '''

    pattern = "([Ee][Cc][=:][\d+].[\d+].[\d+].[n]{0,1}[\d+-])"
    match = re.findall(pattern, x)
    return match

def parse_evidence_fields(x):
    '''
    function to parse evidence codes from swiss prot dat files with the format
    "Evidence={ECO:0000269|PubMed:18375606}"

    returns a dereplicated list of ECO/pubmed pairs
    '''

    pattern = "Evidence=\{(.*?)\}"
    match = re.findall(pattern, x)
    
    r = []
    for m in match:
        r += re.split("[ ,]+", m)

    return list(set(r))

def get_swissprot_records(dat_path):

    # https://biopython.org/docs/1.78/api/Bio.SwissProt.html

    with gzip.open(dat_path) as handle:
        records = SwissProt.parse(handle)

        return records

# test

if __name__ == '__main__':
    query = "BI:77896;   Evidence={ECO:0000269|PubMed:11572992, ECO:0000269|PubMed:29281266,   ECO:0000269|PubMed:30304478, ECO:0000269|PubMed:7592783}; PhysiologicalDirection=left-to-right; Xref=Rhea:RHEA:31576;   Evidence={ECO:0000269|PubMed:11572992};', 'CATALYTIC ACTIVITY: Reaction=H2O + O(6)-methyl-dGTP"
    
    j = parse_evidence_fields(query)
    print("---")
    print(j)