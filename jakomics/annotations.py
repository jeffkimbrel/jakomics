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

    def return_genes(self):
        blast_genes = []
        if hasattr(self, 'blast_hits'):
            for hit in self.blast_hits:
                blast_genes.append(hit.query)

        kofam_genes = []
        if hasattr(self, 'kofam_hits'):
            for kofam in self.kofam_hits:
                kofam_genes.append(kofam.gene)

        return list(set(kofam_genes + blast_genes))
