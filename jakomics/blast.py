from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline, NcbimakeblastdbCommandline
import pandas as pd


class Blast:

    '''
    class for each column of line from blast = 8 or blast+ = 6
    '''

    def __init__(self, line, q, db):
        blastSplit = line.split("\t")

        self.parsed = blastSplit
        self.query = blastSplit[0]
        self.subject = blastSplit[1]
        self.percent = float(blastSplit[2])
        self.alignment_length = int(blastSplit[3])
        self.mismatches = int(blastSplit[4])
        self.gap_openings = int(blastSplit[5])
        self.query_start = int(blastSplit[6])
        self.query_end = int(blastSplit[7])
        self.subject_start = int(blastSplit[8])
        self.subject_end = int(blastSplit[9])
        self.eval = float(blastSplit[10])
        self.bit_score = float(blastSplit[11])
        self.query_file_path = q
        self.db_file_path = db

    def __str__(self):
        return f"<JAKomics BLAST Class>"

    def view(self):
        return [self.query, self.subject, self.percent, self.eval]

    def print_rough_result(self):
        print(self.query, self.subject, self.percent, self.eval, sep="\t")

    def print_full_result(self):
        print(f'{self.query} ({self.query_start}-{self.query_end}) hits {self.subject} ({self.subject_start}-{self.subject_end}) at {self.percent}% and {self.mismatches}/{self.gap_openings} mismatches/gaps (e-value: {self.eval}, score: {self.bit_score}).')
        # print(self.query, self.subject, self.percent, self.alignment_length, self.mismatches,
        # self.gap_openings, self.query_start + "-" + self.query_stop, sep="\t")

    def filter(self, e=10, b=1, p=35):
        passed = False
        if self.eval <= float(e):
            if self.bit_score >= float(b):
                if self.percent >= float(p):
                    passed = True

        return passed

    def result(self):
        return {'gene': self.query,
                'annotation': self.subject,
                'score': self.bit_score,
                'evalue': self.eval}

    def series(self):
        '''
        Return the results of a single blast hit as a pandas dataseries
        '''
        s = pd.Series(data={
            'qseqid': self.query,
            'sseqid': self.subject,
            'pident': self.percent,
            'length': self.alignment_length,
            'mismatch': self.mismatches,
            'gapopen': self.gap_openings,
            'qstart': self.query_start,
            'qend': self.query_end,
            'sstart': self.subject_start,
            'send': self.subject_end,
            'evalue': self.eval,
            'bitscore': self.bit_score
        })

        return s


def test():
    print("blast module loaded correctly")


def make_blast_db(type, db):
    make_blast_db_cl = NcbimakeblastdbCommandline(
        dbtype=type,
        input_file=db)
    make_blast_db_cl()


def run_blast(type, q, db, threads=1, e=0.001, make=False, return_query_results=True, echo=False):
    '''
    type = "prot" or "nucl"
    '''

    if make:
        make_blast_db(type, db)

    if type == 'prot':
        blast_type = NcbiblastpCommandline
    elif type == "nucl":
        blast_type = NcbiblastnCommandline

    blast_cline = blast_type(
        query=q,
        db=db,
        evalue=e,
        outfmt=6,
        num_threads=threads)

    if echo:
        print(blast_cline)

    stdout, stderr = blast_cline()
    raw_results = stdout.split("\n")

    if echo:
        print(stderr)

    results = {}
    if return_query_results:
        for line in [line.strip() for line in raw_results]:
            if len(line) > 0:
                hit = Blast(line, q, db)
                if hit.query in results:
                    results[hit.query].append(hit)
                else:
                    results[hit.query] = [hit]
    else:
        for line in [line.strip() for line in raw_results]:
            if len(line) > 0:
                hit = Blast(line, q, db)
                if hit.subject in results:
                    results[hit.subject].append(hit)
                else:
                    results[hit.subject] = [hit]

    return results


def blast_to_df(blast_results):
    '''
    pass in results from blast.run_blast (a dictionary) and return a
    nicely formatted dataframe suitable for printing.
    '''

    df = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

    for locus_tag, blast_hits in blast_results.items():
        for blast_hit in blast_hits:
            df = df.append(
                blast_hit.series(),
                ignore_index=True)
            df = pd.concat([df, blast_hit.series().to_frame().T],
                            ignore_index=True)

    return df




