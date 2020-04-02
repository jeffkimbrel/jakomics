from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline, NcbimakeblastdbCommandline


class Blast:

    '''
    class for each column of line from blast = 8 or blast+ = 6
    '''

    def __init__(self, line):
        blastSplit = line.split("\t")

        self.query = blastSplit[0]
        self.subject = blastSplit[1]
        self.id = float(blastSplit[2])
        self.alignment_length = int(blastSplit[3])
        self.mismatches = int(blastSplit[4])
        self.gap_openings = int(blastSplit[5])
        self.query_start = int(blastSplit[6])
        self.query_end = int(blastSplit[7])
        self.subject_start = int(blastSplit[8])
        self.subject_end = int(blastSplit[9])
        self.eval = float(blastSplit[10])
        self.bit_score = float(blastSplit[11])

    def print_rough_result(self):
        print(self.query, self.subject, self.id, self.eval, sep="\t")

    def print_full_result(self):
        print(f'{self.query} ({self.query_start}-{self.query_end}) hits {self.subject} ({self.subject_start}-{self.subject_end}) at {self.id}% and {self.mismatches}/{self.gap_openings} mismatches/gaps (e-value: {self.eval}, score: {self.bit_score}).')
        # print(self.query, self.subject, self.id, self.alignment_length, self.mismatches,
        # self.gap_openings, self.query_start + "-" + self.query_stop, sep="\t")


def test():
    print("blast module loaded correctly")


def make_blast_db(type, db):
    make_blast_db_cl = NcbimakeblastdbCommandline(
        dbtype=type,
        input_file=db)
    make_blast_db_cl()


def do_blast(type, q, db, threads=8, e=0.001):
    '''
    type = "prot" or "nucl"
    '''

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

    stdout, stderr = blast_cline()
    raw_results = stdout.split("\n")
    results = {}
    for line in [line.strip() for line in raw_results]:
        if len(line) > 0:
            # print(line)
            hit = Blast(line)
            if hit.query in results:
                results[hit.query].append(hit)
            else:
                results[hit.query] = [hit]

    return results
