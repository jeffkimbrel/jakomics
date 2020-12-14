import re


class GENE:

    def __init__(self, id):
        self.id = id

    def __str__(self):
        return "<JAKomics GENE class>"

    def parse_gbk_location(self, location):

        split = re.split('\[|:|\]|\(|\)', str(location))
        strand = split[4]
        start = ""
        stop = ""

        split[1] = split[1].replace("<", "")
        split[2] = split[2].replace("<", "")
        split[1] = split[1].replace(">", "")
        split[2] = split[2].replace(">", "")

        if strand == '-':
            start = split[2]
            stop = str(int(split[1]) + 1)
        else:
            start = str(int(split[1]) + 1)
            stop = split[2]

        self.start = int(start)
        self.stop = int(stop)
        self.strand = strand

    def result(self, type="product"):
        '''
        method to get results suitable for gator
        '''

        return {'gene': self.id,
                'annotation': getattr(self, type, None),
                'score': getattr(self, "score", None),
                'evalue': getattr(self, "evalue", None)}
