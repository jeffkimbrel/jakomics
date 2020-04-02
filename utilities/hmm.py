class HMM:

    '''
    class for each column of a hmmer3 output
    '''

    def __init__(self, line):
        split = line.split()

        self.gene = split[0]
        self.model = split[3].replace(".hmm", "")
        self.evalue = float(split[6])
        self.c_evalue = float(split[11])
        self.i_evalue = float(split[12])
        self.score = float(split[7])
        self.geneLength = int(split[2])
        self.modelLength = int(split[5])
        self.geneCoordinates = (int(split[17]), int(split[18]))
        self.modelCoordinates = (int(split[15]), int(split[16]))
        self.description = split[22]

        self.passQC = 0

    def alignLength(self):
        return(self.modelCoordinates[1] - self.modelCoordinates[0] + 1)

    def hmmCoverage(self):
        return((self.modelCoordinates[1] - self.modelCoordinates[0] + 1) / self.modelLength)

    def passCAZY(self):
        if self.evalue <= 1e-18 and self.hmmCoverage() >= 0.35:
            self.passQC = '1A'
        elif self.evalue <= 1e-15 and self.hmmCoverage() >= 0.35:
            self.passQC = '1B'
        elif self.alignLength() >= 80 and self.evalue <= 1e-05:
            self.passQC = '2'
        elif self.evalue <= 1e-03 and self.hmmCoverage() >= 0.3:
            self.passQC = '3'

        return(self.passQC)

    def view(self):
        geneCoordinates = str(self.geneCoordinates[0]) + "_" + str(self.geneCoordinates[1])
        modelCoordinates = str(self.modelCoordinates[0]) + "_" + str(self.modelCoordinates[1])
        modelCoverage = str(round(100*self.hmmCoverage(), 2)) + '%'

        print(self.gene, self.model, self.evalue, self.score, self.c_evalue, self.i_evalue, geneCoordinates,
              modelCoordinates, self.alignLength(), modelCoverage, self.description, self.passQC, sep="\t")

    def write(self):
        geneCoordinates = str(self.geneCoordinates[0]) + "_" + str(self.geneCoordinates[1])
        modelCoordinates = str(self.modelCoordinates[0]) + "_" + str(self.modelCoordinates[1])
        modelCoverage = str(round(100*self.hmmCoverage(), 2)) + '%'

        return(self.gene + "\t" + self.model + "\t" + str(self.evalue) + "\t" + str(self.score) + "\t" + str(self.c_evalue) + "\t" + str(self.i_evalue) + "\t" + str(geneCoordinates) + "\t" + str(modelCoordinates) + "\t" + str(self.alignLength()) + "\t" + str(modelCoverage) + "\t" + self.description + "\t" + str(self.passQC) + "\n")


def test():
    print("hmm module loaded correctly")
