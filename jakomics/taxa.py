class RDP:

    '''
    class for RDP fixrank taxonomy
    '''

    def __init__(self, line, threshold):
        self.RDP = line
        self.threshold = threshold
        self.ASV = line['ASV']

        self.D = self.assign_taxonomy(self.RDP['D_p'], self.RDP['D'], 'root')
        self.P = self.assign_taxonomy(self.RDP['P_p'], self.RDP['P'], self.D)
        self.C = self.assign_taxonomy(self.RDP['C_p'], self.RDP['C'], self.P)
        self.O = self.assign_taxonomy(self.RDP['O_p'], self.RDP['O'], self.C)
        self.F = self.assign_taxonomy(self.RDP['F_p'], self.RDP['F'], self.O)
        self.G = self.assign_taxonomy(self.RDP['G_p'], self.RDP['G'], self.F)

    def view(self):
        # print(f'RDP\t{self.ASV}\t{self.RDP["D"]}({self.RDP["D_p"]})\t{self.RDP["P"]}({self.RDP["P_p"]})\t{self.RDP["C"]}({self.RDP["C_p"]})\t{self.RDP["O"]}({self.RDP["O_p"]})\t{self.RDP["F"]}({self.RDP["F_p"]})\t{self.RDP["G"]}({self.RDP["G_p"]})')
        print(self.ASV, self.D, self.P, self.C, self.O, self.F, self.G, sep="\t")

    def assign_taxonomy(self, p, level, up_level):
        p = float(p.strip('%'))/100
        if up_level.startswith("unclassified"):
            return up_level
        elif p >= self.threshold:
            return level
        else:
            return 'unclassified_' + up_level

    def count_unclassified(self, counts):
        if self.D.startswith("unclassified"):
            counts['domain'] += 1
        if self.P.startswith("unclassified"):
            counts['phylum'] += 1
        if self.C.startswith("unclassified"):
            counts['class'] += 1
        if self.O.startswith("unclassified"):
            counts['order'] += 1
        if self.F.startswith("unclassified"):
            counts['family'] += 1
        if self.G.startswith("unclassified"):
            counts['genus'] += 1

        return counts
