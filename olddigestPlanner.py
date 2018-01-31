import re
from copy import copy
import operator
import math



class Enzyme:
    def __init__(self, name, sequence, sites):
        self.name = name
        self.sequence = sequence
        self.sites = sites

    def second_largest(self, numbers):
        count = 0
        m1 = m2 = float('-inf')
        for x in numbers:
            count += 1
            if x > m2:
                if x >= m1:
                    m1, m2 = x, m1
                else:
                    m2 = x
        return m2 if count >= 2 else None

    def pattern(self):
        self.numberBands = len(self.sites)
        self.bandSizes = []

        if self.numberBands > 0:
            tempSites = copy(self.sites)
            for cut in self.sites:
                if len(tempSites) > 2:
                    self.bandSizes.append(max(tempSites) - self.second_largest(tempSites))
                    tempSites.remove(max(tempSites))
                elif len(tempSites) == 2:
                    self.bandSizes.append(max(tempSites) - min(tempSites))
                    tempSites.remove(max(tempSites))
                else:
                    self.bandSizes.append(len(plasmid_sequence) - sum(self.bandSizes))

        self.bandSizes.sort(reverse = True)

    def select(self):
        self.cost = 0
        self.spacing = []

        # Band number
        if self.numberBands > 6:
            self.cost += 5 * self.numberBands
        elif self.numberBands > 3:
            self.cost += self.numberBands
        elif self.numberBands == 1 or self.numberBands == 0:
            self.cost += 250
        elif self.numberBands == 3:
            self.cost += -1

        # Band sizing
        for size in self.bandSizes:
            if size > 10000:
                self.cost += 100
            elif size > 7000:
                self.cost += size/500
            elif size < 500:
                self.cost += 100
            elif size < 1000:
                self.cost += (1000-size)/15

        # Band spacing
        tempBand = copy(self.bandSizes)
        for band in self.bandSizes:
            if len(tempBand) > 2:
                self.spacing.append(max(tempBand) - self.second_largest(tempBand))

                if band > 6500 and self.spacing[-1] < 1000:
                    self.cost += 100
                elif band > 3500 and self.spacing[-1] < 1000:
                    self.cost += self.spacing[-1]/25
                elif band > 1500 and self.spacing[-1] < 100:
                    self.cost += 50

                tempBand.remove(max(tempBand))
            elif len(tempBand) == 2:
                self.spacing.append(max(tempBand) - min(tempBand))
                tempBand.remove(max(tempBand))

                if band > 6500 and self.spacing[-1] < 1000:
                    self.cost += 100
                elif band > 3500 and self.spacing[-1] < 1000:
                    self.cost += self.spacing[-1]/25
                elif self.spacing[-1] < 100:
                    self.cost += 50

    def read_out(self):

        self.bandOutput = self.round_up(self.bandSizes)

        self.output = self.name, self.cost, self.bandOutput

        return self.output

    def round_up(self, bands):

        self.roundedBands = []
        for band in bands:
            self.roundedBands.append(int(round(band, -2)))

        return self.roundedBands


with open('plasmid.txt','r') as f_open:
    plasmid_sequence = f_open.read()

plasmid_sequence = plasmid_sequence.upper()

enzdic = {'HindIII': 'AAGCTT', 'BamHI': 'GGATCC', 'XhoI': 'CTCGAC', 'XmnI': 'GAANNNNTTC', 'XbaI': 'TCTAGA', 'AscI': 'GGCGCGCC', 'KpnI': 'GGTACC', 'PmeI': 'GTTTAAAC', 'ApaI': 'GGGCCC', 'SbfI': 'CCTGCAGG', 'EcoRV': 'GATATC', 'SwaI': 'ATTTAAAT', 'SnaBI': 'TACGTA', 'MluI': 'ACGCGT', 'NarI': 'GGCGCC', 'EcoRI': 'GAATTC', 'BseRI': 'GAGGAG', 'PstI': 'CTGCAG', 'NotI': 'GCGGCCGC', 'SphI': 'GCATGC', 'NdeI': 'CATATG', 'DpnI': 'GATC', 'BclI': 'TGATCA', 'HpaI': 'GTTAAC', 'BsrGI': 'TGTACA', 'AflII': 'CTTAAG', 'XmaI': 'CCCGGG', 'BbsI': 'GAAGAC', 'NcoI': 'CCATGG', 'BspEI': 'TCCGGA', 'AgeI': 'ACCGGT', 'PflMI': 'CCANNNNNTGG', 'SmaI': 'CCCGGG', 'BglII': 'AGATCT', 'SalI': 'GTCGAC', 'ClaI': 'ATCGAT', 'PshAI': 'GACNNNNGTC', 'NsiI': 'ATGCAT', 'PspOMI': 'GGGCCC', 'SacI': 'GAGCTC', 'FspI': 'TGCGCA', 'BsaI': 'GGTCTCN', 'SacII': 'CCGCGG', 'PacI': 'TTAATTAA', 'MfeI': 'CAATTG'}


outputArray = []

for key in enzdic.keys():
    sites = [m.start() for m in re.finditer(enzdic[key], plasmid_sequence)]
    enzyme = Enzyme(key, enzdic[key], sites)
    enzyme.pattern()
    enzyme.select()

    outputArray.append(enzyme.read_out())

outputArray.sort(key=operator.itemgetter(1))

print '------------------------------'
for item in outputArray[0:8]:
    print item[0], item[2], item[1]
print '------------------------------'
