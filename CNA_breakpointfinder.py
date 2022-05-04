import pandas as pd
import numpy as np
import sys

version = "1.0"

def Help():
    s = """CNA_breakpointfinder.py
    Compares the start coordinate of a CNA across each cell in the population to identify unique breakpoints, derived from CNAs.
Breakpointfinder.py can by run by the following input the in terminal:
  python CNA_breakpointfinder.py -b /path/to/bedfile.bed --lb 1.0 --ub 1.0 -o /path/to/outputPrefix
Required arguments are:
  -b,--bedfile              input bedfile,
  -li,--lowerInterval          sets the lower interval for CNA comparison (this will be multiplied by 10^6),
  -ui,--upperInterval          sets the upper interval for CNA comparison (this will be multiplied by 10^6),
  -o,--output               output prefix for bed file.

Author:
Alexander Malling Andersen
Email: Alexander_ma11@hotmail.dk"""

    print(s)
    print('Version ',version)
    exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
    pass


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

    a = arguments()
    a.lowerInterval = 1.0
    a.upperInterval = 1.0

    for i, argument in enumerate(args):

        if argument == '-b' or argument == '--bedfile':
            a.bedfilePath = str(args[i+1])

        elif argument == '-li' or argument == '--lowerInterval':
            a.lowerInterval = float(args[i+1])

        elif argument == '-ui' or argument == '--upperInterval':
            a.upperInterval = float(args[i+1])

        elif argument == '-o' or argument == '--output':
            a.outPrefix = str(args[i+1])

        elif argument == '-h' or argument == '--help' or argument == '--version' or argument == '-v':
            splashHelp()

    #check that required arguments are met
    if not ( ( hasattr( a, 'bedfilePath') or hasattr( a, 'lowerInterval') or hasattr( a, 'upperInterval') ) and  hasattr( a, 'outPrefix') ):
        splashHelp()
    return a

#-------------------------------------------------
def somy2var(somy):

    splitS = somy.split('-')
    return int(splitS[0])

#-------------------------------------------------
class SCNA:
    def __init__(self, chromosome, start, end, var):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.var = var

#-------------------------------------------------
class SingleCellSeq:
    def __init__(self):
        self.SCNAs = []
        self.name = ''

#-------------------------------------------------

#Create list of all cells, as objects, and their CNA features (chromosome, start-coordinate, end-coordinate, copy number)
args = parseArguments(sys.argv[1:])
lowerInterval = args.lowerInterval
upperInterval = args.upperInterval

filename = args.bedfilePath

f = open(filename,'r')

allCells = []
first = True
for line in f:

    if len(line.rstrip()) == 0:
        continue

    if line[0:5] == 'track':
        if not first:
            allCells.append(cell)
        first = False
        cell = SingleCellSeq()
        splitLine = line.rstrip().split('"')
        splitLine2 = splitLine[1].split(' ')
        cell.name = splitLine2[-1:][0][0:splitLine2[-1:][0].index('.')]

    else:
        splitLine = line.rstrip().split('\t')
        chromosome = splitLine[0]
        start = int(splitLine[1])
        end = int(splitLine[2])
        call = int(somy2var(splitLine[3]))
        cell.SCNAs.append(SCNA(chromosome,start,end,call))
allCells.append(cell)

filtered_SCNAs = {}
for cell in allCells:
    for s in cell.SCNAs:
        if (s.chromosome) not in filtered_SCNAs:
            filtered_SCNAs[(s.chromosome)] = []
        if (s.start) != 0:
            filtered_SCNAs[(s.chromosome)].append((s.start,s.end, s.var, cell.name))

#Find unique CNAs
print("There is a total of ", len(allCells)," cells in this bedfile. How many cells can a breakpoint be present in?")
number_cells = float(input())
fraction = number_cells/len(allCells)

a_uniqueCNAsBed = args.outPrefix + '_CellThreshold='+str(number_cells)+'_lowerInterval='+str(args.lowerInterval)+'mb'+'_upperInterval='+str(args.upperInterval)+'mb'+'.bed'
b_uniqueCNAsBed = open(a_uniqueCNAsBed,'w')

for chromosome in filtered_SCNAs:
    for i in filtered_SCNAs[chromosome]:
        overlaps = 0
        for j in filtered_SCNAs[chromosome]:
            if i[0] == j[0] and i[3] == j[3]:
                overlaps += 1
            elif j[0]-(lowerInterval*10**6) < i[0] < j[0]+(upperInterval*10**6):
                    overlaps += 1
            else:
                continue

        if float(overlaps)/len(allCells) <= fraction:
            b_uniqueCNAsBed.write(chromosome + ' ' + str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + '\n')

b_uniqueCNAsBed.close()
