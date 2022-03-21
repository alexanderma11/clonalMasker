import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
    s = """Pileup_plot.py: Creates pileup graphs, inspired from clonalMasker.py, by Michael Boemo.
To plot your CNAs, please do:
  python Pileup_plots.py -b /path/to/BedFile_with_CNAs -o /path/to/outputPrefix

Required arguments are:
  -b,--bedFile              input file, containing CNAs for plotting: Format must be "Chr, Start, End, CNA-Ploidy"
  -c,--cellPloidy          set the modal ploidy of your cells: Default set as Diploid
  -o,--output               output prefix for bed and plot files.

Written by Alexander Malling Andersen, Department of Cancer Genomics and Computational Biology.
The framework was inspired by Michael Boemo, Department of Pathology, University of Cambridge.
Please submit bug reports to alexander_ma11@hotmail.dk"""

    print(s)
    print('Version ',VERSION)
    exit(0)

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
    pass

#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

    a = arguments()
    a.cellPloidy = 2

    for i, argument in enumerate(args):

        if argument == '-b' or argument == '--bedFile':
            a.bedFile = str(args[i+1])

        elif argument == '-c' or argument == '--cellPloidy':
            a.cellPloidy = float(args[i+1])

        elif argument == '-o' or argument == '--output':
            a.outPrefix = str(args[i+1])

        elif argument == '-h' or argument == '--help' or argument == '--version' or argument == '-v':
            splashHelp()

    #check that required arguments are met
    if not ( ( hasattr( a, 'bedFile') or hasattr( a, 'cellPloidy') ) and  hasattr( a, 'outPrefix') ):
        splashHelp()
    return a

#-------------------------------------------------
def chr2len(chromosome):

    if chromosome == 'chr1':
        return 248956422
    elif chromosome == 'chr2':
        return     242193529
    elif chromosome == 'chr3':
        return     198295559
    elif chromosome == 'chr4':
        return     190214555
    elif chromosome == 'chr5':
        return     181538259
    elif chromosome == 'chr6':
        return     170805979
    elif chromosome == 'chr7':
        return     159345973
    elif chromosome == 'chr8':
        return     145138636
    elif chromosome == 'chr9':
        return     138394717
    elif chromosome == 'chr10':
        return     133797422
    elif chromosome == 'chr11':
        return     135086622
    elif chromosome == 'chr12':
        return     133275309
    elif chromosome == 'chr13':
        return     114364328
    elif chromosome == 'chr14':
        return     107043718
    elif chromosome == 'chr15':
        return     101991189
    elif chromosome == 'chr16':
        return     90338345
    elif chromosome == 'chr17':
        return     83257441
    elif chromosome == 'chr18':
        return     80373285
    elif chromosome == 'chr19':
        return     58617616
    elif chromosome == 'chr20':
        return     64444167
    elif chromosome == 'chr21':
        return     46709983
    elif chromosome == 'chr22':
        return     50818468
    elif chromosome == 'chrX':
        return     156040895
    else:
        print('chromosome not recognised')

#-------------------------------------------------
def chr2centromere(chromosome):

    if chromosome == 'chr1':
        return 121700000, 125100000
    elif chromosome == 'chr2':
        return     91800000, 96000000
    elif chromosome == 'chr3':
        return     87800000, 94000000
    elif chromosome == 'chr4':
        return     48200000, 51800000
    elif chromosome == 'chr5':
        return     46100000, 51400000
    elif chromosome == 'chr6':
        return     58500000, 62600000
    elif chromosome == 'chr7':
        return     58100000, 62100000
    elif chromosome == 'chr8':
        return     43200000, 47200000
    elif chromosome == 'chr9':
        return     42200000, 45500000
    elif chromosome == 'chr10':
        return     38000000, 41600000
    elif chromosome == 'chr11':
        return     51000000, 55800000
    elif chromosome == 'chr12':
        return     33200000, 37800000
    elif chromosome == 'chr13':
        return     16500000, 18900000
    elif chromosome == 'chr14':
        return     16100000, 18200000
    elif chromosome == 'chr15':
        return     17500000, 20500000
    elif chromosome == 'chr16':
        return     35300000, 38400000
    elif chromosome == 'chr17':
        return     22700000, 27400000
    elif chromosome == 'chr18':
        return     15400000, 21500000
    elif chromosome == 'chr19':
        return     24200000, 28100000
    elif chromosome == 'chr20':
        return     25700000, 30400000
    elif chromosome == 'chr21':
        return     10900000, 13000000
    elif chromosome == 'chr22':
        return     13700000, 17400000
    elif chromosome == 'chrX':
        return     58100000, 63800000
    else:
        print('chromosome not recognised')

#-------------------------------------------------
def rectanglesOverlap(rec1,rec2):

    rec1_x1 = rec1[0]
    rec1_x2 = rec1[0] + rec1[2]
    rec1_y1 = rec1[1]
    rec1_y2 = rec1[1] + rec1[3]

    rec2_x1 = rec2[0]
    rec2_x2 = rec2[0] + rec2[2]
    rec2_y1 = rec2[1]
    rec2_y2 = rec2[1] + rec2[3]

    overlapX = (rec2_x1 <= rec1_x1 and rec1_x1 <= rec2_x2) or (rec2_x1 <= rec1_x2 and rec1_x2 <= rec2_x2) or (rec1_x1 <= rec2_x1 and rec2_x1 <= rec1_x2) or (rec1_x1 <= rec2_x2 and rec2_x2 <= rec1_x2)
    overlapY = (rec2_y1 <= rec1_y1 and rec1_y1 <= rec2_y2) or (rec2_y1 <= rec1_y2 and rec1_y2 <= rec2_y2) or (rec1_y1 <= rec2_y1 and rec2_y1 <= rec1_y2) or (rec1_y1 <= rec2_y2 and rec2_y2 <= rec1_y2)

    return overlapX and overlapY

#-------------------------------------------------
def rectanglesOverlapSearch(rec1,allRecs):

    foundOverlap = False
    for r in allRecs:
        if rectanglesOverlap(rec1,r):
            foundOverlap = True
    return foundOverlap

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
def somy2var(somy):

    splitS = somy.split('-')
    return int(splitS[0])

#-------------------------------------------------
args = parseArguments(sys.argv[1:])
CNA_bed = args.bedFile
cellPloidy = args.cellPloidy

plt.figure()
chr2offsetStart = {}
startIdx = 0

chromosomeNames = []
for c in range(1,23):
    chromosomeNames.append('chr'+str(c))

chromosomeNames.append('chrX')
colPalette = sns.color_palette("dark",n_colors=23)

for c in chromosomeNames:

    start = startIdx
    end = start + chr2len(c)

    cenStart, cenEnd = chr2centromere(c)
    cenStart += start
    cenEnd += start

    chr2offsetStart[c] = start

    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((start,-500000), chr2len(c), 1000000, alpha=1,facecolor='silver',edgecolor=None,linewidth = 0.01))

    currentAxis.add_patch(Rectangle((cenStart,-500000), cenStart-cenEnd, 1000000, alpha=1,facecolor='darkorange',edgecolor='k',linewidth = 0))

    startIdx = end + 10000000
    if len(c) == 4:
        plt.text(start, -500000, c[-1], fontsize = 3)
    else:
        plt.text(start, -500000, c[-2:], fontsize = 3)
plt.xlim(-10000000,end+10000000)
plt.ylim(-15000000,15000000)

chr2rectangles = {}
for c in chromosomeNames:
    chr2rectangles[c] = []

f = open(CNA_bed,'r')

allCNAs = []
for line in f:
    if len(line.rstrip()) == 0:
        continue

    else:
        cell = SingleCellSeq()
        splitLine = line.rstrip().split('\t')
        chromosome = splitLine[0]
        start = int(splitLine[1])
        end = int(splitLine[2])
        if '-' in splitLine[3]:
            call = somy2var(splitLine[3])
        else:
            call = int(splitLine[3])
        cell.name = splitLine[4]

        cell.SCNAs.append(SCNA(chromosome,start,end,call))
    allCNAs.append(cell)

for cell in allCNAs:

    for s in cell.SCNAs:

        offset = chr2offsetStart[s.chromosome]

        cna_colour = ""
        if s.var > cellPloidy+2:
            cna_colour = "darkred"
        elif s.var == cellPloidy+2:
            cna_colour = "red"
        elif s.var == cellPloidy+1:
            cna_colour = "lightcoral"
        elif s.var == cellPloidy-1:
            cna_colour = "lightsteelblue"
        elif s.var == cellPloidy-2:
            cna_colour = "blue"
        elif s.var < cellPloidy-2:
            cna_colour = "darkblue"
        else:
            continue

        if s.var > cellPloidy:
            yOffset = 1000000
            rect = [offset + s.start, yOffset, s.end-s.start, 200000]
            search = rectanglesOverlapSearch(rect,chr2rectangles[s.chromosome])
            while search:
                yOffset += 500000
                rect = [offset + s.start, yOffset, s.end-s.start, 200000]
                search = rectanglesOverlapSearch(rect,chr2rectangles[s.chromosome])
            chr2rectangles[s.chromosome].append(rect)
            currentAxis.add_patch(Rectangle((rect[0],rect[1]), rect[2], rect[3], alpha=0.5,facecolor=cna_colour,edgecolor=None))
        else:
            yOffset = -1000000
            rect = [offset + s.start, yOffset, s.end-s.start, 200000]
            search = rectanglesOverlapSearch(rect,chr2rectangles[s.chromosome])
            while search:
                yOffset -= 500000
                rect = [offset + s.start, yOffset, s.end-s.start, 200000]
                search = rectanglesOverlapSearch(rect,chr2rectangles[s.chromosome])
            chr2rectangles[s.chromosome].append(rect)
            currentAxis.add_patch(Rectangle((rect[0],rect[1]), rect[2], rect[3], alpha=0.5, facecolor=cna_colour,edgecolor=None))

plt.axis("off")
plt.savefig(args.outPrefix + '_pileup_graph'+'.pdf')
