#--------------------------------------------------------------------------------
# Copyright (c) 2020 Michael A. Boemo (mb915@cam.ac.uk)

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#--------------------------------------------------------------------------------


import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

VERSION="0.1.1"


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
    s = """clonalMasker.py: Uses control cells to filter clonal SCNAs from test cells.
To run clonalMasker.py, do:
  python clonalMasker.py -c /path/to/controlCells.bed -t /path/to/testCells.bed -f 0.01 -o /path/to/outputPrefix
Required arguments are:
  -c,--control              output directory which will be created,
  -t,--test                 output directory which will be created,
  -f,--fraction             mask clonal CNVs that occur in greater than this fraction of control cells,
  -o,--output               output prefix for bed and plot files.
Optional arguments are:
     --celllPloidy          set the ploidy of your cells
     --minOverlap           minimum fraction of overlap between two CNAs above which they are considered the same (default: 0.5),
     --largeGenes           bed file of large genes,
     --geneExpression       bed file of genes with expression,
     --replicationTiming    bed file of replication timing.
Written by Michael Boemo, Department of Pathology, University of Cambridge.
Please submit bug reports to mb915@cam.ac.uk."""

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
    a.minOverlap = 0.5
    a.useLargeGenes = False
    a.computeGeneExpression = False
    a.replicationTiming = False

    for i, argument in enumerate(args):

        if argument == '-c' or argument == '--control':
            a.controlPath = str(args[i+1])

        elif argument == '-t' or argument == '--test':
            a.testPath = str(args[i+1])

        elif argument == '-f' or argument == '--fraction':
            a.fraction = float(args[i+1])

        elif argument == '-o' or argument == '--output':
            a.outPrefix = str(args[i+1])

        elif argument == '--cellPloidy':
            a.cellPloidy = float(args[i+1])

        elif argument == '--minOverlap':
            a.minOverlap = float(args[i+1])

        elif argument == '--largeGenes':
            a.useLargeGenes = True
            a.genePath = str(args[i+1])

        elif argument == '--replicationTiming':
            a.replicationTiming = True
            a.replicationPath = str(args[i+1])

        elif argument == '--geneExpression':
            a.computeGeneExpression = True
            a.geneExpressionPath = str(args[i+1])

        elif argument == '-h' or argument == '--help' or argument == '--version' or argument == '-v':
            splashHelp()

    #check that required arguments are met
    if not ( ( hasattr( a, 'controlPath') or hasattr( a, 'testPath') or hasattr( a, 'fraction') ) and  hasattr( a, 'outPrefix') ):
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
def parseGenes(filename):
    f = open(filename,'r')
    genes_chr2bounds = {}
    for line in f:

        if len(line.rstrip()) == 0:
            continue

        splitLine = line.rstrip().split()
        chromosome = splitLine[0]
        lb = int(splitLine[1])
        ub = int(splitLine[2])
        if chromosome in genes_chr2bounds:
            genes_chr2bounds[chromosome].append((lb,ub))
        else:
            genes_chr2bounds[chromosome] = [(lb,ub)]

    f.close()
    return genes_chr2bounds


#-------------------------------------------------
def parseGeneExpression(filename):
    f = open(filename,'r')
    genes_chr2bounds = {}
    for line in f:

        if len(line.rstrip()) == 0:
            continue

        splitLine = line.rstrip().split()
        chromosome = splitLine[0]
        lb = int(splitLine[1])
        ub = int(splitLine[2])

        if splitLine[3] == '-':
            continue

        expression = float(splitLine[3])
        if chromosome in genes_chr2bounds:
            genes_chr2bounds[chromosome].append((lb,ub,expression))
        else:
            genes_chr2bounds[chromosome] = [(lb,ub,expression)]

    f.close()
    return genes_chr2bounds


#-------------------------------------------------
def parseReplicationTiming(filename):
    f = open(filename,'r')
    chrBin2timing = {}
    for line in f:

        if len(line.rstrip()) == 0:
            continue

        splitLine = line.rstrip().split()

        #pass over bins that don't have timing
        if len(splitLine) < 3:
            continue

        chromosome = splitLine[0]
        timingBin = int(splitLine[1])
        timing = int(splitLine[2])

        chrBin2timing[(chromosome,timingBin)] = timing

    f.close()
    return chrBin2timing

#-------------------------------------------------
def replicationWindows(cells_test,args,chrBin2timing):

    #write the header
    print('chromosome','CNA_start', 'CNA_end', 'breakpoint', 'timing')

    for cell in cells_test:

        for s in cell.SCNAs:

            lb = int(s.start)
            ub = int(s.end)

            for breakpoint in [lb,ub]:

                breakpointBin = int(breakpoint/1000000) + 1

                if (s.chromosome,breakpointBin) in chrBin2timing:

                    print(s.chromosome, s.start, s.end, breakpoint, chrBin2timing[(s.chromosome,breakpointBin)])


#-------------------------------------------------
def geneExpressionWindows(cells_test,args,genes_chr2bounds):

    windowSizes = range(100000,1000001,100000)

    #write a header
    print('chromosome','CNA_start', 'CNA_end', 'breakpoint', 'Side', " ".join([str(el) for el in windowSizes]))

    for cell in cells_test:

        for s in cell.SCNAs:

            if s.chromosome not in genes_chr2bounds:
                continue

            lb = int(s.start)
            ub = int(s.end)

            for breakpoint in [lb,ub]:

                windowSizes_expression = []
                windowSizes_expression_left = []
                windowSizes_expression_right = []

                for ws in windowSizes:

                    windowExpTotal = 0.
                    windowExpTotal_left = 0.
                    windowExpTotal_right = 0.

                    for gene in genes_chr2bounds[s.chromosome]:

                        if (breakpoint - ws < gene[0] < breakpoint + ws) and (breakpoint - ws < gene[1] < breakpoint + ws):

                            windowExpTotal += gene[2]

                        if (breakpoint - ws < gene[0] < breakpoint) and (breakpoint - ws < gene[1] < breakpoint):

                            windowExpTotal_left += gene[2]

                        if (breakpoint < gene[0] < breakpoint + ws) and (breakpoint < gene[1] < breakpoint + ws):

                            windowExpTotal_right += gene[2]

                    windowSizes_expression.append(windowExpTotal)
                    windowSizes_expression_left.append(windowExpTotal_left)
                    windowSizes_expression_right.append(windowExpTotal_right)

                print(s.chromosome, s.start, s.end, breakpoint, 'Total', " ".join([str(el) for el in windowSizes_expression]))
                print(s.chromosome, s.start, s.end, breakpoint, 'Left', " ".join([str(el) for el in windowSizes_expression_left]))
                print(s.chromosome, s.start, s.end, breakpoint, 'Right', " ".join([str(el) for el in windowSizes_expression_right]))


#-------------------------------------------------
def plotGeneHistogram(cells_test,args,genes_chr2bounds):

    print('Chromosome', 'CopyNumber', 'CNAstart', '', 'CNAend', 'directionalDist(kb)')

    directionalDistances = []
    absDistances = []
    for cell in cells_test:

        for s in cell.SCNAs:

            if s.chromosome not in genes_chr2bounds:
                continue

            lb = int(s.start)
            ub = int(s.end)

            minDist = 10000000000
            minDist_left = 10000000000
            minDist_right = 10000000000
            found = False
            found_left = False
            found_right = False

            for gene in genes_chr2bounds[s.chromosome]:

                d1 = lb - gene[0]
                d2 = lb - gene[1]
                delta = [d1,d2]
                idx = np.argmin([abs(d1),abs(d2)])

                if abs(delta[idx]) < abs(minDist_left):
                    minDist_left = delta[idx]
                    found_left = True

            for gene in genes_chr2bounds[s.chromosome]:

                d3 = ub - gene[0]
                d4 = ub - gene[1]
                delta = [d3,d4]
                idx = np.argmin([abs(d3),abs(d4)])

                if abs(delta[idx]) < abs(minDist_right):
                    minDist_right = delta[idx]
                    found_right = True

            if found_left and found_right:
                print(s.chromosome, s.var, s.start, minDist_left, s.end, minDist_right)

            for gene in genes_chr2bounds[s.chromosome]:

                d1 = lb - gene[0]
                d2 = lb - gene[1]
                d3 = ub - gene[0]
                d4 = ub - gene[1]
                delta = [d1,d2,d3,d4]
                idx = np.argmin([abs(d1),abs(d2),abs(d3),abs(d4)])

                if abs(delta[idx]) < minDist:
                    minDist = abs(delta[idx])
                    directionalDist = delta[idx]
                    found = True

            if found:
                directionalDistances.append(directionalDist)
                absDistances.append(minDist)

            if not found:
                print('error')

    plt.figure()
    plt.hist(np.array(directionalDistances)/1000,50)
    plt.xlabel('Distance from Nearest Large Gene (kb)')
    plt.ylabel('Count')
    plt.savefig(args.outPrefix + '_distanceFromLargeGenes_fraction_'+str(clonalThreshold)+'_minOverlap_'+str(args.minOverlap)+'.pdf')
    plt.close()
    print('Signed mean distance from nearest large gene (bp):',np.mean(directionalDistances))
    print('Signed median distance from nearest large gene (bp):',np.median(directionalDistances))
    print('Absolute mean distance from nearest large gene (bp):',np.mean(absDistances))
    print('Absolute median distance from nearest large gene (bp):',np.median(absDistances))


#-------------------------------------------------
def parseBed(filename,clonalSCNAs,minOverlap):

    f = open(filename,'r')
    #print filename
    allCells = []
    first = True
    for line in f:

        if len(line.rstrip()) == 0:
            continue

        if line[0:5] == 'track':

            if not first:
                allCells.append(cell)
                print(cell)
            first = False
            cell = SingleCellSeq()
            print(cell)
            splitLine = line.rstrip().split('"')
            splitLine2 = splitLine[1].split(' ')
            cell.name = splitLine2[-1:][0][0:splitLine2[-1:][0].index('.')]
        else:

            splitLine = line.rstrip().split('\t')
            chromosome = splitLine[0]
            start = int(splitLine[1])
            end = int(splitLine[2])
            call = somy2var(splitLine[3])

            #print chromosome, start, end, call
            if call != cellPloidy:
                if (chromosome,call) not in clonalSCNAs:
                    cell.SCNAs.append(SCNA(chromosome,start,end,call))
                else:
                    isClonal = False
                    for clonalPair in clonalSCNAs[(chromosome,call)]:
                        overlap = max(0, min(clonalPair[1],end) - max(clonalPair[0],start))
                        if float(overlap) / abs(end - start) > minOverlap and float(overlap) / abs(clonalPair[1] - clonalPair[0]) > minOverlap:
                            isClonal = True
                            break
                    if not isClonal:
                        cell.SCNAs.append(SCNA(chromosome,start,end,call))

    allCells.append(cell)
    return allCells


#-------------------------------------------------
def countRes(res):
    count = 0
    for c in range(1,23):
        for i in range(0,chr2len('chr'+str(c)),res):
            count += 1
    for i in range(0,chr2len('chrX'),res):
        count += 1
    return count


#-------------------------------------------------
def columnCmap(resolution,palette):

    row = []
    for c in range(1,23):
        for i in range(0,chr2len('chr'+str(c)),resolution):
            row.append(palette[c])
    for i in range(0,chr2len('chrX'),resolution):
        row.append(palette[c])
    return row


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
def cellToRow(cell, resolution):

    row = []
    for c in range(1,23):
        for i in range(0,chr2len('chr'+str(c)),resolution):

            inSCNA = False
            for s in cell.SCNAs:
                if s.start < i and i < s.end and s.chromosome == 'chr'+str(c):
                    row.append(s.var)
                    inSCNA = True
                    break
            if not inSCNA:
                row.append(2)
    return row


#MAIN---------------------------------------------
args = parseArguments(sys.argv[1:])
clonalThreshold = args.fraction
cellPloidy = args.cellPloidy

#make the data
cells_control = parseBed(args.controlPath, [], args.minOverlap)
cells_test = parseBed(args.testPath, [], args.minOverlap)

#open output files
f_clonalBed = args.outPrefix + '_clonal_fraction_'+str(clonalThreshold)+'_minOverlap_'+str(args.minOverlap)+'.bed'
f_nonClonalBed = args.outPrefix + '_nonClonal_fraction_'+str(clonalThreshold)+'_minOverlap_'+str(args.minOverlap)+'.bed'
fh_clonalBed = open(f_clonalBed,'w')
fh_nonClonalBed = open(f_nonClonalBed,'w')

allSCNAs = {}
for cell in cells_control:
    for s in cell.SCNAs:
        if (s.chromosome,s.var) not in allSCNAs:
            allSCNAs[(s.chromosome,s.var)] = []
        allSCNAs[(s.chromosome,s.var)].append((s.start,s.end,cell.name))

clonalSCNAs = {}
for key in allSCNAs:
    for i, s1 in enumerate(allSCNAs[key]):
        overlaps = 0
        for j in range(0,len(allSCNAs[key])):

            overlap = max(0, min(allSCNAs[key][j][1],s1[1]) - max(allSCNAs[key][j][0],s1[0]))
            if float(overlap) / (s1[1] - s1[0]) > args.minOverlap:
                overlaps += 1

        if float(overlaps)/len(cells_control) > clonalThreshold:

            fh_clonalBed.write(key[0] + ' ' + str(s1[0]) + ' ' + str(s1[1]) + ' ' + str(key[1]) + ' ' + str(s1[2]) + '\n')

            if (key[0],key[1]) not in clonalSCNAs:
                clonalSCNAs[(key[0],key[1])] = [(s1[0],s1[1])]
            else:
                clonalSCNAs[(key[0],key[1])].append((s1[0],s1[1]))

cells_control = parseBed(args.controlPath, clonalSCNAs, args.minOverlap)
cells_test = parseBed(args.testPath, clonalSCNAs, args.minOverlap)


#make the chromosomes
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

plt.xlim(-10000000,end+10000000)
plt.ylim(-15000000,15000000)


chr2rectangles = {}
for c in chromosomeNames:
    chr2rectangles[c] = []

matrix = []
resolution = 2000000 #2 Mb
for cell in cells_test:
    row = np.array(cellToRow(cell, resolution))
    matrix.append(row)
    for s in cell.SCNAs:

        offset = chr2offsetStart[s.chromosome]

        fh_nonClonalBed.write(s.chromosome + ' ' + str(s.start) + ' ' + str(s.end) + ' ' + str(s.var) + ' ' + str(cell.name) + '\n')
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

fh_clonalBed.close()
fh_nonClonalBed.close()

args.outPrefix + '_nonClonal_fraction_'+str(clonalThreshold)+'_minOverlap_'+str(args.minOverlap)+'.bed'

plt.axis("off")
plt.savefig(args.outPrefix + '_pileup_fraction_'+str(clonalThreshold)+'_minOverlap_'+str(args.minOverlap)+'.pdf')

plt.figure()
colPalette = sns.color_palette("dark",n_colors=23)
colColours = columnCmap(resolution,colPalette)
p = sns.clustermap(np.array(matrix),col_colors=colColours,metric='euclidean',method='complete',row_cluster=True,col_cluster=False,cmap=sns.diverging_palette(10,240, as_cmap=True,center="light"),center=2.0,vmin=0,vmax=4)
plt.savefig(args.outPrefix + '_heatmap_fraction_'+str(clonalThreshold)+'.pdf')
plt.close()

if args.useLargeGenes:
    genes_chr2bounds = parseGenes(args.genePath)
    plotGeneHistogram(cells_test,args,genes_chr2bounds)

if args.computeGeneExpression:
    genes_chr2bounds = parseGeneExpression(args.geneExpressionPath)
    geneExpressionWindows(cells_test,args,genes_chr2bounds)

if args.replicationTiming:
    chrBin2timing = parseReplicationTiming(args.replicationPath)
    replicationWindows(cells_test,args,chrBin2timing)
