import ROOT
import math
import optparse
import os, sys
from syncUtils import *
from operator import attrgetter


# define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--input', dest='inFile', type='string', default="ZZ4lAnalysis.root",    help='input file')
    parser.add_option('-o', '--output', dest='outFile', type='string', default="CJLSTevents.txt",    help='output sync file')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if not "." in opt.outFile:
        print "Please use an extension for the output file (e.g. \".txt\")"
        sys.exit()





def loop():

    inFileName = opt.inFile
    outFileName = opt.outFile
    print "Processing file: ",inFileName

    cands = []
    totCounter = 0
    chanCounter = {}

    hfile = ROOT.TFile(inFileName)
    hCounters = hfile.Get('ZZTree/Counters')

    tree = ROOT.TChain()
    tree.Add(inFileName+'/ZZTree/candTree')
    tree.Add(inFileName+'/ZZTree/candTree_failed')

    iEntry = 0
    while tree.GetEntry(iEntry):

        iEntry +=1

        theCand = Candidate(tree)
        cands.append(theCand)

        # if iEntry == 1500: break

    sortedCands = sorted(cands, key=attrgetter('run', 'lumi', 'event'))
    outFile = open(outFileName,"w")
    line = ''
    for aCand in sortedCands:
        line += aCand.printOut()
        line += "\n"

    outFile.write(line)
    outFile.close()


parseOptions()
loop()
