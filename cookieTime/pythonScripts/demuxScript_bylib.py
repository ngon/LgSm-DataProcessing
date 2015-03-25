#! /usr/bin/env python
##########################################################
# Code to demux the GBS raw sequence files using the li- #
# rary number. The input files are read from other files #
##########################################################

import os
import re
import argparse
import gzip
import numpy as np
import sys
import datetime

if (__name__=='__main__'):
    parser = argparse.ArgumentParser(description='Demultiplexing script for GBS')
    parser.add_argument('-l', '--library', metavar='Library', type=str, dest='lib', help='Library to demultiplex', required=True)
    parser.add_argument('-1', '--infirst', metavar='LeftFile', type=str, dest='left', help='File of input filenames (left reads if paired)', required=True)
    parser.add_argument('-2', '--insecond', metavar='RightFile', type=str, dest='right', help='File of right read filenames', required=False, default='')
    parser.add_argument('-s', '--adapters', metavar='AdapterFile', type=str, dest='sidfname', help='Adapter file (tab separated)', required=True)
    parser.add_argument('-z', '--zipped', dest='zipped', help='Input file zipped', action='store_true')
    parser.add_argument('-p', '--pairedEnd', dest='paired', help='Paired end reads', action='store_true')
    parser.add_argument('-q', '--qseqs', dest='qseq', help='Qseq files not sequence files', action='store_true')
    parser.add_argument('-u', '--unzipoutput', dest='unzipout', help='Unzipped output', action='store_true')
    parser.add_argument('-e', '--enzymeSuffix', dest='enzyme', help='Overhang of Enzyme', required=False, default='TGCAG')
    args = parser.parse_args()

def hamming (s1, s2):
    """
    Computes the hamming distance between
    2 strings of bases.
    """
    if len(s1) != len(s2): return 50
    dist = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]: dist += 1
    return dist

def nearestNeigh(readstr, adapters):
    """
    This function calculates the hamming distance 
    between the read and all the adapters. 
    Returns the nearest adapter if the minimium 
    hamming distance is 0 or 1. If not, returns 
    NULL.
    """
    minDist = 50
    nearest = None
    numnearest = 0
    for aid in adapters:
        cd = hamming(readstr[0:len(aid)], aid)
        if cd < 2:
            if cd < minDist:
                minDist = cd
                nearest = aid
                numnearest += 1
            elif cd == minDist:
                numnearest += 1
                nearest = None
    if numnearest > 1: return None
    return nearest

sidfile = open(args.sidfname)
#enzymeSuffix = 'TGCAG' #works only for PstI
enzymeSuffix = args.enzyme
line = sidfile.readline() # gets rid of header
sids = {}
for line in sidfile:
    line = line.strip()
    toks = line.split()
    adap = toks[1].upper()+enzymeSuffix
    if (args.lib == toks[2]):
        sids[adap] = toks[0]
sidfile.close()

fhands_1 = {}
if (args.paired):
    fhands_2 = {}
if (args.unzipout):
    for name in np.unique(sids.values()):
        ftemp = open(name+'.fq', 'w')
        fhands_1[name] = ftemp
        if (args.paired):
            ftemp = open(name+'.fq', 'w')
            fhands_2[name] = ftemp
else:
    for name in np.unique(sids.values()):
        ftemp = gzip.open(name+'.fq.gz', 'w')
        fhands_1[name] = ftemp
        if (args.paired):
            ftemp = gzip.open(name+'.fq.gz', 'w')
            fhands_2[name] = ftemp

cntReads = 0
success = 0
firsts = open(args.left)
if args.right != "":
    seconds = open(args.right)
if (args.qseq):
    fs=[x.strip() for x in firsts.readlines()]
    for fname1 in fs:
        if(args.paired):
            fname3 = seconds.readline().strip()
            if (args.zipped):
                f3=gzip.open(fname3)
            else:
                f3=open(fname3)
        if (args.zipped):
            f1=gzip.open(fname1)
        else:
            f1=open(fname1)
        if (args.paired):
            for l1, l3 in zip(f1, f3):
                l1=l1.strip()
                t1=l1.split()
                l3=l3.split()
                t3=l3.split()
                nearest = nearestNeigh(t1[8], sids)
                if(nearest != None):
                    fhands_1[sids[nearest]].write('@'+':'.join(t1[0:8])+'\n')
                    fhands_1[sids[nearest]].write(t1[8][len(nearest):]+'\n+\n'+t1[9][len(nearest):]+'\n')
                    fhands_2[sids[nearest]].write('@'+':'.join(t3[0:8])+'\n')
                    fhands_2[sids[nearest]].write(t3[8][len(nearest):]+'\n+\n'+t3[9][len(nearest):]+'\n')
        else:
            for l1 in f1:
                l1=l1.strip()
                t1=l1.split()
                nearest = nearestNeigh(t1[8], sids)
                if (nearest != None):
                    fhands_1[sids[nearest]].write('@'+':'.join(t1[0:8])+'\n')
                    fhands_1[sids[nearest]].write(t1[8][len(nearest):]+'\n+\n'+t1[9][len(nearest):]+'\n')
        f1.close()
        if(args.paired):
            f3.close()
else:
#######NOT A QSEQ -- implement part for fq files.
    sname = ''
    curcnt = 0
    if (args.paired):
        print 'Cannot deal with paired end data with fq input file.'
        sys.exit(1)
    files=os.listdir('.')
    fs=[x.strip() for x in firsts.readlines()]
    print 'Processing', len(fs), 'files.'
    for fname1 in fs:
        if (args.zipped):
            f1=gzip.open(fname1)
        else:
             f1=open(fname1)
        for l1 in f1:
            if (curcnt == 0):
                if (len(l1) == 0):
                    break
                l1=l1.strip()
                if l1[0] == '@':
                    if (sname == ''):
                        sname = l1[1:]
                    else:
                        print 'Error after reading read', sname, l1
                        print 'Read', cntReads, 'reads.'
                        sys.exit(2)
                curcnt += 1
            elif (curcnt == 1):
                seq = l1.strip()
                curcnt += 1
            elif (curcnt == 2):
                if l1[0] != '+':
                    print 'Error after reading read -', sname, '- -', l1, '-'
                    print 'Read', cntReads, 'reads.'
                    sys.exit(3)
                curcnt += 1
            elif (curcnt == 3):
                qual = l1.strip()
                nearest = nearestNeigh(seq, sids)
                if (nearest != None):
                    success += 1
                    fhands_1[sids[nearest]].write('@'+sname+'\n')
                    fhands_1[sids[nearest]].write(seq[len(nearest):]+'\n+\n'+qual[len(nearest):]+'\n')
                curcnt = 0
                sname = ''
                cntReads += 1
                if (cntReads % 100000 == 0):
                    print datetime.datetime.now().time(), ': Processed', cntReads, 'reads.'
                    sys.stdout.flush()
        f1.close()
	print "Done processing file."

for name in np.unique(sids.values()):
    fhands_1[name].close()
    if (args.paired):
        fhands_2[name].close()

firsts.close()
if args.right != "":
    seconds.close()
print 'Number of reads processed:', success, '/', cntReads
