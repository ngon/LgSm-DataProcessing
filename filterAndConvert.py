##################################
# Code to filter the variants    #
# based on info score from       #
# impute2. Default keeps vars    #
# with info > 0.3. Also converts #
# them to gemma format with 2 out#
# files - dosage and snpinfo     #
##################################


import sys
import os
import re
import numpy as np
import argparse

if (__name__=='__main__'):
    parser = argparse.ArgumentParser(description='Convert impute2 output files to dosage (gemma) format after filtering')
    parser.add_argument('-f', '--input', metavar='InputFile', type=str, dest='infile', help='Input impute2 gprobs file', required=True)
    parser.add_argument('-o', '--output', metavar='OutputPrefix', type=str, dest='out', help='Output file name', required=True)
    parser.add_argument('-i', '--info', metavar='InfoThreshold', type=float, dest='info', help='Information threshold', required=False, default=0.3)
    args = parser.parse_args()

infile = open(args.infile)
## assumes info filename is infile+"_info" - should
## work if no changes to impute2 file naming scheme
infofile = open(args.infile+"_info")
infoline = infofile.readline()
outfile = open(args.out+".dosage", 'w')
snpfile = open(args.out+".snpinfo", 'w')

cnt = 0
for (inline, infoline) in zip(infile, infofile):
    infotoks = infoline.strip().split()
    if (float(infotoks[6]) < args.info): continue
    toks = inline.strip().split()
    gps = np.array([float(x) for x in toks[5:]])
    gps.resize(len(gps)/3, 3)
    dosages = [str(x) for x in np.round(gps[:,2]*2+gps[:,1], 3)]
    outfile.write(toks[1]+" "+toks[3]+" "+toks[4]+" "+" ".join(dosages)+"\n")
    snpfile.write(infotoks[1]+"\t"+infotoks[2]+"\t"+infotoks[1].split('.')[1]+"\n")
    cnt += 1
    if cnt %10000 == 0:
        outfile.flush()
        snpfile.flush()
        print cnt, "markers done."

infile.close()
infofile.close()
outfile.close()
snpfile.close()
