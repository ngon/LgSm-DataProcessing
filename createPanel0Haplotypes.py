#####################################
# Code to create impute2 panel0 haps#
# and legend files from VCF. Splits #
# by chromosome. THE INPUT GENOTYPES#
# WHICH ARE HETEROZYGOTE ARE SET TO #
# REFERENCE HOMOZYGOTE.             #
#####################################

import sys

if (len(sys.argv) != 3):
    print "Usage: createPanel0Haplotypes.py InputVCF OutputDirectory"
    sys.exit(1)

infile=open(sys.argv[1])

hapfiles = {}
legendfiles = {}
for chrom in range(1,20):
    outfile=open(sys.argv[2]+"/"+"chr"+str(chrom)+".hap", "w")
    hapfiles["chr"+str(chrom)] = outfile
    outfile=open(sys.argv[2]+"/"+"chr"+str(chrom)+".legend", "w")
    outfile.write("rs#\tpos\tref\talt\n")
    legendfiles["chr"+str(chrom)] = outfile

outfile=open(sys.argv[2]+"/"+"chrX.hap", "w")
hapfiles["chrX"] = outfile
outfile=open(sys.argv[2]+"/"+"chrX.legend", "w")
outfile.write("rs#\tpos\tref\talt\n")
legendfiles["chrX"] = outfile

cnt = 0
for line in infile:
    if line[0] == "#": continue
    toks = line.strip().split()
    if toks[0] not in legendfiles: continue
    if len(toks[3]) != 1 or len(toks[4]) != 1: continue
    legendfiles[toks[0]].write("ail."+toks[0]+"."+toks[1]+"\t"+toks[1]+"\t"+toks[3]+"\t"+toks[4]+"\n")
    curhapfile = hapfiles[toks[0]]
    for geno in toks[9:]:
        if geno == "1/1":
            curhapfile.write("1 ")
        else:
            curhapfile.write("0 ")
    curhapfile.write("\n")
    cnt += 1
    if cnt%100000 == 0:
        print cnt, "markers processed."

for fl in hapfiles.values():
    fl.close()

for fl in legendfiles.values():
    fl.close()
            
