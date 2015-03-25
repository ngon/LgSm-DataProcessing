import sys
import re
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert vcf file to impute2 genotype file')
    parser.add_argument('-v', '--vcf', metavar='vcfFile', type=str, dest='vcf', help='Input vcf file name', required=True)
    parser.add_argument('-c', '--chr', metavar='chrom', type=str, dest='chrom', help='Quality cutoff for SNPs', required=True)
    parser.add_argument('-o', '--out', metavar='outFile', type=str, dest='out', help='Output geno file name', required=False, default='test')
    parser.add_argument('-l', '--header', metavar='numHeaderLines', type=int, dest='header', help='Number of header lines', required=False, default=115)
    parser.add_argument('-q', '--qual', metavar='qualScore', type=float, dest='qual', help='Quality cutoff for SNPs', required=False, default=0)
    parser.add_argument('-g', '--glcol', metavar='PLColumn', type=int, dest='glcol', help='Likelihood column (0-indexed)', required=False, default=3)
    args = parser.parse_args()
    vcf = args.vcf
    out = args.out
    qual = args.qual
    chrom = args.chrom
    GLCOL = args.glcol
    if (chrom not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']):
        print 'Invalid chromosome.'
        sys.exit(1)
else:
    print 'Can be used only from command line.'
    sys.exit(1)

f = open(vcf, 'r')
for i in range(args.header):
    line = f.readline()

o = open(out, 'w')

line = f.readline()[0:-1]
samps = line.split()[9:]
nsamps = len(samps)

for line in f:
    line = line.strip()
    toks = line.split()
    if (toks[0] != 'chr'+chrom): continue
    if (toks[6] != '.' and toks[6] != 'PASS'): continue
    if (float(toks[5]) < qual): continue
    if (len(toks[3]) != 1 or len(toks[4]) != 1): continue
    sys.stdout.flush()
    o.write(toks[0]+'\tail.'+ toks[0]+'.'+toks[1]+'\t'+toks[1]+'\t'+toks[3]+'\t'+toks[4]+'\t')
    o.flush()
    for lstr in toks[9:]:
        geno = lstr.split(':')[0]
        if geno != './.' and geno != '.':
            gls = np.array([float(x) for x in (lstr.split(':')[GLCOL]).split(',')])
        else:
            gls = np.array([1.,1.,1.])
        gls = 10**(gls/-10.)
        gls = gls/np.sum(gls)
        o.write('\t'.join([str(round(x,3)) for x in gls]))
        o.write('\t')
    o.write('\n')
    o.flush()
f.close()
o.close()
