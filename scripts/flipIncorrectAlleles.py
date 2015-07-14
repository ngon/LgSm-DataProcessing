#! /usr/bin/env python

import sys
import os

## Function to change strand
def changeStrand(base):
    if base == 'A':
        return 'T'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    return base

## Check that the correct number of arguments is given
if (len(sys.argv) != 4):
    print "Usage: flipIncorrectAlleles.py Reference.fasta Input.vcf Output.vcf"
    sys.exit(1)
    
## Reading in the command line arguments
reffile = sys.argv[1]
vcffile = sys.argv[2]
outfile = sys.argv[3]

## SCRIPT PARAMS
keepflipped = False

## Printing the command line arguments
print reffile, vcffile, outfile

## Initialize the variables for reading the sequence dictionary.
## Not storing the sequence in a pickle file - no significant 
## savings in terms of time or memory
seqDict = {}
curChr = ""
chrString = ""
## Open the sequence file
refin = open(reffile)
## For loop for lines in the sequence file
for line in refin:
    line = line.strip()
    ## Check the line to see if it starts with a >
    ## then it is a line with the chromosome name
    if line.startswith('>'):
        ## If this is not the first time you are
        ## seeing a chromosome name, then you need
        ## to store the previous chromsome's sequence
        if (curChr != ""):
            seqDict[curChr] = chrString
            print "\nDone with chr", curChr
        ## Extract chromosome name and initialize 
        ## the sequence dictionary's key for this chr.
        ## Reset the sequence string for this chr since
        ## we are starting a new chromosome.
        curChr = line[1:]
        seqDict[curChr] = ""
        chrString = ""
        print "Starting chr", curChr
    else:
        ## Sequence line. Just store the line in the
        ## string for the current chromosome.
        chrString += line
        ## Track how many bases have been read. For every
        ## 10M bases output a dot as progress indicator.
        if (len(chrString)%1e7 == 0):
            print "\b.",
            sys.stdout.flush()
## End for loop

## Deal with the last chromosome, it has not been
## written into the seqDict dictionary yet. 
if (chrString != ""):
    seqDict[curChr] = chrString
    print "Done with chr", curChr
refin.close()

## Read the vcf file for the LGSM in mm10 - LGSM.mm10.vcf
## Open the input file and output file.
vcfin = open(vcffile)
out   = open(outfile, "w")
cnt = 0
## For loop for lines in the vcf file.
for line in vcfin:
    ## If the line starts with comment character, just output it 
    ## to the output file.
    if line.startswith('#'):
        out.write(line)
        continue
    cnt += 1
    ## Separate the lines by white space to get tokens
    tokens = line.strip().split()
    ## Get the reference allele from the fasta sequence at this
    ## position - convert it to uppercase to make sure that it
    ## matches the vcf file.
    refallele = seqDict[tokens[0]][int(tokens[1])-1].upper()
    ## Check if the reference allele in the fasta sequence
    ## matches refernce allele indicated by the vcf file.
    if refallele == tokens[3]: # the vcf and fasta ref alleles match
        out.write(line)
    elif not keepflipped:
        continue
    elif refallele == tokens[4]: # the vcf alt and fasta ref match
        temp = tokens[3]
        tokens[3] = tokens[4]
        tokens[4] = temp
        out.write('\t'.join(tokens))
        out.write('\n')
        print "flipped alleles:", tokens[0], tokens[1]
    elif refallele == changeStrand(tokens[3]): # the vcf ref allele in on the opposite strand and matches fasta ref
        tokens[3] = refallele
        tokens[4] = changeStrand(tokens[4])
        out.write('\t'.join(tokens))
        out.write('\n')
        print "flipped strand:", tokens[0], tokens[1]
    else: # fasta reference allele does not match either allele in the vcf file.
        print "Neither ref nor alt allele match the mm10 reference:", tokens[3], tokens[4], refallele
    ## Print progress bar.
    if (cnt%100000 == 0):
        print cnt, "sites done."

## Close the input vcf and output files.
vcfin.close()
out.close()
