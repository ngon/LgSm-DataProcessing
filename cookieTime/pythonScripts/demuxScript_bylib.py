#! /usr/bin/env python
##########################################################
# Code to demux the GBS raw sequence files using the li- #
# rary number. The input files are read from other files #
##########################################################

# You could import using the grammar 'from sys import argv' so you could type
# argv instead of sys.argv with each use... but know that this decreases
# readability and could results in name clashes. 

import os           # misc OS interfaces
import re           # regular expression operations
import argparse     # parser for cmd line options, args and sub-cmds
import gzip         # support for gzip files
import numpy as np  # NUMerical PYthon for scientific computing
                    # np = alias for numpy. purpose: avoid namespace conflicts.
import sys          # system-specific params and fxns
import datetime     # basic date and time types


#### TOP LEVEL CODE  ####
# __    __: double underscores designate special method names used by python
# __main__: this module represents the otherwise anonymous scope in which the
#           interpreter's main program executes. i.e. cmds read from standard
#           input, a script file, or an interactive prompt. the section below 
#           is teling the interpreter what envi it should run the script from.
# __name__: return the object represented as a string

# When the interpreter reads a source file, it executes all of the code within.
# Before execution, it will define a few variables. Here the interpreter is 
# running the source file as the main program . It sets __name__ to have the 
# value '__main__'. If this file is being imported from another module, then
# __name__ will be set to the module's name. 

# After setting up the special variables, python executes the import statement
# and loads libraries. Then it evaluates the def blocks, creating (1) function
# objects and (2) the variables that represent the functions. It will then read
# the if statement and see that __name__ does equal '__main__', so it will 
# execute the block shown there.

# One reason to do this is that sometimes you'll write a module (a .py file)
# where it can be executed directly. Alternatively, it could be imported and 
# used in another module. By doing the __main__ check, you can have that code
# execute only when you want to run the module DIRECTLY, and prevent it from
# being executed when someone just wants to import your module and call your
# functions themselves. If the code is being imported, the various function
# and class defs will be imported, but the main() code won't run. 

if (__name__=='__main__'):  # 'is __name__ being run standalone by the user?'

    # get ArgumentParser from the argparse module. ArgumentParser will hold all
    # the info needed to parse the cmd line into Python data types. 
    
    # argparse makes it easy to write user-friendly cmdline interfaces. The 
    # script defines what args it needs and argparse figures out how to parse 
    # those args out of sys.argv. argparse also generates help and useage msgs 
    # automatically and issues errors when users give the program invalid args. 
    
    parser = argparse.ArgumentParser(description='Demultiplexing script for GBS')

    # add_argument defines how a single cmdline arg should be passed. e.g.:
    
    # a list of option strings (flags or names): -l and --library. optional
    # args are identified by the - prefix and the remaining args are assumed to
    # be positional. this isn't relevant here since all the args use '-', but if
    # there was >parser.add_argument('foo') defined before the 1st arg then it'd
    # assume that 'foo' is required to be first. then it'd spit an error if you 
    # typed > python -l 'libraryName' ... in that case you'd have to type
    # > python 'something' -l 'libraryName' ...
    
    # metavar: the arg name to be displayed in usage messages (the default is
      # whatever is specified in 'dest'.
      
    # type: type (int? str?) to which the cmdline arg should be converted
    
    # help: a brief description of what each arg does
    parser.add_argument('-l', '--library', metavar='Library', type=str, 
                        dest='lib', help='Library to demultiplex', 
                        required=True)
                        
    # dest: the name of the attribute to be added to the object returned by
      # parse_args()
    parser.add_argument('-1', '--infirst', metavar='LeftFile', type=str, 
                        dest='left', 
                        help='File of input filenames (left reads if paired)',
                        required=True)

    # default: the value produced if the arg is absent from the cmdline
    parser.add_argument('-2', '--insecond', metavar='RightFile', type=str, 
                        dest='right', help='File of right read filenames', 
                        required=False, default='')
    
    parser.add_argument('-s', '--adapters', metavar='AdapterFile', type=str, 
                        dest='sidfname', help='Adapter file (tab separated)', 
                        required=True)
                        
    # action: type of action taken when this arg is encountered at the cmdline
    # there are a limited number of action types:
        # 'store': stores the argument's value (this is the default)
        # 'store_const': stores the value specified by the const keyword arg 
          # (default = None. this type is often used with optional args that 
          # specify some sort of flag. 
        # 'store_true' and 'store_false': special cases of 'store_const'. They
          # create default values of True and False, respectively.
        # 'append': stores a list and appends each arg value to the list. This 
          # is useful for allowing an option to be specified multiple times. 
        # 'append_const': stores a list and appends the value specified by const 
          # (useful when multiple args need to store constants to the same list)
        # 'count': counts the number of times a keyword arg occurs (useful for
          # increasing verbosity levels).
        # 'help': prints a complete help message for all options in the current
          # parser and then exits. 
        # 'version': this expects a version= keyword arg in the add_argument 
          # call. it prints version info and exits when invoked.   
        # CUSTOM: you can definte and add custom actions using Action. 
    
    parser.add_argument('-z', '--zipped', dest='zipped', 
                        help='Input file zipped', action='store_true')
    
    parser.add_argument('-p', '--pairedEnd', dest='paired', 
                        help='Paired end reads', action='store_true')
    
    parser.add_argument('-q', '--qseqs', dest='qseq', 
                        help='Qseq files not sequence files', 
                        action='store_true')
    
    parser.add_argument('-u', '--unzipoutput', dest='unzipout', 
                        help='Unzipped output', action='store_true')
    
    parser.add_argument('-e', '--enzymeSuffix', dest='enzyme', 
                        help='Overhang of Enzyme', required=False, 
                        default='TGCAG')

    # other parameters for add_argument() include:
    # nargs: normally a single cmdline arg is associated with a single arg. 
      # nargs associates a different num of cmdline args with a single action.
      # options for nargs include integers, ?, *, + (see pydoc for more info) 
    # const: constant value required by some action and nargs selections
    # choices: a container of allowable values for the argument

    # parse_args converts arg strings to objects and assigns them as attributes
    # of the namespace. it returns the populated namespace.                         
    args = parser.parse_args()
 
#### END TOP LEVEL CODE ####


# DOCSTRINGS
# triple-quoted strings in the functions below
# >>> help() fetches the docStrings attribute of the function in ()
# Convention: a multi-line string where the first line starts with a capital
# letter and ends with a dot. The 2nd line is blank. Detailed information starts
# from the third line. 



# Hamming distance between 2 strings of equal length is the num of positions
# at which the corresponding symbols are different. E.g. it measure the min
# num of substitutions/errors required to change one string into the other.
def hamming (s1, s2):
    """
    Computes the hamming distance between
    2 strings of bases.
    """
    if len(s1) != len(s2): return 50    # 50= #bases in common adapter seq?
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

# sidfname is the adapter file name (-s arg above)
sidfile = open(args.sidfname)
#enzymeSuffix = 'TGCAG' #works only for PstI
# enzyme suffix is given after -e; default is 'TGCAG'
enzymeSuffix = args.enzyme
line = sidfile.readline() # gets rid of header
sids = {}   # sid stands for sample id?
for line in sidfile:
    line = line.strip() # line is a string with whitespace removed
    toks = line.split() # alias for splitting a string into subunits
                        # ask shyam wtf is up with 'toks'
    adap = toks[1].upper()+enzymeSuffix # make subunit[1] (the idx) an uppercase
                                        # string & cat(TGCAG) to the end
    if (args.lib == toks[2]):   # the toks[2] must be the library name ('AIL50')
        sids[adap] = toks[0]    # toks[0] is the mouse id
sidfile.close()

# the function below handles paired end reads and zipped files
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
firsts = open(args.left)  # this will be the only files for SE reads
if args.right != "":      # e.g. 'if paired end' there will be righthand reads
    seconds = open(args.right)
if (args.qseq):           # handles qseq files
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
