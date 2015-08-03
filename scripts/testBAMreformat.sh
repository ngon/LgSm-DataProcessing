AddOrReplaceReadGroups 
I= input bam 
O= output bam 
ID= 46352 # read group id
SM= 46352 # sample name
PL= illumina # read group platform
LB= merged # read group lib <- shyam has it listed as the bam directory
PU= unit1 # run barcode
CREATE_INDEX= true

# 1. extract seq identifiers for control seqs and regular files to see how many
# unique machine IDs there are. also get run and lane info for each sample.
cd /group/palmer-lab/AIL/GBS

for FILE in /group/palmer-lab/AIL/GBS/bams/controls/*rg.requaled.bam; do
    samtools view $FILE | head -1 | cut -f1 >> controlSeqIds.txt
    basename $FILE >> controlFileNames.txt
done

for FILE in /group/palmer-lab/AIL/GBS/bams/*rg.requaled.realign.bam; do
    samtools view $FILE | head -1 | cut -f1 >> ailSeqIds.txt
    basename $FILE >> ailFileNames.txt
done

for FILE in /group/palmer-lab/AIL/GBS/bams/reps/*rg.requaled.realign.bam; do
    samtools view $FILE | head -1 | cut -f1 >> repsSeqIds.txt
    basename $FILE >> repsFileNames.txt
done

# combine the 3 files w/ cat, then merge them into one file with 2 cols
cut -f1 ailSeqIds.txt | paste ailFileNames.txt - > ailSeqIds

# in vim, replace all colons with tabs and join Casava FC IDs & lane# with a
# dash - so that the 2 elements comprise a single field. 
cut -f3 ailSeqIds | sort -n -u # unique run IDs
cut -f4 ailSeqIds | sort -u # unique flow cell lanes


### For all HW1 machines:
# 1. remove the 3rd (C6A... flowcell ID) field
# 2. replace ' '1:N:O  with :0:1 (1 means SE reads; on HiSeq, control field is 
#    always 0; Y=read is filtered (did not pass) and N=otherwise.

 HWI-D00422:152:C6A0UANXX:3:2216:13634:67132 1:N:0 # this is casava format
 7LYMFP1:386:2:2108:15829:46353:0:1 # this Illumina format
 7LYMFP1:386:2:2108:15829:46353 1:N:0:index # this is the format I want

# 1. sed with eregex: delete FCID sequence
sed -r 's/C6A[0-9][0-9A-Z]ANXX://g' ailSeqIds.txt > rmfcid.test

# 2. change :0:1$ in lines !^HW to 1:N:0 so we can just append barcodes to 
#    the end of each line.
sed -r 's/:0:1$/ 1:N:0/g' rmfcid.apnd.test > rmfcid.apnd.test1 


# works 
egrep -c '^HWI-' ailSeqIds.txt
egrep -c 'HWI-[0-9A-Z]' ailSeqIds.txt
egrep -c 'HWI-[0-9A-Z].' ailSeqIds.txt
egrep -c 'HWI-[0-9A-Z].[0-9A-X].' ailSeqIds.txt
# doesn't work
egrep -c '^HWI[-0-9A-Z]:[0-9A-X]:' ailSeqIds.txt 
egrep -c 'HWI-[0-9A-Z]:' ailSeqIds.txt








 



















