AddOrReplaceReadGroups 
I= input bam 
O= output bam 
ID= 46352 # read group id
SM= 46352 # sample name
PL= illumina # read group platform
LB= merged # read group lib <- shyam has it listed as the bam directory
PU= unit1 # run barcode
CREATE_INDEX= true

# extract seq identifiers for control seqs and regular files to see how many
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
cit -f4 ailSeqIds | sort -u # unique flow cell lanes








