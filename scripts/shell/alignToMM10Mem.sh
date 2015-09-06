#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N alignerGBS

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=24:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 4 cores on a single node.
#PBS -l nodes=1:ppn=4

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=10gb

### Inform the scheduler that it is an array job and tell it to limit the
### number of simultaneous jobs being run.

### Merge the output and error streams
#PBS -j oe

#################
# Job Execution #
#################

module load bwa

# Script to align using bwa mem
# INFILE is a text file containing a list of files to process; one file per line.
INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/unmapped.txt | tail -1`
DIR='/group/palmer-lab/AIL/GBS/fastqs/'
OUTDIR='/group/palmer-lab/AIL/GBS/bams/'
BASENAME=`basename $INFILE .fq.gz`

RG="@RG\tID:$BASENAME\tSM:$BASENAME\tPL:illumina\tLB:$OUTDIR\tPU:unit1"
REFGENOME='/group/palmer-lab/reference_genomes/mouse/mm10.fa.gz'
SAMSORT='/home/sgopalakrishnan/tools/picard-tools-1.126/picard.jar SortSam'
JAVA='/usr/bin/java'
STDIN='/dev/stdin'

echo "Running alignment for $INFILE"
if [ ! -s $OUTDIR/${BASENAME}.bam ]; then
  bwa mem -M -R "$RG" $REFGENOME $INFILE | $JAVA -Xmx4g -jar $SAMSORT INPUT=$STDIN OUTPUT=$OUTDIR/$BASENAME.rg.bam SORT_ORDER=coordinate CREATE_INDEX=true
fi

exitstatus=( ${PIPESTATUS[@]} )
if [ ${exitstatus[0]}  != 0 ]; then
  echo "ERROR: Failed bwa mem with exit code "${exitstatus[0]}
  exit ${exitstatus[0]}
elif [ ${exitstatus[1]} != 0 ]; then
  echo "ERROR: Could not sort and make the bam file from sam file "${exitstatus[1]}
  exit ${exitstatus[1]}
else
  echo "Successfully mapped the reads and then made, sorted and indexed the bam file"
fi
