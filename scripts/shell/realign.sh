#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N Realign

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=24:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 2 cores on a single node.
#PBS -l nodes=1:ppn=2

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=3gb

### Merge the output and error streams
#PBS -j oe

#################
# Job Execution #
#################

# Load the approprite applications
module load java
module load gatk/2.7-4

INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/bam.ail.requaled.list | tail -1`
BASE=`basename $INFILE .bam`
DIR=`dirname $INFILE`
echo "Running realignment for: $INFILE" 
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/2.7-4/GenomeAnalysisTK.jar -T IndelRealigner -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -I $INFILE -targetIntervals /group/palmer-lab/reference_genomes/mouse/WT.realign.intervals -o $DIR/$BASE.realign.bam
echo "Done running realignment."

