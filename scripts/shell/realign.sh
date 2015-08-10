#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N realign.ail.2

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=01:30:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 1 core on a single node.
#PBS -l nodes=1:ppn=1

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
module load gatk/3.3-0

INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/bams/realign.ail.2.list | tail -1`
BASE=`basename $INFILE .bam`
DIR=`dirname $INFILE`
echo "Running realignment for: $INFILE" 
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar -T IndelRealigner -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -I $INFILE -targetIntervals /group/palmer-lab/reference_genomes/mouse/LG_SM.mm10.realign.intervals -o $DIR/$BASE.ra.bam
echo "Done running realignment."

