#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N IndexBam

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=24:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 4 cores on a single node.
#PBS -l nodes=1:ppn=4

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=1gb

### Merge the output and error streams
#PBS -j oe

#################
# Job Execution #
#################

# Load the approprite applications
module load samtools

INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/ail.toindex.list | tail -n 1` 
BASE=`basename $INFILE .bam`
cd /group/palmer-lab/AIL/GBS/bams
samtools index $BASE.bam $BASE.bai
echo "Done generating index for $BASE"

