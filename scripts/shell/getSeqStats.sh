#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N vcf.read.stats

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=02:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 4 cores on a single node.
#PBS -l nodes=1:ppn=1

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=2gb

### Merge the output and error streams
#PBS -j oe

#################
# Job Execution #
#################

module load vcftools/0.1.11
cd /group/palmer-lab/AIL/GBS/vars
echo "Running command:" `tail -n +$PBS_ARRAYID /group/palmer-lab/AIL/GBS/vars/getSeqStats.vcf.cmds | head -1`
tail -n +$PBS_ARRAYID /group/palmer-lab/AIL/GBS/vars/getSeqStats.vcf.cmds | head -1 | sh
echo "Done running command."

