#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N requal.ctrl

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=24:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 4 cores on a single node.
### NOTE: Shyam originally set ppn=4. I changed it to 1 at the recommendation of
### HPC personnel. NMG
#PBS -l nodes=1:ppn=1

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

# Run the commands
CMD=`head -$PBS_ARRAYID /group/palmer-lab/AIL/code/requal.ctrl.conversion.cmds | tail -1`
echo "Running command $CMD"
$CMD
echo "Done running command."

