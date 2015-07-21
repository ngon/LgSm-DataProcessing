#! /bin/bash

########################
# Scheduler directives #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N impute.ail

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=48:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate 4 cores on a single node.
#PBS -l nodes=1:ppn=4

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=8gb

### Inform the scheduler that it is an array job and tell it to limit the
### number of simultaneous jobs being run.
##PBS -t 1-20

### Merge the output and error streams
#PBS -j oe

#################
# Job Execution #
#################

echo "Running command:" `tail -n +$PBS_ARRAYID /group/palmer-lab/AIL/code/ail.impute.cmds | head -1`
tail -n +$PBS_ARRAYID /group/palmer-lab/AIL/code/ail.impute.cmds | head -1 | sh
echo "Done running command."

