#! /bin/bash

#### scheduler directives ####

#PBS -N plotManhattanEmp
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -j oe


#### execute job #####

module load R/3.1.0
cd /group/palmer-lab/AIL/qtlmapping/output
R CMD BATCH /group/palmer-lab/AIL/LgSm-DataProcessing/scripts/plotManhattanGemma.R