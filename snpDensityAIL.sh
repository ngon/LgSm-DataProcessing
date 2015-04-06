#! /bin/bash

#### scheduler directives ####

#PBS -N snpDensityAIL
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=120gb
#PBS -j oe


#### execute job #####

module load R/3.1.0
cd /group/palmer-lab/AIL/GBS/genoSummaries/
R CMD BATCH snpdensity.R