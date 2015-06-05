#! /bin/bash

#### scheduler directives ####

#PBS -N append.cM.emp
#PBS -S /bin/bash
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -j oe


#### execute job #####

module load R/3.1.0
cd /group/palmer-lab/AIL/qtlmapping
R CMD BATCH /group/palmer-lab/AIL/LgSm-DataProcessing/scripts/oldScripts/append.cM.R
