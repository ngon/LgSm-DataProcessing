############# Code to compute the Kinship matrix for each chromosome separately ##################
PROJ_HOME='/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing'
NSAMPS=`head -1 $PROJ_HOME/genotypes/ail.chr1.dosage | wc -w`
let NSAMPS=NSAMPS-3
echo "Analyzing with $NSAMPS samples"

## Make dummy phenotype file
#rm dummy.pheno
#for samp in `seq 1 $NSAMPS`; do 
#    echo "1" >> dummy.pheno
#done 

## Compute the kinship matrices
for chr in {1..1}; do 
echo "#! /bin/bash
#$ -cwd
#$ -l h_vmem=5g
#$ -j y

echo Starting chr $chr
cd $PROJ_HOME/genotypes
cat `ls $PROJ_HOME/genotypes | sort -k2.4,2n -t. | grep -v $chr\\\. | tr -s \"\n\" \" \"` > $PROJ_HOME/kinship/notChr$chr.dosage
cd $PROJ_HOME/kinship 
/home/shyamg/bin/gemma -g notChr$chr.dosage -p dummy.pheno -gk 1 -o notChr$chr
rm -f notChr$chr.dosage
echo Done with chr $chr." | qsub
done
