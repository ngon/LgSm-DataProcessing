#! /bin/bash

PROJ_HOME='/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing'
NSAMPS=`head -1 $PROJ_HOME/genotypes/ail.chr1.filtered.dosage | wc -w`
let NSAMPS=NSAMPS-3
echo "Analyzing with $NSAMPS samples"

cd $PROJ_HOME/kinship

## Make dummy phenotype file
if [ -e dummy.pheno ]; then 
    rm dummy.pheno
fi
for samp in `seq 1 $NSAMPS`; do 
    echo "1" >> dummy.pheno
done 

flags="-b y -cwd -j y -l h_vmem=6g"
if [ ! -e /home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing/kinship/allChrs.cXX.txt ]; then
    echo "
#!/bin/bash

#$ -j y 
#$ -cwd 
#$ -l h_vmem=10g 
#$ -N kinship 

/home/shyamg/bin/gemma -g $PROJ_HOME/genotypes/ail.allchr.filtered.dosage.gz -gk 1 -p dummy.pheno -maf 0 -o allChrs
mv output/allChrs.cXX.txt .
"| qsub;
    flags=$flags" -hold_jid kinship"
fi

cd $PROJ_HOME
traits=(cpp8.t sens act1.t act2.t act3.t act4.t act5.t act8.t avg.ppi startle glucose wild)
for trait in ${traits[@]}; do 
    col=`grep -w $trait $PROJ_HOME/pheno.names.txt | cut -f1`
    qsub $flags /home/shyamg/bin/gemma -g $PROJ_HOME/genotypes/dummy.herit.dosage -p $PROJ_HOME/phenos.allgeno.txt -a $PROJ_HOME/snpinfo/dummy.herit.snpinfo -k $PROJ_HOME/kinship/allChrs.cXX.txt -c covariates/$trait.covs -lmm 2 -o $trait.heritability -n $col
done
