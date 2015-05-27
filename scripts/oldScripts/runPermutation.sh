#! /bin/bash

#### Code to run permutations for the trait of interest. 

trait=$1
numperms=$2

PROJ_HOME="/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing"

if [ ! -e $PROJ_HOME/kinship/identity.kinship ]; then
    nsamps=`head -1 $PROJ_HOME/genotypes/ail.chr1.filtered.dosage | wc -w`
    let nsamps=nsamps-3
    echo "Number of samples: $nsamps"
    ssh spudling25 "Rscript -e 'write.table(diag($nsamps), file=\"$PROJ_HOME/kinship/identity.kinship\", sep=\" \", quote=F, row.names=F, col.names=F)'"
fi

colnum=`grep -w $trait $PROJ_HOME/pheno.names.txt | cut -f1`

#### Create permutations directory structure
if [ ! -e $PROJ_HOME/permutations ]; then
    mkdir $PROJ_HOME/permutations
    mkdir $PROJ_HOME/permutations/covariates
    mkdir $PROJ_HOME/permutations/logs
fi

# use ssh to run on a spudling so that we do not have to ensure that it runs no matter what machine we run it on 
ssh spudling25 "Rscript -e 'a <- read.table(\"$PROJ_HOME/phenos.allgeno.txt\"); covs <- read.table(\"$PROJ_HOME/covariates/$trait.covs\"); chosen.trait <- a[,$colnum]; out.perms <- sapply(1:$numperms, function(x) { nona <- which(!is.na(chosen.trait)); ord<-sample(nona); temp.covs <- covs; temp.covs[nona, ] <- covs[ord,]; write.table(temp.covs, file=paste0(\"$PROJ_HOME/permutations/covariates/$trait.perm\", x, \".covs\"), sep=\" \", quote=F, row.names=F, col.names=F); temp.trait <- chosen.trait; temp.trait[nona] <- chosen.trait[ord]; temp.trait }); write.table(out.perms, file=\"$PROJ_HOME/permutations/$trait.permuted.txt\", sep=\" \", quote=F, row.names=F, col.names=F);'"

echo "#! /bin/bash

#$ -N $trait.perms
#$ -cwd
#$ -j y 
#$ -l h_vmem=3g
#$ -t 1:$numperms
#$ -o $PROJ_HOME/permutations/logs/$trait.perm\$TASK_ID.log
#$ -tc 50

cd $PROJ_HOME/permutations
/home/shyamg/bin/gemma -g $PROJ_HOME/genotypes/ail.allchr.filtered.dosage.gz -k $PROJ_HOME/kinship/identity.kinship -a $PROJ_HOME/snpinfo/ail.allchrs.snpinfo -p $PROJ_HOME/permutations/$trait.permuted.txt -n \$SGE_TASK_ID -c $PROJ_HOME/permutations/covariates/$trait.perm\${SGE_TASK_ID}.covs -o $trait.perm\${SGE_TASK_ID} -lmm 2 -maf 0.05" | qsub;

