#! /bin/bash
#$ -cwd
#$ -j y

/home/shyamg/bin/gemma -g genotypes/ail.chr1.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr1.cXX.txt -a snpinfo/ail.chr1.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr1 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr2.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr2.cXX.txt -a snpinfo/ail.chr2.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr2 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr3.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr3.cXX.txt -a snpinfo/ail.chr3.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr3 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr4.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr4.cXX.txt -a snpinfo/ail.chr4.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr4 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr5.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr5.cXX.txt -a snpinfo/ail.chr5.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr5 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr6.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr6.cXX.txt -a snpinfo/ail.chr6.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr6 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr7.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr7.cXX.txt -a snpinfo/ail.chr7.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr7 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr8.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr8.cXX.txt -a snpinfo/ail.chr8.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr8 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr9.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr9.cXX.txt -a snpinfo/ail.chr9.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr9 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr10.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr10.cXX.txt -a snpinfo/ail.chr10.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr10 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr11.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr11.cXX.txt -a snpinfo/ail.chr11.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr11 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr12.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr12.cXX.txt -a snpinfo/ail.chr12.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr12 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr13.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr13.cXX.txt -a snpinfo/ail.chr13.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr13 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr14.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr14.cXX.txt -a snpinfo/ail.chr14.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr14 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr15.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr15.cXX.txt -a snpinfo/ail.chr15.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr15 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr16.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr16.cXX.txt -a snpinfo/ail.chr16.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr16 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr17.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr17.cXX.txt -a snpinfo/ail.chr17.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr17 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr18.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr18.cXX.txt -a snpinfo/ail.chr18.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr18 -n 3
/home/shyamg/bin/gemma -g genotypes/ail.chr19.filtered.dosage -p phenos.allgeno.txt -k kinship/notChr19.cXX.txt -a snpinfo/ail.chr19.snpinfo -c covariates/sens.covs -lmm 2 -maf 0.05 -o sens.chr19 -n 3
