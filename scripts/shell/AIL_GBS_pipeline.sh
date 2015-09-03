#################################################
# Code to do the pipeline for the AIL animals.  #
# Various steps include alignment, requaling,   #
# realignment, variant calling, imputation etc. #
#################################################

## !!! SUPER IMPORTANT !!! ###
## WAIT FOR THE COMMANDS IN EACH SECTION TO FINISH RUNNING.
## you can monitor them using qstat and qstat -u $USER

#########################
#    Alignment steps    #
#########################
## Create the file with the list of fastqs that need to be aligned.
ls /group/palmer-lab/AIL/GBS/fastqs/*fq.gz > /group/palmer-lab/AIL/GBS/fastqs.ail.list
## MANUALLY EDIT THIS FILE to remove founders, F1s and other irrelevant fastqs

## submit the alignment script using the fastqs.ail.list file
## Submit it as a task array

## IMPORTANT: CRI doesn't support the syntax 1-n%n for array jobs, but the pps cluster
## probably does, which is why Shyam wrote it this way. It's a good idea to submit jobs
## in batches, but I think it has to be done manually on the CRI system. -- NMG
numfiles=`wc -l /group/palmer-lab/AIL/GBS/fastqs.ail.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numfiles}%200 /group/palmer-lab/AIL/code/alignToMM10Mem.sh


#########################
#  BAM Correction steps #
#########################
## These steps are needed to get all the 
## bam files to have a standard quality 
## encoding to ensure smooth processing

## Find out the bam quality encoding
## View the top 100000 reads of the bam file
## and check the quality string for certain
## characters to find quality coding.
## Refer to page on fastq in wikipedia for differen
## quality codings and their meanings
## http://en.wikipedia.org/wiki/FASTQ_format
## different encoding - we want to distinguish
## only Illumina 1.3 or 1.5 (phred+64) from
## sanger or illumina 1.9 (phred+33)
module load samtools
for bamfile in /group/palmer-lab/AIL/GBS/*bam; do          ### should probably be /GBS/bams/*bam ###
    n64=`samtools view $i | cut -f11 | head -100000 | grep -c [a-h]` 
    n32=`samtools view $i | cut -f11 | head -100000 | grep -c [0-9]` 
    ns64=`samtools view $i | cut -f11 | head -100000 | grep -c [\<\=\>?]`
    if [ $n32 -eq 0 ]; then 
	if [ $ns64 -ne 0 ]; then 
	    echo $bamfile Solexa - ignoring file because Solexa is obsolete
	else 
	    echo /group/palmer-lab/shared_code/convertQualScoreBam $bamfile >> /group/palmer-lab/AIL/code/requal.conversion.cmds 
	fi
    else 
	if [ $n64 -eq 0 ]; then 
	    echo mv $bamfile ${bamfile/.bam/.requaled.bam} >> /group/palmer-lab/AIL/code/requal.conversion.cmds
	else 
	    echo $bamfile Impossible $n32 $n64 $ns64 - so ignoring file.
	fi
    fi
done
## Run the conversion scripts
numfiles=`wc -l /group/palmer-lab/AIL/code/requal.conversion.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines}%200 /group/palmer-lab/AIL/code/requal.sh

#########################
## 13 Jul 2015 - Natalia
## BAM correction for control samples (LG, SM, F1)
## I didn't have execution permission on the original convertQualScoreBam file so I had to make a copy
## which I called convertQual_tmp in order to run it. Otherwise I would have had to download the C version
## of samtools, compile it, and then compile the convertQualScoreBam script again (it needs the bam.h lib
## from samtools/C, and possibly other packages...)
module load samtools
for bamfile in /group/palmer-lab/AIL/GBS/bams/controls/*bam; do 
    n64=`samtools view $i | cut -f11 | head -100000 | grep -c [a-h]` 
    n32=`samtools view $i | cut -f11 | head -100000 | grep -c [0-9]` 
    ns64=`samtools view $i | cut -f11 | head -100000 | grep -c [\<\=\>?]`
    if [ $n32 -eq 0 ]; then 
	if [ $ns64 -ne 0 ]; then 
	    echo $bamfile Solexa - ignoring file because Solexa is obsolete
	else 
	    echo /group/palmer-lab/shared_code/convertQual_tmp $bamfile >> /group/palmer-lab/AIL/code/requal.ctrl.conversion.cmds 
	fi
    else 
	if [ $n64 -eq 0 ]; then 
	    echo mv $bamfile ${bamfile/.bam/.requaled.bam} >> /group/palmer-lab/AIL/code/requal.ctrl.conversion.cmds
	else 
	    echo $bamfile Impossible $n32 $n64 $ns64 - so ignoring file.
	fi
    fi
done
## Run the conversion scripts
numfiles=`wc -l /group/palmer-lab/AIL/code/requal.ctrl.conversion.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmer-lab/AIL/code/requal.sh

#########################################################################################################
## Indel realignment steps - SHYAM 									#
#########################################################################################################
## In this case, we are using the target list created from the						#
## WT indel list. In case and if we ever get a list from HL						#
## that target list needs to be created again. 								#
## THAT STEP IS NOT INCLUDED HERE - SEE GATK RealignerTargetCreator					#
## Using the known indels - from the Wellcome Trust							#
## data, we can recalibrate the reads to avoid a lot							#
## false positives SNPs in the final data.							       	#
## ls /group/palmer-lab/AIL/GBS/bams/*.requaled.bam > /group/palmer-lab/AIL/GBS/bam.ail.requaled.list  	#
## Run the realignment code  									       	#
## numlines=`wc -l /group/palmer-lab/AIL/GBS/bam.ail.requaled.list | sed "s/^\W\+//" | cut -f1 -d " "` 	#
## qsub -t 1-${numlines}%200 /group/palmerlab/AIL/code/realign.sh                                      	#
#########################################################################################################


#####################################################
## 14 Jul 2015 - Natalia 
## Indel realignment steps using Lawson data

#####################################################
## RealignerTargetCreator - GATK 3.3-0
## Obtain list of intervals to perform realignment using Lawson indels as input.

# 1. Convert /AIL/knownSNPs/Lawson/LG.Indels & SM.Indels to LG_SM_Indels.vcf (Thank you April)
# 	April's script for doing this is in /AIL/LgSm-DataProcessing/scripts/make_indel_vcf.R

# 2. Sort vcf indel and SNP files using /shared_code/sortByRef.pl
#    Remember to restore VCF-style headers to the newly generated file.
perl /group/palmer-lab/shared_code/sortByRef.pl /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf /group/palmer-lab/reference_genomes/mouse/mm10.fasta.fai > LG_SM_sortedIndels.vc
perl /group/palmer-lab/shared_code/sortByRef.pl /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.lawson.vcf /group/palmer-lab/reference_genomes/mouse/mm10.fasta.fai > LGSM_mm10.orderedN.vcf
	
# 3. Run RealignerTargetCreator
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta -known:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf -o LG_SM.mm10.realign.intervals

#####################################################
## IndelRealigner - GATK 3.3-0
## Realign F1, LG and SM controls around known indels

# 1. Make .list files for GATK
## /group/palmer-lab/AIL/GBS/bams/realign.ail.1.list
## /group/palmer-lab/AIL/GBS/bams/realign.ail.2.list
## /group/palmer-lab/AIL/GBS/bam.ail.ctrl.list

# 2. Run IndelRealigner
numlines=`wc -l /group/palmer-lab/AIL/GBS/bam.ail.ctrl.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmerlab/AIL/code/realign.sh

####################################################
## CallableLoci - GATK 3.3-0
## What areas of the genome are considered callable in F1, LG and SM controls?

# 1. Make list file after shortening filenames
ls /group/palmer-lab/AIL/GBS/bams/controls/*rep[0-9].bam /group/palmer-lab/AIL/GBS/bams/controls/*rep[0-9][0-9].bam > /group/palmer-lab/AIL/GBS/bam.realigned.ctrl.list

# 2. Run CallableLoci (default settings)
numlines=`wc -l /group/palmer-lab/AIL/GBS/bam.realigned.ctrl.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmerlab/AIL/code/callable.loci.sh

###################################################
## AddOrReplaceReadGroups (Picard)
## BQSR is applied per lane; each lane should be a read group. The >@RG line in the bam file headers should 
## include the MACHINE:RUN:LANE fields from the sequence identifier. This should be identical to the PU field
## because Platform Unit takes precedence to RG if specified. The following file contains the information
## required to fill in the appropriate fields in the RG header for each bam file in bams:
## 		/group/palmer-lab/AIL/GBS/bams/RG_info.RData
## It contains the RG/PU fields, barcode index, sample ID, and input/output filenames. 
## I used an on-the-fly R script to make the commands to run Picard-AddOrReplaceReadGroups for each sample.
## Using the files 
##		/group/palmer-lab/AIL/code/addRG.sh 
##		/group/palmer-lab/code/addRepRGall.cmds 
## I changed the RG fields by running:
## 		qsub addRG.sh
## PU and RG are the same. LB = barcode index. SM = sample ID. CN = FGF for 152/154, GILAD for all others.
## IMPORTANT: I did this on all files with repN in the filename, not on the merged bam files (see Shyam's code
## below). Some samples were sequenced more than once and have filename 22222.rep0, 22222.rep1... etc. Shyam
## had merged files belonging to the same sample and kept the originals in the 'reps' folder. Since I want to
## run BQSR, and since BQSR is applied to each separate ReadGroup rather than each individual sample, I used
## the original 'rep' files in this application and in the ones that follow. The idea is that rep0 and rep1 for
## a given sample would have been sequenced in separate lanes and hence would have different RGs - they should
## not be merged before their RGs are changed to reflect that.

##################################################
## BaseRecalibrator - GATK 3.3-0
## BaseRecalibrator spits an error when I try to run more than two flowcells at a time. Since the command is
## applied to each lane separately, it doesn't really matter if I don't run all flowcells in one batch. It is
## faster (<24hrs per run) to submit a separate job for each flowcell, but there are multiple scripts involved.
## I made a bam.list file for each run (located in /bams/run290, e.g.) and fed them to GATK separately using 
## this command (for example):

## GENERATE FIRST PASS RECALIBRATION TABLE
## java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar 
## -T BaseRecalibrator 
## -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta 
## -I /group/palmer-lab/AIL/GBS/bqsr152.list 
## -knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf 
## -knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf
## -o bqsr152.table

## The files containing these commands are:
## /group/palmer-lab/AIL/GBS/bams/bqsr152.list		/group/palmer-lab/AIL/code/bqsr152.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr154.list		/group/palmer-lab/AIL/code/bqsr154.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr175.177.list 	/group/palmer-lab/AIL/code/bqsr175.177.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr183.list		/group/palmer-lab/AIL/code/bqsr183.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr195.list		/group/palmer-lab/AIL/code/bqsr195.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr290.list		/group/palmer-lab/AIL/code/bqsr290.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr291.list		/group/palmer-lab/AIL/code/bqsr291.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr381.list		/group/palmer-lab/AIL/code/bqsr381.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr386.list		/group/palmer-lab/AIL/code/bqsr386.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr388.list		/group/palmer-lab/AIL/code/bqsr388.sh
## /group/palmer-lab/AIL/GBS/bams/bqsr407.list		/group/palmer-lab/AIL/code/bqsr402.sh

## GENERATE SECOND PASS RECALIBRATION TABLE
## java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
## -T BaseRecalibrator
## -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta 
## -I /group/palmer-lab/AIL/GBS/bqsr152.list 
## -knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf 
## -knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson.LG_SM_Indels.vcf
## -BQSR bqsr152.table
## -o bqsr152b.table

## GENERATE PLOTS AND KEEP A COPY OF THE CSV
## java -Xmx4g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
## -T AnalyzeCovariates
## -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta 
## -before bqsr152.table
## -after bqsr152_b.table
## -csv result.bqsr152.csv  # optional
## -plots plot.bqsr152.pdf

#################################################
## PrintReads and Apply BQSR - GATK 3.3-0
## Use the BaseRecalibration table data to recalibrate QUAL scores in the input bam. The new QUAL
## score is the sum of the global difference between reported QUAL scores and the empirical QUAL +
## the quality bin-specific shift + the cycle*quality and dinucleotide*quality effects.
## PrintReads takes 1 or more bam files and input and outputs a single processed bam file. 

## java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
## -T PrintReads
## -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta 
## -I /group/palmer-lab/AIL/GBS/bqsr152.list 
## -BQSR bqsr152.table
## -o output.bam


##############################
# Merge bams for same sample #
##############################
## Merge the bam files for the samples that have 
## more than one legitimate bam file. 
## Find the samples that have more than one bam file 
## using 'rep' as the indicator for multiple bam files.
## Store the original bam files in the reps directory
if [ ! -d /group/palmer-lab/AIL/GBS/bams/reps ]; then 
  mkdir /group/palmer-lab/AIL/GBS/bams/reps
fi
ls /group/palmer-lab/AIL/GBS/bams/*rep*requaled.realign.bam | cut -f1 -d. | sort | uniq | xargs -i basename {} > /group/palmer-lab/AIL/GBS/ail.repeatsamples.list
for sample in `cat /group/palmer-lab/AIL/GBS/ail.repeatsamples.list`; do
  inputfiles=`ls /group/palmer-lab/AIL/GBS/bams/${sample}*bam | xargs -i echo -n "I={} " | sed "s/bams/bams\/reps/g"`
  mv /group/palmer-lab/AIL/GBS/bams/${sample}.rep* /group/palmer-lab/AIL/GBS/bams/reps/ # move the repfiles to the reps directory
  echo java -Xmx2g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar MergeSamFiles $inputfiles O=/group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bam SO=coordinate AS=true CREATE_INDEX=true 
done > /group/palmer-lab/AIL/code/merge.repeats.cmds
## Run the repeat commands
numlines=`wc -l /group/palmer-lab/AIL/code/merge.repeats.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines}%200 merge.ailbams.sh 
## Change the read group information for the merged bam file
## GATK assumes that each RGSM is a separate sample
## Each individual file has its own read group
## WAIT FOR MERGE TO FINISH BEFORE RUNNING RG CHANGE
for sample in `cat /group/palmer-lab/AIL/GBS/ail.repeatsamples.list`; do
  echo java -Xmx2g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar AddOrReplaceReadGroups I=/group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bam O=/group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.newRG.bam ID=$sample SM=$sample PL=illumina LB=merged PU=unit1 CREATE_INDEX=true
done > /group/palmer-lab/AIL/code/merge.changerg.cmds
## Run the RG change commands
numlines=`wc -l /group/palmer-lab/AIL/code/merge.changerg.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines}%200 merge.changerg.sh
## WAIT FOR THE END OF THE PREVIOUS SECTION BEFORE RUNNING THE MV LOOP
for sample in `cat /group/palmer-lab/AIL/GBS/ail.repeatsamples.list`; do
  mv /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.newRG.bam /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bam
  mv /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.newRG.bai /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bai
done

###########################
# Generate known vcf file #
###########################
## The vcf file with SNPs known to be different 
## between LG and SM. There are SNPs where the 
## reference allele in the VCF and the reference
## fasta are not matching. These SNPs are removed
## or flipped.
module load python
/group/palmer-lab/AIL/knownSNPs/flipIncorrectAlleles.py /group/palmer-lab/reference_genomes/mouse/mm10.fasta /group/palmer-lab/AIL/knownSNPs/LGSM.mm10.vcf /group/palmer-lab/AIL/knownSNPs/LGSM.mm10.flipremoved.vcf 

###########################
#   Variant calling step  #
###########################
## Using the list of known variants, estimate the genotypes 
## for the samples in the list of samples. This assumes that
## we have less than 2048 total samples, if not, we need to
## split into more than 2 parts. 
halfnumbams=`ls /group/palmer-lab/AIL/GBS/bams/*requaled.realign.bam | wc -l | sed "s/^\W\+//" | cut -f1 -d " "`
let halfnumbams=halfnumbams/2
ls /group/palmer-lab/AIL/GBS/bams/*requaled.realign.bam | head -$halfnumbams > /group/palmer-lab/AIL/GBS/bam.ail.1.list
let halfnumbams=halfnumbams+1
ls /group/palmer-lab/AIL/GBS/bams/*requaled.realign.bam | tail -n +$halfnumbams > /group/palmer-lab/AIL/GBS/bam.ail.2.list
for setnum in {1..2}; do
  for chrom in {1..19} X; do 
    echo java -Xmx9g -jar /apps/software/GenomeAnalysisTK/2.7-4/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -I /group/palmer-lab/AIL/GBS/bam.ail.$setnum.list -stand_call_conf 30.0 -stand_emit_conf 30.0 --heterozygosity 0.005 -L chr$chrom -gt_mode GENOTYPE_GIVEN_ALLELES -alleles /group/palmer-lab/AIL/knownSNPs/LGSM.mm10.flipremoved.vcf -o /group/palmer-lab/AIL/GBS/vars/ail.$setnum.chr$chrom.known.vcf
  done
done > /group/palmer-lab/AIL/code/ail.callvars.cmds
## Run the variant call command
numlines=`wc -l /group/palmer-lab/AIL/code/ail.callvars.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-$numlines /group/palmer-lab/AIL/code/runUG.ail.sh

#### To do for positive controls ####
# Check genotype concordance for reps
#java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar -T GenotypeConcordance -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -eval data -truth data -o filename -sites filename 


############################
# VCF merging for the sets #
############################
## Merge the vcf for all the samples
## across the two sets of samples
for chrom in {1..19} X; do
  echo java -Xmx2g -jar \$GATK -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -T CombineVariants --variant /group/palmer-lab/AIL/GBS/vars/ail.1.chr$chrom.known.vcf --variant /group/palmer-lab/AIL/GBS/vars/ail.2.chr$chrom.known.vcf -o /group/palmer-lab/AIL/GBS/vars/ail.chr$chrom.known.vcf -genotypeMergeOptions UNIQUIFY
done > /group/palmer-lab/AIL/code/ail.mergevcfs.cmds
## Run the merge commands
numlines=`wc -l /group/palmer-lab/AIL/code/ail.mergevcfs.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmer-lab/AIL/code/merge.ailvcfs.sh

##########################################
# Convert vcf to genprobs for imputation #
##########################################
## Use the mkPreImputeGenoFile to make preimpute files
if [ ! -d /group/palmer-lab/AIL/GBS/preimpute ]; then
  mkdir /group/palmer-lab/AIL/GBS/preimpute
fi
## Submit the preimpute commands using an inline script
echo "#! /bin/bash
#PBS -N convertImpute
#PBS -S /bin/bash
#PBS -l mem=3g
#PBS -t 1-20
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe

module load python/2.7.6
if [ \$PBS_ARRAYID -eq 20 ]; then
  chrom=X
else 
  chrom=\$PBS_ARRAYID
fi
python /group/palmer-lab/shared_code/makeImputeGenoFile.py -v /group/palmer-lab/AIL/GBS/vars/ail.chr\$chrom.known.vcf -c \$chrom -l 97 -g 4 -o /group/palmer-lab/AIL/GBS/preimpute/ail.chr\$chrom.preimpute.geno
" | qsub 
## Make the haplotype and legend files for panel 0 markers for impute2
if [ ! -d /group/palmer-lab/AIL/knownSNPs/imputeHaplotypes ]; then 
  mkdir /group/palmer-lab/AIL/knownSNPs/imputeHaplotypes
fi
module load python
python /group/palmer-lab/shared_code/createPanel0Haplotypes.py /group/palmer-lab/AIL/knownSNPs/LGSM.mm10.flipremoved.vcf /group/palmer-lab/AIL/knownSNPs/imputeHaplotypes

###########################################
# Run imputation using the preimpute geno #
###########################################
## make the imputation directory
if [ ! -d /group/palmer-lab/AIL/GBS/imputed ]; then
  mkdir /group/palmer-lab/AIL/GBS/imputed
fi
## make the imputation commands
## Constants that are used in the imputation process.
## window - size of region to impute - in base pairs
## buffer - size of region surrounding the window, used to improve imputation (in kb)
## numhaps - number of phasing haplotypes
## popsize - effective population size
window=2000000
buffer=500
numhaps=200
popsize=100
for chrom in {1..19} X; do 
  chrom=chr$chrom
  length=`grep -w $chrom /group/palmer-lab/reference_genomes/mouse/chrLengths_mm10_chrnames.txt | tr -s " " " " | cut -f2 -d " "`
  start=1
  end=0
  while [ $start -lt $length ]; do 
    let end=start+window-1
    if [ $end -gt $length ]; then
      end=$length
    fi
    echo /group/palmer-lab/tools/impute_v2.3.2_x86_64_static/impute2 -g /group/palmer-lab/AIL/GBS/preimpute/ail.$chrom.preimpute.geno -m /group/palmer-lab/reference_genomes/mouse/map/$chrom.mm10.impute2.map -l /group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/$chrom.legend -h /group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/$chrom.hap -int $start $end -Ne $popsize -buffer $buffer -allow_large_regions -k $numhaps -pgs_prob -prob_g -o /group/palmer-lab/AIL/GBS/imputed/$chrom.$start
    let start=start+window
  done
done > /group/palmer-lab/AIL/code/ail.impute.cmds
## Run the impute commands
numlines=`wc -l /group/palmer-lab/AIL/code/ail.impute.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines}%500 run.ailimpute.sh

##########################################################################################
# Run imputation using the preimpute geno -- IMPUTE MISSING GENOS WITHIN SAMPLES ONLY    #
# Preparing files for Natalia to map QTL with empirical SNPs only - preliminary analysis #
# to get a list of top SNPs for eQTL mapping                                             #
##########################################################################################
## make the imputation directory
if [ ! -d /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical ]; then
  mkdir /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical
fi
## make the imputation commands
## Constants that are used in the imputation process.
## window - size of region to impute - in base pairs
## buffer - size of region surrounding the window, used to improve imputation (in kb)
## numhaps - number of phasing haplotypes
## popsize - effective population size
window=2000000
buffer=500
numhaps=200
popsize=100
for chrom in {1..19} X; do 
  chrom=chr$chrom
  length=`grep -w $chrom /group/palmer-lab/reference_genomes/mouse/chrLengths_mm10_chrnames.txt | tr -s " " " " | cut -f2 -d " "`
  start=1
  end=0
  while [ $start -lt $length ]; do 
    let end=start+window-1
    if [ $end -gt $length ]; then
      end=$length
    fi
    echo /group/palmer-lab/tools/impute_v2.3.2_x86_64_static/impute2 -g /group/palmer-lab/AIL/GBS/preimpute/ail.$chrom.preimpute.geno -m /group/palmer-lab/reference_genomes/mouse/map/$chrom.mm10.impute2.map -int $start $end -Ne $popsize -buffer $buffer -allow_large_regions -k $numhaps -pgs_prob -prob_g -o /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical/$chrom.$start
    let start=start+window
  done
done > /group/palmer-lab/AIL/code/ail.impute.cmds
## Run the impute commands
numlines=`wc -l /group/palmer-lab/AIL/code/ail.impute.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines}%500 run.ailimpute.sh


#########################################
# Filter, convert to dosage and combine #
#########################################
if [ ! -d /group/palmer-lab/AIL/GBS/dosage/ ]; then
  mkdir /group/palmer-lab/AIL/GBS/dosage
fi
## Combine within chromosomes
for chrom in {1..19} X; do 
  echo "Starting chromosome $chrom."
  cat `ls /group/palmer-lab/AIL/GBS/imputed/chr$chrom.*[0-9] | sort -k2,2n -t.` > /group/palmer-lab/AIL/GBS/dosage/chr$chrom
  head -1 /group/palmer-lab/AIL/GBS/imputed/chr1.10000001_info > /group/palmer-lab/AIL/GBS/dosage/chr${chrom}_info
  cat `ls /group/palmer-lab/AIL/GBS/imputed/chr$chrom.*info | sort -k2,2n -t.` | grep -v snp >> /group/palmer-lab/AIL/GBS/dosage/chr${chrom}_info
  echo "Done with chromosome $chrom."
done
## Filter on info.
echo "#! /bin/bash
#PBS -N imputeFilter
#PBS -S /bin/bash
#PBS -l mem=16g
#PBS -t 1-20
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe

module load python/2.7.6
if [ \$PBS_ARRAYID -eq 20 ]; then
  chrom=X
else 
  chrom=\$PBS_ARRAYID
fi
python /group/palmer-lab/shared_code/filterAndConvert.py -f /group/palmer-lab/AIL/GBS/dosage/chr\$chrom -o /group/palmer-lab/AIL/GBS/dosage/chr\$chrom.filtered
echo Done with \$chrom
" | qsub
## Generate the leave-one-out kinship matrices
## Generate the input files - dont do this for X, 
## can use one of the autosomal things for X
for chrom in {1..19}; do 
  cat `ls /group/palmer-lab/AIL/GBS/dosage/chr[0-9]*.filtered.dosage | grep -vw chr$chrom |sort -k7.4,7n -t/` > /group/palmer-lab/AIL/GBS/dosage/chrNot$chrom.filtered.dosage
  echo "Done with $chrom"
done
## Generate dummy phenotype file
let numsamps=`head -1 /group/palmer-lab/AIL/GBS/dosage/chr1.filtered.dosage | wc -w`-3
rm -f /group/palmer-lab/AIL/GBS/dosage/dummy.pheno
while [ $numsamps -gt 0 ]; do 
  echo 1 >> /group/palmer-lab/AIL/GBS/dosage/dummy.pheno
  let numsamps=numsamps-1
done
## Run gemma to generate the LOO GRMs
echo "#! /bin/bash
#PBS -N computeKinship
#PBS -S /bin/bash
#PBS -l mem=20g
#PBS -t 1-19
#PBS -l walltime=12:00:00
#PBS -l nodes=1
#PBS -j oe

chrom=\$PBS_ARRAYID
module load gemma
cd /group/palmer-lab/AIL/GBS/dosage/
gemma -g /group/palmer-lab/AIL/GBS/dosage/chrNot\$chrom.filtered.dosage -p /group/palmer-lab/AIL/GBS/dosage/dummy.pheno -gk 1 -o chrNot\$chrom -maf 0.4 
echo Done estimating kinship matrices \$chrom
" | qsub
## clean up the temporary files - WAIT TILL THE gemma processes are done
rm /group/palmer-lab/AIL/GBS/dosage/dummy.pheno
rm /group/palmer-lab/AIL/GBS/dosage/chrNot*dosage
## move the kinship files to the correct place
if [ ! -d /group/palmer-lab/AIL/qtlmapping ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping
fi
if [ ! -d /group/palmer-lab/AIL/qtlmapping/kinship ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping/kinship
fi
mv /group/palmer-lab/AIL/GBS/dosage/output/chrNot* /group/palmer-lab/AIL/qtlmapping/kinship


##########################################################################################
# Filter, convert to dosage and combine -- IMPUTE MISSING GENOS WITHIN SAMPLES ONLY      #
# Preparing files for Natalia to map QTL with empirical SNPs only - preliminary analysis #
# to get a list of top SNPs for eQTL mapping                                             #
##########################################################################################

#########################################
if [ ! -d /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical ]; then
  mkdir /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical
fi
## Combine within chromosomes
for chrom in {1..19} X; do 
  echo "Starting chromosome $chrom."
  cat `ls /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical/chr$chrom.*[0-9] | sort -k2,2n -t.` > /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr$chrom
  head -1 /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical/chr1.10000001_info > /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr${chrom}_info
  cat `ls /group/palmer-lab/AIL/GBS/imputed/onlyEmpirical/chr$chrom.*info | sort -k2,2n -t.` | grep -v snp >> /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr${chrom}_info
  echo "Done with chromosome $chrom."
done
## Filter on info.
echo "#! /bin/bash
#PBS -N imputeFilter
#PBS -S /bin/bash
#PBS -l mem=16g
#PBS -t 1-20
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe

module load python/2.7.6
if [ \$PBS_ARRAYID -eq 20 ]; then
  chrom=X
else 
  chrom=\$PBS_ARRAYID
fi
python /group/palmer-lab/shared_code/filterAndConvert.py -f /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr\$chrom -o /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr\$chrom.filtered
echo Done with \$chrom
" | qsub
## Generate the leave-one-out kinship matrices
## Generate the input files - dont do this for X, 
## can use one of the autosomal things for X
for chrom in {1..19}; do 
  cat `ls /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr[0-9]*.filtered.dosage | grep -vw chr$chrom |sort -k7.4,7n -t/` > /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrNot$chrom.filtered.dosage
  echo "Done with $chrom"
done
## Generate dummy phenotype file
let numsamps=`head -1 /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr1.filtered.dosage | wc -w`-3
rm -f /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/dummy.pheno
while [ $numsamps -gt 0 ]; do 
  echo 1 >> /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/dummy.pheno
  let numsamps=numsamps-1
done
## Run gemma to generate the LOO GRMs
echo "#! /bin/bash
#PBS -N computeKinship
#PBS -S /bin/bash
#PBS -l mem=20g
#PBS -t 1-19
#PBS -l walltime=12:00:00
#PBS -l nodes=1
#PBS -j oe

chrom=\$PBS_ARRAYID
module load gemma
cd /group/palmer-lab/AIL/GBS/dosage/
gemma -g /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrNot\$chrom.filtered.dosage -p /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/dummy.pheno -gk 1 -o chrNot\$chrom -maf 0.2 
echo Done estimating kinship matrices \$chrom
" | qsub
## clean up the temporary files - WAIT TILL THE gemma processes are done
rm /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/dummy.pheno
rm /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrNot*dosage
## move the kinship files to the correct place
if [ ! -d /group/palmer-lab/AIL/qtlmapping/ ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping/
fi
if [ ! -d /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical
fi
mv /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/output/chrNot* /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical

#


################################
# Mapping of traits using GBS  #
################################
if [ ! -d /group/palmer-lab/AIL/qtlmapping/covariates ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping/covariates
fi
## Make the file with the genotyped samples list
head -10000 /group/palmer-lab/AIL/GBS/vars/ail.chr1.known.vcf | grep CHROM | cut -f10- | sed "s/\.variant2\{0,1\}//g" | tr -s "\t" "\n" > /group/palmer-lab/AIL/GBS/dosage/genotyped.samples.txt
## Run the R script to make all the qtl gemma commands and
## appropriate covariate file.
module load R
R CMD BATCH /group/palmer-lab/AIL/code/make.ail.qtl.scripts.R 
## Run the gemma script
numlines=`wc -l /group/palmer-lab/AIL/code/gemma.alltraits.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t1-${numlines}%125 /group/palmer-lab/AIL/code/run.gemma.sh

######################################
# Compute heritabilities using gemma #
######################################
if [ ! -d /group/palmer-lab/AIL/qtlmapping/heritabilities ]; then
  mkdir /group/palmer-lab/AIL/qtlmapping/heritabilities
fi
## Make the genomewide kinship matrix
cat `ls /group/palmer-lab/AIL/GBS/dosage/chr[0-9]*.filtered.dosage |sort -k7.4,7n -t/` > /group/palmer-lab/AIL/GBS/dosage/chrAll.filtered.dosage
echo "Done with concating dosages."
## Generate dummy phenotype file
let numsamps=`head -1 /group/palmer-lab/AIL/GBS/dosage/chr1.filtered.dosage | wc -w`-3
rm -f /group/palmer-lab/AIL/qtlmapping/heritabilities/dummy.pheno
while [ $numsamps -gt 0 ]; do 
  echo 1 >> /group/palmer-lab/AIL/qtlmapping/heritabilities/dummy.pheno
  let numsamps=numsamps-1
done
## Run gemma to generate the LOO GRMs
echo "#! /bin/bash
#PBS -N computeKinship
#PBS -S /bin/bash
#PBS -l mem=20g
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -j oe

module load gemma
cd /group/palmer-lab/AIL/qtlmapping/heritabilities
gemma -g /group/palmer-lab/AIL/GBS/dosage/chrAll.filtered.dosage -p /group/palmer-lab/AIL/qtlmapping/heritabilities/dummy.pheno -gk 1 -o chrAll -maf 0.4 
echo Done estimating kinship matrices using all autosomes
" | qsub
## Move the kinship files to the correct place
mv /group/palmer-lab/AIL/qtlmapping/heritabilities/output/chrAll.cXX.txt /group/palmer-lab/AIL/qtlmapping/heritabilities/kinship
rm -rf /group/palmer-lab/AIL/qtlmapping/heritabilities/output 
## Make the dummy dosage file.
head -100 /group/palmer-lab/AIL/GBS/dosage/chrAll.filtered.dosage > /group/palmer-lab/AIL/qtlmapping/heritabilities/dummy.dosage
## Run the R script to make all the heritability gemma commands
## This is done in the previous section. Cmds files made at the 
## same time as the gemma command file. If it does not exist throw
## throw a warning.
if [ ! -s /group/palmer-lab/AIL/code/gemma.herit.cmds ]; then
  echo "The list of commands file - gemma.herit.cmds - is not present or"
  echo "is empty. Please generate this file using the make.ail.qtl.scripts.R"
  echo "script file."
else
  ## Run the gemma script
  numlines=`wc -l /group/palmer-lab/AIL/code/gemma.herit.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
  qsub -t1-${numlines}%125 /group/palmer-lab/AIL/code/run.herit.sh
fi

#######################################################################
##### Make genome wide kinship matrix for onlyEmpirical SNPs      #####
##### So Natalia can get genetic correlations for multivar traits #####
#######################################################################

## Make the genomewide kinship matrix
cat `ls /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr[0-9]*.filtered.dosage |sort -k7.4,7n -t/` > /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrAll.filtered.dosage
echo "Done with concating dosages."
## Generate dummy phenotype file
let numsamps=`head -1 /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr1.filtered.dosage | wc -w`-3
rm -f /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical/dummy.pheno
while [ $numsamps -gt 0 ]; do 
  echo 1 >> /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical/dummy.pheno
  let numsamps=numsamps-1
done
## Run gemma to generate the genome-wide GRM
echo "#! /bin/bash
#PBS -N computeGWkinship
#PBS -S /bin/bash
#PBS -l mem=20gb
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -j oe

module load gemma
cd /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical
gemma -g /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrAll.filtered.dosage -p /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical/dummy.pheno -gk 1 -o chrAll -maf 0.2 
echo Done estimating kinship matrices using all autosomes
" | qsub


