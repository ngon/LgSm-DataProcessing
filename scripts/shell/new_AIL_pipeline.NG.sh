## -----------------------------------------------------------------------------------------------
## TITLE: AIL GBS Pipeline
## AUTHOR: N Gonzales & S Gopalakrishnan
## DATE: September 3, 2015
## DESCRIPTION: Alignment, requaling, realignment, merging, variant calling, imputation, GWAS, etc.


## WAIT FOR THE COMMANDS IN EACH SECTION TO FINISH RUNNING.



## ALIGNMENT --------------------------------------------------------------------------------------
## Input: demultiplexed, gzipped fastq files
## Requirements: bwa

## Create the file with the list of fastqs that need to be aligned.
ls /group/palmer-lab/AIL/GBS/fastqs/*fq.gz > /group/palmer-lab/AIL/GBS/fastqs.ail.list

## Submit it as a task array in a shell script
## example: /group/palmer-lab/AIL/code/alignToMM10Mem.sh (author: SG)
numfiles=`wc -l /group/palmer-lab/AIL/GBS/fastqs.ail.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numfiles} /group/palmer-lab/AIL/code/alignToMM10Mem.sh


## STANDARDIZE PHRED QUALITY SCORES (PL) -----------------------------------------------------------
## Input: aligned, coord-sorted .bam files
## Requirements: samtools, convertQualScoreBam

## Notes
## I didn't have execution permission on the original convertQualScoreBam file so I had to make a copy
## which I called convertQual_tmp in order to run it.

## Determine whether the file is encoded in Illumina 1.3 or 1.5 (phred+64) or Illumina 1.9 (phred+33)
## by looking at the first 100k reads of each file. Convert to phred+33 format.
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
## example: /group/palmer-lab/AIL/code/requal.sh
numfiles=`wc -l /group/palmer-lab/AIL/code/requal.ctrl.conversion.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmer-lab/AIL/code/requal.sh

## Check success by running
module load samtools
samtools view {filename} | cut -f 11 | head -10

## INDEL REALIGNMENT ------------------------------------------------------------------------------
## Input: Aligned, coord-sorted, requaled .bams
## Requirements: GATK v.3, interval list from Heather Lawson's LG/SM WGS data, mm10 dictionary

## HOW THE REQUIRED FILES WERE GENERATED
## April converted text files of LG and SM variants into VCF files for indels and SNPs using
## /AIL/LgSm-DataProcessing/scripts/make_indel_vcf.R. NG coordinate-sorted them and added VCF-style
## headers using /shared_code/sortByRef.pl. Examples:
perl /group/palmer-lab/shared_code/sortByRef.pl /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf /group/palmer-lab/reference_genomes/mouse/mm10.fasta.fai > LG_SM_sortedIndels.vc
perl /group/palmer-lab/shared_code/sortByRef.pl /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.lawson.vcf /group/palmer-lab/reference_genomes/mouse/mm10.fasta.fai > LGSM_mm10.orderedN.vcf
## NG ran GATK: RealignerTargetCreator to make interval list
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta -known:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf -o LG_SM.mm10.realign.intervals

## GATK: IndelRealigner
## Make list files (input) for GATK
numlines=`wc -l /group/palmer-lab/AIL/GBS/bam.ail.ctrl.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmerlab/AIL/code/realign.sh


## EDIT READ GROUP INFO ---------------------------------------------------------------------------
## Input: Aligned, coord-sorted, requaled, realigned .bams
## Requirements: picard: AddOrReplaceReadGroups, R, samtools, RG_info.RData

## Notes
## Base recalibration (next step) is applied per lane (read group = RG = lane). RG info was lost when
## samples were demuxed, so we have to add this information back in. The '@RG' line in the bam file
## headers should include the MACHINE:RUN:LANE fields from the sequence identifier. This should be
## identical to the PU (platform unit) field because in GATK v3, PU takes precedence over RG if both
## are specified. /group/palmer-lab/AIL/GBS/bams/readGroups/RG_info.RData contains information for
## all fields (RG, PU, LB, CN, SM). Picard commands were generated from this information using a for
## loop in R. LB = barcode index. SM = sample ID. CN = FGF for run152/154, GILAD for all others.

## Using the files
##		/group/palmer-lab/AIL/bams/readGroup/addRG.sh
##		/group/palmer-lab/AIL/GBS/bams/readGroup/addRepRGall.cmds
## I changed the RG fields by running
## 		qsub -t 1-N addRG.sh

## Validate files using
module load samtools
samtools view -H {filename}

## BASE RECALIBRATION -----------------------------------------------------------------------------
## Input: Input: Aligned, coord-sorted, requaled, realigned, regrouped .bams
## Requirements: GATK v.3, samples.list files, R (ggplot2, )

## Notes
## Run each flow cell separately because GATK runs exponentially slower as the number of read groups
## increases. Samples for each flow cell should be in a separate folder (run152, run381... etc). It
## is a good idea to validate bam files with picard: ValidateSamFile before running BQSR as a sanity
## check to make sure everything worked properly in previous steps. This process is very fast.
## There are 4 steps to Base Recalibration.

## 1. GATK: BaseRecalibrator (Round 1)
## Generating a first-pass recalibration table requires a command with this structure:
java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
-T BaseRecalibrator
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-I /group/palmer-lab/AIL/GBS/bqsr152.list
-knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf
-knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LG_SM_Indels.vcf
-o bqsr152.table

## 2. GATK: BaseRecalibrator (Round 2)
## Generate second-pass recalibration table
java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
-T BaseRecalibrator
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-I /group/palmer-lab/AIL/GBS/bqsr152.list
-knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf
-knownSites:VCF /group/palmer-lab/AIL/knownSNPs/Lawson.LG_SM_Indels.vcf
-BQSR bqsr152.table
-o bqsr152b.table

## 3. GATK: AnalyzeCovariates
## Make plots comparing bases before and after BQSR. Saving the csv is optional.
java -Xmx4g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
-T AnalyzeCovariates
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-before bqsr152.table
-after bqsr152_b.table
-csv result.bqsr152.csv
-plots plot.bqsr152.pdf

## 4. GATK: PrintReads
## Use the BaseRecalibration table data to recalibrate QUAL scores in each bam. The new QUAL
## score is the sum of the global difference between reported QUAL scores and the empirical QUAL +
## the quality bin-specific shift + the cycle*quality and dinucleotide*quality effects.
## PrintReads takes 1 or more bam files and input and outputs a single processed bam file.
java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar
-T PrintReads
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-I /group/palmer-lab/AIL/GBS/bqsr152.list
-BQSR bqsr152.table
-o output.bam


## MERGE BAMS FOR REP SAMPLES ---------------------------------------------------------------------
## Input: Aligned, coord-sorted, requaled, realigned, regrouped bams that were sequenced >1x
## Requirements: picard

## Notes
## Some samples were sequenced on more than one flow cell lane. Merge them before feeding them to
## the variant caller.

## Create a directory for repeat samples
if [ ! -d /group/palmer-lab/AIL/GBS/bams/reps ]; then
    mkdir /group/palmer-lab/AIL/GBS/bams/reps
fi
## Find samples with 'rep' in the filename and move them to the rep directory.
ls /group/palmer-lab/AIL/GBS/bams/run*/rep*.bam | cut -f1 -d. | sort | uniq | xargs -i basename {} > /group/palmer-lab/AIL/GBS/ail.repeatsamples.list
## Before running this script, test the command within 'inputfiles' using the shell function 'find'.
## Make sure the file extensions match your input.
for sample in `cat /group/palmer-lab/AIL/GBS/ail.repeatsamples.list`; do
    inputfiles=`ls /group/palmer-lab/AIL/GBS/bams/${sample}*bam | xargs -i echo -n "I={} " | sed "s/bams/bams\/reps/g"`
    mv /group/palmer-lab/AIL/GBS/bams/${sample}.rep* /group/palmer-lab/AIL/GBS/bams/reps/
    echo java -Xmx2g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar MergeSamFiles $inputfiles O=/group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bam SO=coordinate AS=true CREATE_INDEX=true
done > /group/palmer-lab/AIL/code/merge.repeats.cmds

## Run the merge commands
numlines=`wc -l /group/palmer-lab/AIL/code/merge.repeats.cmds | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} merge.ailbams.sh

## Wait for the previous commands to finish before running the mv loop
for sample in `cat /group/palmer-lab/AIL/GBS/ail.repeatsamples.list`; do
    mv /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.newRG.bam /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bam
    mv /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.newRG.bai /group/palmer-lab/AIL/GBS/bams/${sample}.rg.requaled.realign.bai
done


## RE-REALIGN MERGED SAMPLES -------------------------------------------------------------------------
## Input: Aligned, coord-sorted, requaled, realigned, regrouped bams that were merged
## Requirements: GATK v.3, interval list from Heather Lawson's LG/SM WGS data, mm10 dictionary

## Note: This is optional, but it ensures consistent alignment across all RGs for a sample.

## GATK: IndelRealigner
## Make list files (input) for GATK
numlines=`wc -l /group/palmer-lab/AIL/GBS/merged.reps.list | sed "s/^\W\+//" | cut -f1 -d " "`
qsub -t 1-${numlines} /group/palmerlab/AIL/code/realign.sh


## PER-SAMPLE VARIANT CALLING ---------------------------------------------------------------------
## Input: Aligned, coord-sorted, requaled, realigned, regrouped bams, 1 file per sample.
## Requirements: GATK v.3

INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/vars/alleleCalls.list | tail -1`
BASENAME=`basename $INFILE`
OUTDIR="/group/palmer-lab/AIL/GBS/vars/rawAllelesHC/run152"

 java -Xmx9g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
 -T HaplotypeCaller
 -R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
 -I $INFILE
 -o $OUTDIR/$BASENAME.rawVars.params.g.vcf
 --dbsnp /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf
 --emitRefConfidence GVCF
 --output_mode EMIT_ALL_SITES
 --genotyping_mode DISCOVERY
 --minPruning 1
 --minDanglingBranchLength 1
 --max_alternate_alleles 2
 --minReadsPerAlignmentStart 1
 --heterozygosity 0.005

## MERGE gVCFs IN BATCHES OF 200 -------------------------------------------------------------------
## Input: Per-sample gVCFs and corresponding .list file
## Requirements: GATK, bcftools

## GATK: CombineGVCFs
## /vars/rawAllelesHC/mergeGVCF.sh
java -Xmx3g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jarqsub
-T CombineGVCFs
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-o allelesValidationCohort.g.vcf
--variant validationGvcfs.list

## Set variant IDs
module load bcftools/1.2
bcftools annotate --set-id +'ail\.%CHROM\.%POS' $INFILE

## GENOTYPE CALLING
java -Xmx6g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T GenotypeGVCFs
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-o genotypesValidationCohort.vcf
-variant allelesValidationCohort.g.vcf

## VARIANT ANNOTATION
## VARIANT FILTERING (VQSR) FOR SNPS AND INDELS


## NOTES FROM BROAD INST MAR 15 WORKSHOP -----------------------------------------------------------

## Joint discovery ensures info is available for every site rather than for just vars.
## Database of likelihoods: prob of a geno being HOM-ALT or HET for each individual
## Result of GenoCalling is a squared off matrix for analysis. A column for each site, variant, and
## then the genotype for each sample along with phred-scaled likelihoods given NGS data in
## hom/het/hom-alt format.
## If no coverage at a particular site, joint genotyping will ignore it
##
## To determine active regions, HC traverses the genome with a sliding window, counting mismatches,
## indels and soft clips. Once it reaches a threshold (gathers enough evidence for variation in a
## region), it stops and assembles plausible haplotypes in the active region using De Bruijn-like
## graph assembly. The reads at each locus are chopped up into smaller sized kmers. When kmers
## overlap, this defines an edge in the De Bruijn graph. Edges are weighted according to how many
## reads support the overlap. Haplotypes are assembled by traversing the graph.
## Min pruning: gets rid of edges that don't have enough reads (weights) supporting overlap.
## Haplotypes are then realigned to the reference (reassembly).
##
## HMM: What is the likelihood of the haplotype given the reads? Empirical gap penalties are derived
## from BQSR. Base mismatch penalties are derived from base quality scores.
## Then each sample is genotyped (2 steps). First get likelihoods of alleles given the reads and
## then apply it to a Bayesian model. Take the highest probability of haplotypes given the reads
## (result of the HMM) that contain the allele (for each variant position). For example, suppose
## the variant is either a T or a gap, and there are 4 possible haplotypes. Only Hap 1 and Hap 3
## contain a gap. HC takes Hap 3, e.g., because given that there are 2 reads,it has the highest
## probability associated with it.
## find the probability of each genotype given the reads.
## bonus perk of haplotype calling: free physical phasing.
## PID - phase identifier; PID - phase genotype. If phased, you see a pipe | in the gVCF.
## For multisample analysis, you need likelihoods of all possible genotypes at all sites (Ref conf
## model) - even nonvariant sites.
##
## gVCF files include a <non-ref> allele annotation for potentially non-called alternative alleles.
## start, end and site level annotations. PL (genotype confidence)
##
## VQSR NOTES
## Main idea is that by building a model of what true genetic variation looks like, one can rank-
## order variants based on their likelihood of being real.
##
## Use truth sites to estimate sensitivity
## Ti/Tv ratio is used to estimate the proportion of true:false positives in the callset (based on
## expectation of what the Ti/Tv ratio should be)
## you want to pick the tranche that has mostly true positives and only a few false positives.
## The first tranche (from the bottom, with the highest value of truth
## sensitivity but usually the lowest values of novel Ti/Tv) is exceedingly specific but less
## sensitive, and each subsequent tranche in turn introduces additional true
## positive calls along with a growing number of false positive calls.
## a higher priority than one can dip down into lower tranches. These tranches are applied to the
## output VCF using the FILTER field. you can choose to use some of the filtered records or only use
## the PASSing records.
## VQSR only uses annotations that are available for all of the samples.
##
## Math details: uses the variational Bayes EM algorithm to fit a mixture of Gaussians on a training
## set of annotated variants. each gaussian gets a weight according to the density of points and the
## parameters are optimized in a bayesian iterative process.
## MAP probability = maximum a-posteriori probability
##














