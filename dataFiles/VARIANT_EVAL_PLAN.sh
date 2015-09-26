Evaluating output of haplotype caller

SNPs only
Exclude indels with GATK-SelectVariants

Venn Diagrams of sites covered by various resources.
- SNPs in Lawson vs. Cheverud: they should not be very different. It would be interesting if they
  were. In that case the known SNP file should be the Intersection of the 2 sets
- If concordant as expected, use Union(Lawson + Cheverud + MUGA) as KNOWN

- KNOWN vs. GBS
- For all subsets: check sex.
    - subset I[nbred]: LG and SM samples from GBS and LAWSON BAMS, F1 from GBS
        * Concordance between LG/LG, SM/SM, and F1/F1 GBS files (there are 2 biological reps for
          each genetic background)
            > % homozygous sites in F1s, % heterozygous sites in LG and SM
        * Concordance between KNOWN and LAWSON BAMS, LAWSON BAMS and Subset I
            > Running HaplotypeCaller on LAWSON BAMS to see how its genotype calls match up to
              the ones filtered and considered True by HL. This could be useful as a sanity check
              for recalibrating variants. If the KNOWN set is used as the truth set, most of these
              calls should pass. HLs paper will have a tally of how many SNPs were filtered out of
              the original BAMs, which will be useful. I can use her KNOWN SNPs as the truth set
              for recalibrating bases with BQSR and then run those samples through the same process
              as the GBS BAMs.
    - subset A[ils genotyped on MUGA]: comparison of GBS and MUGA data
        * Overall coverage distribution at sites overlapping with KNOWN variants
        * Overall coverage distribution at sites not overlapping with KNOWN variants.
        * PerSample: % concordant genotypes between sites covered by both GBS and MUGA
        * PerSample: total DP (and DP at each allele if heterozygous) at each site
        * PerSample: DP by quality score at each site
        * PerSample: % heterozygous, hom-ref, and hom-alt sites
        * PerSite: average depth and quality score
        * PerSite: average heterozygosity
        * PerSibPair: compare genetic kinship between sibs within and between genotyping platforms
            > i.e. Sib1,2-MUGA, Sib1,2-GBS, Sib1-GBS,Sib2-MUGA, Sib1-MUGA,Sib2-GBS
            > kinship coefficient and %IBS(0) SNPs
            > also compare kinship and IBS(0) for a few unrelated pairs
    - subset F[amily of AILs from F39-43]: GBS
        * Overall coverage distribution at sites overlapping with KNOWN variants
        * Overall coverage distribution at sites not overlapping with KNOWN variants.
        * PerSample: total DP (and DP at each allele if heterozygous) at each site
        * PerSample: DP by quality score at each site
        * PerSample: % heterozygous, hom-ref, and hom-alt sites
        * PerSite: average depth and quality score
        * PerSite: average heterozygosity
        * Make a small pedigree file (it has to include all samples in the gVCF)
        * Mendelian error checking with GATK
            > % variants violating Mendelian expectations
            > it is possible to check for de Novo variants - could be useful to see how this
              function is working and in a few cases could provide additional evidence for
              discovery calls.
        * PhaseByTransmission
        * Theoretical vs. Realized Kinship
        * Annotate these sites with Inbreeding Coefficient

Error analysis
- What is the read depth?
- What is the average base quality?
- What is the Ts/Tv ratio?
- What is the nucleotide context?
- e.g. for...
    * sites where more than 2 alternative alleles are acknowledged (see NDA field - these are
      present but not necessarily called
    * Mendelian errors identified from Subset F
    * discordant SNPs identified from Subsets I and A


### PREPARING HEATHERS BAM FILES FOR GATK PROCESSING ---------------------------------------------

# BAM files don't have the proper contig labels (i.e. '1' is listed instead of 'chr1')
> module load samtools
> samtools view -H SM_rmdup_sorted_hl.bam > SM_header
> samtools view -H LG_rmdup_sorted_hl.bam > LG_header
# vim: %s/SN:/SN:chr/g
> samtools reheader SM_header SM_rmdup_sorted_hl.bam > SM.rehead.bam
> samtools reheader LG_header LG_rmdup_sorted_hl.bam > LG.rehead.bam
# create list of unique reads in Lawson bams.
> module load samtools
> samtools view SM_rmdup_sorted_hl.bam | cut -f1 | sort -u > sm.readNames
> samtools view LG_rmdup_sorted_hl.bam | cut -f1 | sort -u > lg.readNames
# find unique read groups (first 4 fields of readName)
> cut -f1-4 -d ":" sm.readNames | uniq > sm.uniq.RGs
> cut -f1-4 -d ":" lg.readNames | uniq > lg.uniq.RGs
# make read list files for filtering. takes lots of memory so i'll process one RG in chunks to
# see how it goes. sm file is > 365m reads, lg is > 457m reads
head -1000000 sm.readNames | grep '^HWI-ST157:570:H11KYADXX:1' > sm.run570.1
sed -n '1000001,10000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:1' >> sm.run570.1
sed -n '10000001,45000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:1' >> sm.run570.1
sed -n '45000001,50000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:1' >> sm.run570.1
sed -n '50000001,70000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:1' >> sm.run570.1
sed -n '50000001,70000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:2' > sm.run570.2
sed -n '70000001,120000000p' sm.readNames | grep '^HWI-ST157:570:H11KYADXX:2' >> sm.run570.2
sed -n '70000001,120000000p' sm.readNames | grep '^HWI-ST157:571:H12JYADXX:1' >> sm.run571.1
sed -n '120000001,150000000p' sm.readNames | grep '^HWI-ST157:571:H12JYADXX:1' >> sm.run571.1
sed -n '120000001,170000000p' sm.readNames | grep '^HWI-ST157:571:H12JYADXX:2' >> sm.run571.2
sed -n '170000001,200000000p' sm.readNames | grep '^HWI-ST157:571:H12JYADXX:2' >> sm.run571.2
sed -n '170000001,250000000p' sm.readNames | grep '^HWI-ST157:572:H0EKHADXX:1' >> sm.run572.1
sed -n '170000001,250000000p' sm.readNames | grep '^HWI-ST157:572:H0EKHADXX:2' >> sm.run572.2
sed -n '250000000,300000000p' sm.readNames | grep '^HWI-ST157:572:H0EKHADXX:2' >> sm.run572.2
sed -n '250000000,300000000p' sm.readNames | grep '^HWI-ST157:572:H0EKHADXX:2' >> sm.run572.2
sed -n '250000000,320000000p' sm.readNames | grep '^HW-ST997:375:H76K3ADXX:1' >> sm.run375.1
sed -n '250000000,320000000p' sm.readNames | grep '^HW-ST997:375:H76K3ADXX:2' >> sm.run375.2
sed -n '320000000,365018164p' sm.readNames | grep '^HW-ST997:375:H76K3ADXX:2' >> sm.run375.2
# number of SM reads in each lane
#    52211323 sm.run375.1
#    50419795 sm.run375.2
#    50750778 sm.run570.1
#    47632603 sm.run570.2
#    39590570 sm.run571.1
#    39047266 sm.run571.2
#    42841471 sm.run572.1
#    42524360 sm.run572.2
#    365 018 166 total (2 extra. one of the files may have 2 blank lines)

# overlap between Lawson and Cheverud LGSM VCFs
bedtools intersect -a /group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf
-b /group/palmer-lab/AIL/knownSNPs/LGSM.mm10.flipremoved.vcf -wa -wb
> /group/palmer-lab/AIL/GBS/vars/siteOverlap_LGSM

# how many snps in lawson data
sed -n '71,9456361p' LGSM.mm10.orderedN.vcf | cut -f2 | uniq -u | wc -l       # 4790363
# how many snps just on the x chr in lawson data
sed -n '9273737,9456361p' LGSM.mm10.orderedN.vcf | cut -f2 | uniq -u | wc -l  # 127847
# how many snps in jim's data (autosomes only)
sed -n '71,4297119p' LGSM.mm10.flipremoved.vcf | cut -f2 | uniq -u | wc -l    # 4297049

# filter reads
java -Xmx3g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar FilterSamReads I=SM.rehead.bam
FILTER=includeReadList
READ_LIST_FILE=sm.run571.1
O=/group/palmer-lab/AIL/knownSNPs/Lawson/sm.run571.1.bam
SORT_ORDER=coordinate
CREATE_INDEX=true

# reformat RG line in headers
cmds <- c()
for (i in 1:25){
cmds <- c(cmds, paste0("java -Xmx2g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar AddOrReplaceReadGroups INPUT=/group/palmer-lab/AIL/GBS/bams/", rgRedo[i,6], " OUTPUT=", rgRedo[i,1], " RGID=", rgRedo[i,2], " RGSM=", rgRedo[i,4], " RGPL=ILLUMINA RGLB=", rgRedo[i,5], " RGPU=", rgRedo[i,2], " RGCN=FGF CREATE_INDEX=true"))
}

# validate bams
java -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar ValidateSamFile INPUT='/group/palmer-lab/AIL/GBS/bams/run152/50736.rep0.forBQSR.ra.bam' OUTPUT=/group/palmer-lab/AIL/GBS/bams/run152/report


###### COMBINE GVCFS ------------------------------------------------------------------------------
## /vars/rawAllelesHC/mergeGVCF.sh
java -Xmx3g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T CombineGVCFs
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-o allelesValidationCohort.g.vcf
--variant validationGvcfs.list

## Set variant IDs (/vars/rawAllelesHC/setVarIDs.sh)
cd /group/palmer-lab/AIL/GBS/vars/rawAllelesHC/
module load bcftools/1.2
INFILE=`head -$PBS_ARRAY_ID /group/palmer-lab/AIL/GBS/vars/rawAllelesHC/gigaMiceGVCFs.list | tail -1`
BASE=`head -$PBS_ARRAY_ID /group/palmer-lab/AIL/GBS/vars/rawAllelesHC/gigaMiceGVCFs.base | tail -1`
echo "Creating IDs for each variant: "
bcftools annotate --set-id +'ail\.%CHROM\.%POS' $INFILE -o $BASE.id.g.vcf
echo "Done naming variants."


## CALL GENOTYPES ---------------------------------------------------------------------------------
java -Xmx6g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T GenotypeGVCFs
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
--variant micetoGenotype.list
-hets 0.005
-stand_call_conf 20.0
-stand_emit_conf 20.0
-maxAltAlleles 2
-A GCContent -A HardyWeinberg
-o genotypesValCohortSet1.vcf
### first pass done for 70 samples currently available


## ANNOTATE VARIANTS ------------------------------------------------------------------------------
## Annotations present in the cohort vcf (LowQual Filter is in place): AD,DP,GQ,GT,MIN_DP,PGT,PID,
## PL, RGQ, SB, AC,AF(allele freq for each alt allele),AN, BaseQRankSum,ClippingRankSum,DS,END,FS,
## HaplotypeScore, InbreedingCoeff, MLEAC, MLEAF, MQ, MQRankSum,NDA,QD, ReadPosRankSum,SOR
##  DP=10;MLEAC=2,0;MLEAF=1.00,0.00;MQ=58.67;NDA=2  GT:AD:DP:GQ:PGT:PID:PL:SB

 java -Xmx6g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T VariantAnnotator
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
# for ped: -A PossibleDeNovo -A InbreedingCoeff
-o annotGenotypesValidationCohort.vcf


## VQSR -------------------------------------------------------------------------------------------
# downloaded Wellcome Trust dbSNP SNP and indel VCFs to knownSNPs/dbSNP. the contigs in this file
# don't match the ones in our reference dictionary.

java -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar UpdateVcfSequenceDictionary
INPUT='./mgp.v5.merged.snps_all.dbSNP142.vcf'
OUTPUT='./mgp.v5.allSNPs.update.dbSNP142.vcf'
SEQUENCE_DICTIONARY='/group/palmer-lab/reference_genomes/mouse/mm10.dict'

# one liner that adds 'chr' to chromosome numbers so picard can sort the file.
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' mgp.v5.allSNPs.update.dbSNP142.vcf > tmp && mv tmp mgp.v5.allSNPs.update.dbSNP142.vcf
sed -i 's/MT/chrM/g' mgp.v5.allSNPs.update.dbSNP142.vcf

java -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar SortVcf
INPUT='./mgp.v5.allSNPs.update.dbSNP142.vcf'
OUTPUT='./mgp.v5.allSNPs.updateSorted.dbSNP142.vcf'
SEQUENCE_DICTIONARY='/group/palmer-lab/reference_genomes/mouse/mm10.dict'

# alternative (above isn't working): download mm138 build from NCBI instead.
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/snp138.txt.gz
gunzip snp138.txt.gz
/apps/software/bcftools/1.2/bin/vcfutils.pl ucscsnp2vcf snp138.txt > snp138.ucsc.vcf


java -Xmx6g -jar
/apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T VariantRecalibrator
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-input genotypesValidationCohort.vcf
-resource:lawson,known=false,training=true,truth=true,
    prior=15.0,/group/palmer-lab/AIL/knownSNPs/Lawson/LGSM.mm10.orderedN.vcf
-resource:cheverud,known=false,training=true,truth=true,
    prior=12.0,/group/palmer-lab/AIL/knownSNPs/LGSM.mm10.flipremoved.vcf
-an SB -an MQRankSum -an ReadPosRankSum -an LowMQ
-an NCC -an GQ_MEAN
#-tranche
#-allPoly (trust all input training sets contain only polymorhpic sites to speed up computation)
-mode SNP
-recalFile snpValCohort.recal
-tranchesFile snpValCohort.tranches
-rscriptFile snpValCohort.plots.R

   6035271 vcf_chr_1.vcf
   5311788 vcf_chr_2.vcf
   4951557 vcf_chr_3.vcf
   4430973 vcf_chr_4.vcf
   4669216 vcf_chr_5.vcf
   4434987 vcf_chr_6.vcf
   4212636 vcf_chr_7.vcf
   4055978 vcf_chr_8.vcf
   3823719 vcf_chr_9.vcf
   3534874 vcf_chr_11.vcf
   3597096 vcf_chr_12.vcf
   3652121 vcf_chr_13.vcf
   3644278 vcf_chr_14.vcf
   3260020 vcf_chr_15.vcf
   2941055 vcf_chr_16.vcf
   2750530 vcf_chr_17.vcf
   2810474 vcf_chr_18.vcf
   1882674 vcf_chr_19.vcf
  3267153 vcf_chr_X.vcf
     5595 vcf_chr_Y.vcf
73 271 995 TOTAL

73271995 - (20*14) = 73 271 715
77373261
# Apply Recalibration

3 272 720



## SUBSET DATA ------------------------------------------------------------------------------------

## exclude indels
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T SelectVariants
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
-selectTypeToExclude INDEL
-o snpGenosValidationCohort.vcf

## extract giga samples and snps
## gigaMice that are sibs: 46134-5 (F50), 48369-70 (F51), 51240,42 (F53), 52123-4 (F55)
## java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T SelectVariants
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
--selectTypeToInclude SNP
--keepIDs /group/palmer-lab/AIL/GBS/vars/gigaSnps.list
--sample_file /group/palmer-lab/AIL/GBS/vars/gigaMice.list
-o snpValidationGigaMice.vcf

## extract littleped samples and snps
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T SelectVariants
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
--selectTypeToInclude SNP
--keepIDs /group/palmer-lab/AIL/GBS/vars/ -----------fillin
--sample_file /group/palmer-lab/AIL/GBS/vars/ ---------fillin
-o snpValidationGigaMice.vcf





######################## other ----------------------------------------------------------
# count aligned reads
samtools view -c -q 1 $INFILE >> alignedReads_MQgrthan0
samtools view -c -q 1 $INFILE >> alignedReads_MQgrthan20

# depth of coverage stats
INFILE=`head -$PBS_ARRAYID /group/palmer-lab/AIL/GBS/bams/bams.list | tail -1`
BASE= `basename $INFILE .rg.requaled.realign.bam`
DIR= `dirname INFILE`
echo "running depth for: $INFILE"
## note: this awk line calculates coverage over regions but does not account for regions that were
## not covered with a read. Then it divides all the sum with the number of rows that were reported.
## but the positions without a read hit are not included. Basically - it gives a too high result.
samtools depth -f -Q 20 $INFILE |  awk '{sum+=$3; sumsq+=$3*$3} END { echo "$BASE"; print "MEAN = ",
sum/NR; print "SD = ",sqrt(sumsq/NR - (sum/NR)**2)}'
## the denominator needs to be the size of genome where the BAM file was generated against, rather
## than the number of bases covered at least once.
## get genome size from bam file: 2730871774
samtools view -H *bamfile* | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
## corrected awk line
samtools depth -f -Q 20 $INFILE |  awk '{sum+=$3; sumsq+=$3*$3} END { echo "$BASE"; print "MEAN = ",
sum/2730871774; print "SD = ",sqrt(sumsq/2730871774 - (sum/2730871774)**2)}'


# output is in bams/depthStats
# mean in line 2 and sd in line 3 of each outfile. to combine them:
dirlist <- list.files()
depthStats <- c()
for (file in dirlist) {
    tmp <- read.table(file, header=F, skip=1, nrows=2)[c(1,3)]
    depthStats <- rbind(depthStats, tmp)
    }

# length depthStats is even
meanDepth <- depthStats[seq(1, nrow(depthStats)-1, by=2),2]
sdDepth <- depthStats[seq(2, nrow(depthStats), by=2),2]
depthStats <- data.frame(meanDepth=meanDepth, sdDepth=sdDepth)

# histogram
pdf("meanDepth_sampsSeqdOnce.pdf")
hist(depthStats$meanDepth, breaks=seq(0,0.34, 0.02), plot=TRUE,
     main="Mean per-base read depth for single runs\n (Genome-wide)")
dev.off()












# could also use bedtools: depth
module load bedtools
bedtools genomecov -ibam $INFILE -g $REF_GENOME -bga -max 30



