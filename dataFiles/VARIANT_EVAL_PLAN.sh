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
module load samtools
samtools view -H SM_rmdup_sorted_hl.bam > SM_header
samtools view -H LG_rmdup_sorted_hl.bam > LG_header
# vim: %s/SN:/SN:chr/g
samtools reheader SM_header SM_rmdup_sorted_hl.bam > SM.rehead.bam
samtools reheader LG_header LG_rmdup_sorted_hl.bam > LG.rehead.bam #####

# create list of unique reads in Lawson bams.
module load samtools
samtools view SM_rmdup_sorted_hl.bam | cut -f1 | sort -u > sm.readNames
samtools view LG_rmdup_sorted_hl.bam | cut -f1 | sort -u > lg.readNames

# separate by read group (first 4 fields of readName)

# filter reads
java -Xmx2g -jar /group/palmer-lab/tools/picard-tools-1.126/picard.jar FilterSamReads
I=SM.rehead.bam
FILTER=includeReadList
READ_LIST_FILE=
O=/group/palmer-lab/AIL/knownSNPs/Lawson/SM.runNNN.bam
SORT_ORDER=coordinate
ASSUME_SORTED=false
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

## Set variant IDs
module load bcftools/1.2
bcftools annotate --set-id +'ail\.%CHROM\.%POS' $INFILE


## CALL GENOTYPES ---------------------------------------------------------------------------------
java -Xmx6g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T GenotypeGVCFs
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
--variant allelesValidationCohort.g.vcf
-o genotypesValidationCohort.vcf


## ANNOTATE VARIANTS ------------------------------------------------------------------------------
java -Xmx6g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T VariantAnnotator
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
-A GCContent -A HardyWeinberg
# for ped: -A PossibleDeNovo -A InbreedingCoeff
# # MAKE PED
-o annotGenotypesValidationCohort.vcf


## VQSR -------------------------------------------------------------------------------------------



## SUBSET DATA ------------------------------------------------------------------------------------

## exclude indels
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
-T SelectVariants
-R:REFSEQ /group/palmer-lab/reference_genomes/mouse/mm10.fasta
-V genotypesValidationCohort.vcf
-selectTypeToExclude INDEL
-o snpGenosValidationCohort.vcf

## extract giga samples and snps
java -Xmx2g -jar /apps/software/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar
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










