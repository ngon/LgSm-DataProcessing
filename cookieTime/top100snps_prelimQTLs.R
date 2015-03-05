
act2 <- read.table("./cmteMeeting/act2.t.best100.txt", sep="\t", header=T)
act3 <- read.table("./cmteMeeting/act3.t.best100.txt", sep="\t", header=T)
act8 <- read.table("./cmteMeeting/act8.t.best100.txt", sep="\t", header=T)
ppi <- read.table("./cmteMeeting/avg.ppi.best100.txt", sep="\t", header=T)
wild <- read.table("./cmteMeeting/wild.best100.txt", sep="\t", header=T)

act2[act2$chr == 3,]
# top snp on chr 4 is in Ubxn10
# methylphenidate listed as an interacting chemical
# expressed in testes and neonate cerebellum/cortex
## theres a phenotype listed by MGI in the range of the first and last markers
## listed for chr3.. the phenotype is abnormal locomotor activation
## nearest gene is Hnf4g

act3
# in the region of chr 15 there are a few genes
# one is ADCY8
# reacts with cocaine, morphine, ethanol, ketamine
# Efr3a
# Oc90
# Kcnq3
# Asap1/Shag1

act8
# all top hits in chr 4
# there's a SNP in this region that's been associted with abnormal locomotor coordination, abnormal gait, and seizures.(Slc9a1, oprd1)

ppi
## all hits on chr 10
# behavior associated: abnormal gait, abnormal motor capabilities/coordination

wildness
# hits on chr1 and 8
#chr 1 - in gene called Aff3 - locomotor coordination comes up as does abnormal learning/memory/conditioning(spatial learning), rediced LTP


