lgindel<- read.delim("~/post_doc/AIL/LG.Indels", header=T)
lgindel$deletion <- grepl("-", lgindel$Variant)

lgindel <- lgindel[,-3] #get rid of ref alelle

library("BSgenome.Mmusculus.UCSC.mm10")

get_base <- function(indel){getSeq(Mmusculus, indel$Chr, start=indel$Position, width=1,as.character=T)}

lgindel$Reference <- get_base(lgindel)

#split variants by /
temp <- strsplit(as.character(lgindel$Variant), "/")
lgindel[,6:7] <- matrix(unlist(temp), ncol=2, byrow=TRUE)


#strip + or - sign
lgindel[,8]<-gsub(pattern="\\+|\\-",replacement= "",x=lgindel$V6, ignore.case=F)
lgindel[,9]<-gsub(pattern="\\+|\\-",replacement= "",x=lgindel$V7, ignore.case=F)

lgindel <- lgindel[,c(-3,-6,-7)]

#split the matrix to capture those indels with different alleles
temp1 <- lgindel[,-5]
temp2 <- lgindel[,-6]
colnames(temp1)[5] <- c("X")
colnames(temp2)[5] <- c("X")
lgindel <- rbind(temp1, temp2)

#remove duplicate values
lgindel <- unique(lgindel)

#remove reference alleles
lgindel <- lgindel[lgindel$X != "*",]

#concatenate reference and variant
lgindel$X2 <- paste(lgindel$Reference, lgindel$X, sep="")

#swap reference and variant for deletions
lgindel  <- transform(lgindel, Reference = ifelse(deletion == T, X2, Reference), X2 = ifelse(deletion == T, Reference, X2))

#clean up 
lgindel <- lgindel[,c(-3,-5)]
colnames(lgindel) <- c("#CHROM", "POS", "REF", "ALT")
lgindel$ID <- rep("*",nrow(lgindel))
lgindel$QUAL <- rep("100", nrow(lgindel))
lgindel$FILTER <- rep("PASS", nrow(lgindel))
lgindel$INFO <- rep("*", nrow(lgindel))
lgindel$FORMAT <- rep("GT",nrow(lgindel))
lgindel$LG <- rep("1/1", nrow(lgindel))
lgindel$SM <- rep("0/0", nrow(lgindel))
lgindel <- lgindel[,c(1,2,5,3,4,6,7,8,9,10,11)]