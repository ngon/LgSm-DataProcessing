## Figures for thesis cmte mtg 15 Oct 2015

# VennDiagram package
install.packages('VennDiagram')
library(VennDiagram)

# Load RNAID table and subset it
hiprna<-rnaid[which(rnaid$tissue == "hip"),] # N = 280
strrna<-rnaid[which(rnaid$tissue == "str"),] # N = 254
pfcrna<-rnaid[which(rnaid$tissue == "pfc"),] # N = 255

# Find the size of intersection between data sets
length(intersect(hiprna$id, strrna$id)) # 206
length(intersect(pfcrna$id, strrna$id)) # 182
length(intersect(hiprna$id, pfcrna$id)) # 201
pfcStr <- intersect(pfcrna$id, strrna$id)
length(intersect(pfcStr, hiprna$id)) # 155

# Draw venn diagram
draw.triple.venn( area1=280, area2=254, area3=255, n12=206, n23=182, n13=201, n123=155,
                 category=c("HIP", "STR", "PFC"), lty=c("blank", "blank", "blank"), lwd=c(0,0,0),
                 col=c("dodgerblue", "green", "yellow"), fill=c("dodgerblue", "green", "yellow"),
                 alpha=rep(0.7, 3), cex=rep(0.9, 7), fontfamily=rep("sans", 7),
                 cat.fontfamily=rep("sans", 3), cat.cex=rep(0.9,3))




