library(plyr); 
library(reshape); 
library(ggplot2)
source("SummarySE.R")
source("multiplot.R")

#########1. Read and reformat data ############
###############################################

# read and format ppidata
data <- read.table(file="ppiJul2014.csv", sep=",", header=TRUE, na.strings="NA")
data$BOX <- as.factor(data$BOX)

# read and format cpp and other covariate data
cpp <- read.table(file="Jul2014.csv", sep=",", header=TRUE, na.strings="NA")
cpp$sire.age <- as.numeric(cpp$sire.age)
cpp$batch <- as.factor(cpp$batch)
# take sire.age, ppi.mo, batch and a few other traits
covs <- cbind(cpp[1], cpp[4], cpp[10], cpp[12], cpp[36], 
              cpp[58], cpp[72], cpp[93], cpp[102])
sens <- covs[7] - covs[6] # create sensitization variable
names(sens)[1] <-"sens"
covs <- cbind(covs, sens)

# read in and format info data frame
info <- read.table(file="ppi.info.csv", sep=",", header=TRUE, na.strings="NA")
info$ppi.age <- as.numeric(info$ppi.age)
info$gen <- as.factor(info$gen)
ppinfo <- info[complete.cases(info),] # removes five samples with NA entries
ppinfo <- cbind(ppinfo[1:2], ppinfo[7:10]) # keep only necessary columns

# order data by first col (id)
data <- data[do.call(order, data),]
covs <- covs[do.call(order, covs ),]
ppinfo <- ppinfo[do.call(order, ppinfo),]

# match up variables by id
ppinfo0 <- ppinfo[ppinfo$id %in% covs$id,]
covs0 <- covs[covs$id %in% ppinfo$id,]
cov.ppi <- cbind(covs0, ppinfo0)
names(cov.ppi)[11]<- "id.2"

# remove samples with errors and error column
cov.ppi <- cov.ppi[!cov.ppi$id == "54343",] 
cov.ppi <- cov.ppi[!cov.ppi$id == "54350",]
cov.ppi <- cov.ppi[c(-16)]

# cleanup
rm(covs, covs0, ppinfo, ppinfo0, cpp)

# match up data and full cov set
cov <- cov.ppi[cov.ppi$id %in% data$ID,]
d <- data[data$ID %in% cov.ppi$id,]
ppi <- cbind(cov, d)

#cleanup
rm(d, cov, cov.ppi, data)

# for plotting
data <-ppi
data <-ppi[complete.cases(ppi),]
data <-data[!data$BOX == "5",]
data <-data[!data$BOX == "NA",]
data$BOX=droplevels(data$BOX)

ppidata <- data

############# ppi by chamber #########

ppi3= 
      ggplot(data=data, aes(x=BOX,y=PPI3, color=BOX))+
      geom_boxplot() +
      ggtitle("3dB PPI")+
      xlab("PPI chamber") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

ppi6= 
      ggplot(data=data, aes(x=BOX,y=PPI6, color=BOX))+
      geom_boxplot() +
      ggtitle("6dB PPI")+
      xlab("PPI chamber") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

ppi12= 
      ggplot(data=data, aes(x=BOX,y=PPI12, color=BOX))+
      geom_boxplot() +
      ggtitle("12dB PPI")+
      xlab("PPI chamber") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

start= 
      ggplot(data=data, aes(x=BOX,y=startle, color=BOX))+
      geom_boxplot() +
      ggtitle("120 dB Startle")+
      xlab("PPI chamber") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

habit= 
      ggplot(data=data, aes(x=BOX,y=habituation, color=BOX))+
      geom_boxplot() +
      ggtitle("Habituation")+
      xlab("PPI chamber") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

############# ppi by month #########

ppi= 
      ggplot(data=data, aes(x=ppi.mo,y=avg_ppi, color=ppi.mo))+
      geom_boxplot() +
      ggtitle("Average PPI")+
      xlab("Month") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

start= 
      ggplot(data=data, aes(x=ppi.mo,y=startle, color=ppi.mo))+
      geom_boxplot() +
      ggtitle("120 dB Startle")+
      xlab("Month") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

habit= 
      ggplot(data=data, aes(x=ppi.mo,y=habituation, color=ppi.mo))+
      geom_boxplot() +
      ggtitle("Habituation")+
      xlab("Month") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

############# ppi by batch #########

ppi= 
      ggplot(data=data, aes(x=batch,y=avg_ppi, color=batch))+
      geom_boxplot() +
      ggtitle("Average PPI")+
      xlab("Batch") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

start= 
      ggplot(data=data, aes(x=batch,y=startle, color=batch))+
      geom_boxplot() +
      ggtitle("120 dB Startle")+
      xlab("Batch") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

habit= 
      ggplot(data=data, aes(x=batch,y=habituation, color=batch))+
      geom_boxplot() +
      ggtitle("Habituation")+
      xlab("Batch") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

############# ppi by gen #########

ppi= 
      ggplot(data=data, aes(x=gen,y=avg_ppi, color=gen))+
      geom_boxplot() +
      ggtitle("Average PPI")+
      xlab("Generation") +
      ylab("% PPI") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

start= 
      ggplot(data=data, aes(x=gen,y=startle, color=gen))+
      geom_boxplot() +
      ggtitle("120 dB Startle")+
      xlab("Generation") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")

habit= 
      ggplot(data=data, aes(x=gen,y=habituation, color=gen))+
      geom_boxplot() +
      ggtitle("Habituation")+
      xlab("Generation") +
      ylab("Startle amplitude") +
      theme(axis.title.x = element_text(colour="black", size=18),
            axis.text = element_text(colour="black", size=12),
            axis.title.y = element_text(colour="black", size=18),
            plot.title = element_text(colour="black", size=20),
            legend.position = "none")