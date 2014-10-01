## violin plots of trait distributions

## get mouse IDs
mice <- read.table("./genotyped.samples.txt", sep="\t", header=F)
names(mice)[1] <- "id"

## determine which mice have been phenotyped
imputeData <- merge(x=data, y=mice, by="id", all.y=TRUE)
sum(!is.na(imputeData$act1.t)) ## 514 have phenotypes

## data frame of mice both phenotyped and genotyped
qtldata <- imputeData[!is.na(imputeData$act1.t),]

## sample size for each generation
## there are a total of 253 females and 261 males
n50 <- qtldata[qtldata$gen == 50,] # 128 
n51 <- qtldata[qtldata$gen == 51,] # 159 
n52 <- qtldata[qtldata$gen == 52,] # 156 
n53 <- qtldata[qtldata$gen == 53,] # 71 

hist(qtldata$act2.t) # good
hist(qtldata$act4.t) # better


library(reshape2)

## subsample of 514 phenotyped mice
qtlmelt = as.data.frame(cbind(qtldata$id, qtldata$act1.t, qtldata$act2.t, qtldata$act3.t, 
                qtldata$act4.t, qtldata$act5.t, qtldata$act8.t, qtldata$cpp8.t,
                qtldata$sens, qtldata$avg.ppi, qtldata$startle, qtldata$glucose, 
                qtldata$wild))
names(qtlmelt) <- c("id", "act1", "act2", "act3", "act4", "act5", "act8", "cpp", "sens",
                    "ppi", "startle", "glucose", "wild")

qtlmeltAct = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(3,5))
qtlmeltAct$value = as.numeric(as.character(qtlmeltAct$value))

qtlmeltAct2 = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(2,7))
qtlmeltAct2$value = as.numeric(as.character(qtlmeltAct2$value))

qtlmeltAct35 = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(4,6))
qtlmeltAct35$value = as.numeric(as.character(qtlmeltAct35$value))

qtlmeltppi = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(10))
qtlmeltppi$value = as.numeric(as.character(qtlmeltppi$value))

qtlmeltstr = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(11))
qtlmeltstr$value = as.numeric(as.character(qtlmeltstr$value))

qtlmeltcpp = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(8))
qtlmeltcpp$value = as.numeric(as.character(qtlmeltcpp$value))

qtlmeltsens = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(9))
qtlmeltsens$value = as.numeric(as.character(qtlmeltsens$value))

qtlmeltglu = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(12))
qtlmeltglu$value = as.numeric(as.character(qtlmeltglu$value))

qtlmeltWild = melt.data.frame(qtlmelt, id.vars = c(1), measure.vars=c(13))
qtlmeltWild$value = as.numeric(as.character(qtlmeltWild$value))


### all phenos
allmelt = as.data.frame(cbind(data$id, data$act1.t, data$act2.t, data$act3.t, 
                              data$act4.t, data$act5.t, data$act8.t, data$cpp8.t,
                              data$sens, data$avg.ppi, data$startle, data$glucose, 
                              data$wild))
names(allmelt) <- c("id", "act1", "act2", "act3", "act4", "act5", "act8", "cpp", "sens",
                    "ppi", "startle", "glucose", "wild")


allmeltAct = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(3,5))
allmeltAct$value = as.numeric(as.character(allmeltAct$value))

allmeltAct35 = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(4,6))
allmeltAct35$value = as.numeric(as.character(allmeltAct35$value))

allmeltAct2 = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(2,7))
allmeltAct2$value = as.numeric(as.character(allmeltAct2$value))

allmeltppi = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(10))
allmeltppi$value = as.numeric(as.character(allmeltppi$value))

allmeltstr = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(11))
allmeltstr$value = as.numeric(as.character(allmeltstr$value))

allmeltcpp = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(8))
allmeltcpp$value = as.numeric(as.character(allmeltcpp$value))

allmeltsens = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(9))
allmeltsens$value = as.numeric(as.character(allmeltsens$value))

allmeltglu = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(12))
allmeltglu$value = as.numeric(as.character(allmeltglu$value))

allmeltWild = melt.data.frame(allmelt, id.vars = c(1), measure.vars=c(13))
allmeltWild$value = as.numeric(as.character(allmeltWild$value))


## violin plots
library(ggplot2)

act24.violin <- 
        ggplot(data=qtlmeltAct, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltAct, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("distance (cm)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))

act35.violin <- 
        ggplot(data=qtlmeltAct35, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltAct35, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("distance (cm)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))

       
act18.violin <- 
        ggplot(data=qtlmeltAct2, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltAct2, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("distance (cm)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))
        
        

ppi.violin <- 
        ggplot(data=qtlmeltppi, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltppi, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("percent ppi") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))

str.violin <- 
        ggplot(data=qtlmeltstr, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltstr, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("startle") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))


cpp.violin <- 
        ggplot(data=qtlmeltcpp, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltcpp, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("time (sec)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))


sens.violin <- 
        ggplot(data=qtlmeltsens, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltsens, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("day 4 - day2 distance (cm)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))


glu.violin <- 
        ggplot(data=qtlmeltglu, aes(factor(variable), value)) +
        geom_violin(fill="darkslateblue", alpha=0.4) +
        geom_violin(data=allmeltglu, aes(factor(variable), value),
                    fill="white", alpha=0.3) +
        xlab(" ") +
        ylab("glucose levels (mg/dL)") +
        theme(axis.text.x=element_text(size=14),
              axis.text.y=element_text(size=12),
              axis.title.y=element_text(size=14))

# wild.violin <- 
#         ggplot(data=qtlmeltWild, aes(factor(variable), value)) +
#         geom_violin(fill="blue", alpha=0.3) +
#         geom_violin(data=allmeltWild, aes(factor(variable), value),
#                     fill="darkgreen", alpha=0.5) +
#         xlab(" ") +
#         ylab("number of escapes") +
#         theme(axis.text.x=element_text(size=14),
#               axis.text.y=element_text(size=12),
#               axis.title.y=element_text(size=14))





