### all of the ppi entries for which the sum.nostim trials is greater than ~120
### were measured from box 3. 

# read complete data set
data <-read.table(file="./dataFiles/allData.txt", sep="\t", 
                  header=T, na.strings="NA")
library(plyr)
library(ggplot2)

## break ppi.weight into three classes (15-25g, 25-35g, 35-45g)
breaks <- seq(15, 50, 10)
weightbins <- cut(data$ppi.weight, breaks)
data$weightbins <- as.numeric(weightbins)


## ppi by weight class
ppi <- cbind(data$avg.ppi, data$weightbins)
ppi <- ppi[complete.cases(ppi),]
ppi <- as.data.frame(ppi)
names(ppi) <- c("avg.ppi", "weightbins")

ppibyWeight=  ggplot(data=ppi, aes(x=as.factor(weightbins),y=avg.ppi))+
        geom_boxplot() +
        ggtitle("Average PPI")+
        xlab("Weight class") +
        ylab("% PPI") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")

## ppi by box
ppib <- cbind(data$avg.ppi, data$ppi.box)
ppib <- ppib[complete.cases(ppib),]
ppib <- as.data.frame(ppib)
names(ppib) <- c("avg.ppi", "box")
ppib$box <-as.factor(ppib$box)
b5<-ppib[ppib$box == "5",]
ppib <- ppib[-43,]

ppibybox=  ggplot(data=ppib, aes(x=as.factor(box),y=avg.ppi))+
        geom_boxplot() +
        ggtitle("Average PPI")+
        xlab("Testing chamber") +
        ylab("% PPI") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")


### startle by weight class
startle <- cbind(data$startle, data$weightbins)
startle <- startle[complete.cases(startle),]
startle <- as.data.frame(startle)
names(startle) <- c("startle", "weightbins")

strbyWeight=  ggplot(data=startle, aes(x=as.factor(weightbins),y=startle))+
        geom_boxplot() +
        ggtitle("Startle")+
        xlab("Weight class") +
        ylab("Startle response") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")


gluc <- cbind(data$glucose, data$weightbins)
gluc <- gluc[complete.cases(gluc),]
gluc <- as.data.frame(gluc)
names(gluc) <- c("gluc", "weightbins")
gluc$weightbins <- as.factor(gluc$weightbins)

levels(gluc$weightbins)[levels(gluc$weightbins)=="1"] <- "15-25g"
levels(gluc$weightbins)[levels(gluc$weightbins)=="2"] <- "25-35g"
levels(gluc$weightbins)[levels(gluc$weightbins)=="3"] <- "35-45g"

gluWeight=  ggplot(data=gluc, aes(x=as.factor(weightbins),y=gluc))+
        geom_boxplot() +
        xlab("Weight class") +
        ylab("Blood glucose levels (mg/dL)") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")

glucose=  ggplot(data=data, aes(x=as.factor(sex),y=glucose))+
        geom_boxplot() +
        xlab("Sex") +
        geom_text(data=NULL, x=0.65, y=275,label="LG/J", col="red")+
        geom_text(data=NULL, x=0.65, y=260,label="SM/J", col="dodgerblue")+
        scale_y_continuous(breaks=c(50,100,150,200,250), limits=c(0,300))+
        geom_pointrange(x=0.75, y=156, ymax=197.5, ymin=114.5, col="red", linetype=2)+ #lg 
        geom_pointrange(x=1.25, y=177, ymax=221.7, ymin=132.3, col="dodgerblue", linetype=2)+ #sm f
        geom_pointrange(x=1.75, y=196, ymax=234.8, ymin=157.2, col="red", linetype=2)+ #lg m
        geom_pointrange(x=2.25, y=231, ymax=273.3, ymin=188.7, col="dodgerblue", linetype=2)+ #sm m
        ylab("Blood glucose levels (mg/dL)") +
       
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")


## tail length
data$tail <- as.numeric(as.character(data$tail))
data[data$tail < 5,]

tails=  ggplot(data=data, aes(x=as.factor(sex),y=tail))+
        geom_boxplot() +
        xlab("Sex") +
        ylab("Tail length (cm)") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")

## wildness
hist(data$wild, ylim=c(0,20), main="Wildness", 
     xlab="Number of escapes during CPP",
     ylab="Number of mice out of 939")



wild=  ggplot(data=data, aes(x=as.factor(sex),y=tail))+
        geom_boxplot() +
        xlab("Sex") +
        ylab("Tail length (cm)") +
        theme(axis.title.x = element_text(colour="black", size=18),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=18),
              plot.title = element_text(colour="black", size=20),
              legend.position = "none")


