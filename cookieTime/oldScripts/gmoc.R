
########## CPP DATA F50 - F55 ###########
setwd("C:/Users/Administrator/Desktop/Scripts/AIL CPP interim analyses/gmoc2014")
library(plyr); 
library(reshape); 
library(ggplot2)
load("")

data = read.table(file="cppGluData.csv", sep=",", header=TRUE)
data = data[-(877:991),]
data = data[,-(3:11)]
data = data[-491,]
# variables: cpp8.p.corrected, cpp1.p.corrected, cpp.diff.corrected, 
# cc, glu.weight, glucose, act1.t, cpp8.t, cpp1.t
cppdata = data


ids = read.table(file="ids.csv", sep=",", header=FALSE)

gbsdata = subset(cppdata, cppdata$id %in% ids$V1)

#########################################

## recode variables
#cppdata$id = as.character(cppdata$id)
#cppdata$cpp.box = as.factor(cppdata$cpp.box)
#cppdata$cage = as.factor(cppdata$cage)
#cppdata$sex = as.factor(cppdata$sex)
#cppdata$sex = droplevels(cppdata$sex)
#cppdata$gen = as.factor(cppdata$gen)
#cppdata$d5.weight = as.numeric(cppdata$d5.weight)
#cppdata$d8.weight = as.numeric(cppdata$d8.weight)

## remove empty rows
#cppdata = cppdata[-819,];
#cppdata = cppdata[-818,]

## create new variables
#cppdata$sen = cppdata$act4.t - cppdata$act2.t
#cppdata$sc.chg = cppdata$sc8.t - cppdata$sc1.t

## data dimensions
#dim(cppdata); names(cppdata) # 817 x 96 df

## save file
#write.table(cppdata, file="cpp54c.csv", sep=",", col.names=T, row.names=F)




## exploratory analysis
## t-tests, both sexes, totals
Tcpp = t.test(cpp8.t, cpp1.t, data=cppdata, paired=TRUE) # *  7.207826e-54
Tactcpp = t.test(act8.t, act1.t, data=cppdata, paired=TRUE) # 0.008825248
## t-tests, both sexes, first 5 min
Tcpp5 = t.test(cpp8.1, cpp1.1, data=cppdata, paired=TRUE) # * 4.204007e-92
Tactcpp5 = t.test(act8.1, act1.1, data=cppdata, paired=TRUE) #  1.33736e-39

attach(gbs)
## t-tests, both sexes, totals for GBS data set
Tcpp2 = t.test(cpp8.t, cpp1.t, data=gbs, paired=TRUE) # * 4.159446e-21 
Tactcpp2 = t.test(act8.t, act1.t, data=gbs, paired=TRUE) # 0.8173993
## t-tests, both sexes, first 5 min for GBS data set
Tcpp3 = t.test(cpp8.1, cpp1.1, data=gbs, paired=TRUE) # * 2.344559e-39
Tactcpp3 = t.test(act8.1, act1.1, data=gbs, paired=TRUE) #  6.282493e-16


### Boxplots
attach(cppdata)
methcpp = cbind(id, cpp1.t, cpp8.t)
methcpp = as.data.frame(methcpp)
methcpp = rename(methcpp, c(cpp1.t = "Day 1", cpp8.t = "Day 8"))
methcpp = melt.data.frame(methcpp, id.vars = c(1), measure.vars=c(2,3))
methcpp = rename(methcpp, c(variable= "Trial", value="Time"))
methcpp = na.omit(methcpp)
methcpp$Time = as.numeric(as.character(methcpp$Time))
methcpp.sum = summarySE(methcpp, measurevar="Time", groupvars=c("Trial"))

dodge = position_dodge(width=0.9)



bp = ggplot(data=methcpp, aes(x=Trial, y=Time)) +
    geom_boxplot(position="dodge", notch=TRUE, outlier.size=2, outlier.shape=16, width=.4) +
    ylab("Seconds spent on Meth-paired side") +
    ggtitle("CPP for 1 mg/kg Meth \n in F50-55 AIL")+
  scale_y_continuous(breaks=c(0,250,500,750,1000,1250,1500))+
  
   stat_summary(data = methcpp, aes(x=Trial, y=Time), 
                 fun.y=mean, colour="darkgreen", geom="point", 
                 shape=16, size=3,show_guide = FALSE) +
  
  geom_linerange(data=methcpp.sum, aes(x=Trial, y=se),
                 ymax=(methcpp.sum$Time + methcpp.sum$se), ymin=(methcpp.sum$Time-methcpp.sum$se)
                 , linetype=1, size=1, colour="green")+
  
  geom_hline(yintercept =900, linetype=2)+
  
  theme(plot.title= element_text(size=18), 
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position= "none")


png(filename="./cpp.alldata.png", 
    height=420, width=280)
bp
dev.off()


### cpp data in only gbs data

gbs2 = cbind(id, cpp1.t, cpp8.t)
gbs2 = as.data.frame(gbs2)
gbs2 = rename(gbs2, c(cpp1.t = "Day 1", cpp8.t = "Day 8"))
gbs2 = melt.data.frame(gbs2, id.vars = c(1), measure.vars=c(2,3))
gbs2 = rename(gbs2, c(variable= "Trial", value="Time"))
gbs2 = na.omit(gbs2)
gbs2$Time = as.numeric(as.character(gbs2$Time))
gbs=gbs2

gbs.sum = summarySE(gbs, measurevar="Time", groupvars=c("Trial"))


bp2 = ggplot(data=gbs, aes(x=Trial, y=Time)) +
  geom_boxplot(position="dodge", notch=TRUE, outlier.size=2, outlier.shape=16, width=.4) +
  ylab("Seconds spent on Meth-paired side") +
  ggtitle("CPP for 1 mg/kg Meth \n in F50-51 AIL")+
  scale_y_continuous(breaks=c(0,250,500,750,1000,1250, 1500,1750))+
  
  stat_summary(data = gbs, aes(x=Trial, y=Time), 
               fun.y=mean, colour="darkgreen", geom="point", 
               shape=16, size=3,show_guide = FALSE) +
  
  geom_linerange(data=gbs.sum, aes(x=Trial, y=se),
                 ymax=(gbs.sum$Time + gbs.sum$se), ymin=(gbs.sum$Time-gbs.sum$se)
                 , linetype=1, size=1, colour="green")+
  
   geom_hline(yintercept =900, linetype=2)+
  
  theme(plot.title= element_text(size=18), 
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position= "none")




png(filename="./cpp.gbsdata.png", 
    height=420, width=280)
bp2
dev.off()



####### glucose gbs data
### cpp data in only gbs data

gbs2 = cbind(id, cpp1.t, cpp8.t)
gbs2 = as.data.frame(gbs2)
gbs2 = rename(gbs2, c(cpp1.t = "Day 1", cpp8.t = "Day 8"))
gbs2 = melt.data.frame(gbs2, id.vars = c(1), measure.vars=c(2,3))
gbs2 = rename(gbs2, c(variable= "Trial", value="Time"))
gbs2 = na.omit(gbs2)
gbs2$Time = as.numeric(as.character(gbs2$Time))
gbs=gbs2

indivs = rep("F50-51", 387)
gbsdata$indiv = indivs

gbs.sum = summarySE(gbs, measurevar="Time", groupvars=c("Trial"))
sum <- c(99.35401, 28.71842, 1.459839)
names(sum) <- c("mean", "sd", "se")

bp2 = ggplot(data=gbsdata, aes(y=glucose, x=gbsdata$indiv))+
  geom_boxplot(position=dodge, notch=TRUE, outlier.size=2, outlier.shape=16, width=.3) +
  ylab("Glucose (mg/dL)") +
  xlab(" ")+
  ggtitle("Fasting blood glucose \n levels (4h)") +
  
  scale_y_continuous(breaks=c(0,25,50,75,100,125,150,175,200))+
    
  stat_summary(fun.y=mean, colour="blue", geom="point", 
               shape=16, size=3,show_guide = FALSE) +
  

  theme(plot.title= element_text(size=18), 
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position= "none")

png(filename="./glu.gbsdata.png", 
    height=420, width=280)
bp2
dev.off()
