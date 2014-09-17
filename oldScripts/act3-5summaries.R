source("./oldScripts/SummarySE.R")
library(plyr)
library(ggplot2)

# read in correct data frame
load(cppData, file="./CPPandPPI.11sep14.RData")
data <- cppData

# create activity df for exploratory analysis
act <- cbind(data$id, data$gen, data$cage, data$sex, data$cpp.age,
                  data$cpp.box, data[c(37:43, 59:93)])
names(act)[1:6] <- c("id", "gen", "cage", "sex", "cpp.age", "cpp.box") 



####################   ACTIVITY*SEX PLOTS   ###################
####################     DAY 1, 3, 5, 8     ###################

act <- act[complete.cases(act),] #1080 complete cases

test1= summarySE(act, measurevar=c("act1.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act1.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act1.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act1.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act1.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act1.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act1.1"="Distance"))
test2=rename(test2, c("act1.2"="Distance"))
test3=rename(test3, c("act1.3"="Distance"))
test4=rename(test4, c("act1.4"="Distance"))
test5=rename(test5, c("act1.5"="Distance"))
test6=rename(test6, c("act1.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD1act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
        #ggtitle("Day 1 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(1000,2500))+
        theme_tufte+
        theme(
              legend.position = "none")



test1= summarySE(act, measurevar=c("act8.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act8.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act8.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act8.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act8.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act8.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act8.1"="Distance"))
test2=rename(test2, c("act8.2"="Distance"))
test3=rename(test3, c("act8.3"="Distance"))
test4=rename(test4, c("act8.4"="Distance"))
test5=rename(test5, c("act8.5"="Distance"))
test6=rename(test6, c("act8.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD8act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
        #ggtitle("Day 8 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(1000,2500))+
        theme_tufte+
        theme(
              legend.position = c(0.9,0.9))
             


test1= summarySE(act, measurevar=c("act3.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act3.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act3.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act3.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act3.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act3.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act3.1"="Distance"))
test2=rename(test2, c("act3.2"="Distance"))
test3=rename(test3, c("act3.3"="Distance"))
test4=rename(test4, c("act3.4"="Distance"))
test5=rename(test5, c("act3.5"="Distance"))
test6=rename(test6, c("act3.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD3act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
        #ggtitle("Day 3 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(100,450))+
        theme_tufte+
        theme(legend.position = "none")
            




test1= summarySE(act, measurevar=c("act5.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act5.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act5.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act5.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act5.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act5.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act5.1"="Distance"))
test2=rename(test2, c("act5.2"="Distance"))
test3=rename(test3, c("act5.3"="Distance"))
test4=rename(test4, c("act5.4"="Distance"))
test5=rename(test5, c("act5.5"="Distance"))
test6=rename(test6, c("act5.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD5act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
        #ggtitle("Day 5 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(100,450))+
        theme_tufte+
        theme(legend.position = "none")