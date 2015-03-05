source("./oldScripts/SummarySE.R")
library(plyr)
library(ggplot2)

# read in correct data frame
load(cppData, file="./CPPandPPI.11sep14.RData")
data <- cppData



####################   ACTIVITY*SEX PLOTS   ###################
####################   DAY 2 - DAY 4 - SZN  ###################

act <- cbind(data$id, data$sex, data[c(59:65, 73:79)])
names(act)[1:2] <- c("id", "sex")
act <- act[complete.cases(act),] #1106 complete cases

test1= summarySE(act, measurevar=c("act4.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act4.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act4.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act4.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act4.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act4.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act4.1"="Distance"))
test2=rename(test2, c("act4.2"="Distance"))
test3=rename(test3, c("act4.3"="Distance"))
test4=rename(test4, c("act4.4"="Distance"))
test5=rename(test5, c("act4.5"="Distance"))
test6=rename(test6, c("act4.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD4act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+ 
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
        #ggtitle("Day 4 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(300,550))+
        theme_tufte+
        theme(legend.position = "none")
             
   


test1= summarySE(act, measurevar=c("act2.1"), groupvars="sex")
test2= summarySE(act, measurevar=c("act2.2"), groupvars="sex")
test3= summarySE(act, measurevar=c("act2.3"), groupvars="sex")
test4= summarySE(act, measurevar=c("act2.4"), groupvars="sex")
test5= summarySE(act, measurevar=c("act2.5"), groupvars="sex")
test6= summarySE(act, measurevar=c("act2.6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1=rename(test1, c("act2.1"="Distance"))
test2=rename(test2, c("act2.2"="Distance"))
test3=rename(test3, c("act2.3"="Distance"))
test4=rename(test4, c("act2.4"="Distance"))
test5=rename(test5, c("act2.5"="Distance"))
test6=rename(test6, c("act2.6"="Distance"))

bindata= rbind(test1, test2, test3, test4, test5, test6)
bindata$sex= as.factor(bindata$sex)

ailD2act= ggplot(data=bindata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=0.6) +
        geom_point() +
       # ggtitle("Day 2 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(300,500))+
        theme_tufte+
        theme(legend.position = "none")
              

############# sensitization data ###############
bin1 = act[10] - act[3]
bin2 = act[11] - act[4]
bin3 = act[12] - act[5]
bin4 = act[13] - act[6]
bin5 = act[14] - act[7]
bin6 = act[15] - act[8]
szn = cbind(bin1,bin2,bin3,bin4,bin5,bin6, act$sex)
names(szn) <- c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "sex")

test1= summarySE(szn, measurevar=c("bin1"), groupvars="sex")
test2= summarySE(szn, measurevar=c("bin2"), groupvars="sex")
test3= summarySE(szn, measurevar=c("bin3"), groupvars="sex")
test4= summarySE(szn, measurevar=c("bin4"), groupvars="sex")
test5= summarySE(szn, measurevar=c("bin5"), groupvars="sex")
test6= summarySE(szn, measurevar=c("bin6"), groupvars="sex")
test1$Time=c(5, 5)
test2$Time=c(10, 10)
test3$Time=c(15, 15)
test4$Time=c(20, 20)
test5$Time=c(25, 25)
test6$Time=c(30, 30)
test1= rename(test1, c("bin1"="Distance"))
test2=rename(test2, c("bin2"="Distance"))
test3=rename(test3, c("bin3"="Distance"))
test4=rename(test4, c("bin4"="Distance"))
test5=rename(test5, c("bin5"="Distance"))
test6=rename(test6, c("bin6"="Distance"))

szdata= rbind(test1, test2, test3, test4, test5, test6)
szdata$sex= as.factor(szdata$sex)

sensMF= ggplot(data=szdata, aes(x=Time, y=Distance, colour=sex))+
        geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
        geom_line(lwd=1) +
        geom_point() +
        ggtitle("Locomotor sensitization \n Day4 - Day2 activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(10,100))+
        theme_bw()+
        theme(axis.title.x = element_text(colour="black", size=14),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=14),
              plot.title = element_text(colour="black", size=16),
              legend.position = c(0.2, 0.25))
         
        








