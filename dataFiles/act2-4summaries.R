source("./oldScripts/SummarySE.R")
library(plyr)
library(ggplot2)

# read in correct data frame
load(cppData, file="./CPPandPPI.11sep14.RData")
data <- cppData

# create activity df for exploratory analysis
actWcovs <- cbind(data$id, data$gen, data$cage, data$sex, data$cpp.age,
             data$cpp.box, data[c(37:43, 59:93)])
names(actWcovs)[1:6] <- c("id", "gen", "cage", "sex", "cpp.age", "cpp.box") 

# paired t tests for activity on meth and saline
senst <- t.test(actWcovs$act2.t, actWcovs$act4.t, paired=TRUE)
senst$p.value
salLoct <- t.test(actWcovs$act3.t, actWcovs$act5.t, paired=TRUE)
salLoct$p.value

# create sens variable
actWcovs$sens <- actWcovs$act4.t - actWcovs$act2.t

# create male and female activity dfs
actF <- actWcovs[actWcovs$sex == "F",]
actM <- actWcovs[actWcovs$sex == "M",]

# t tests on male and female activity (all are p 0.02(sz) or lower)
d1at <- t.test(actF$act1.t, actM$act1.t, paired=F)
d1at$p.value
d2at <- t.test(actF$act2.t, actM$act2.t, paired=F)
d2at$p.value
d3at <- t.test(actF$act3.t, actM$act3.t, paired=F)
d3at$p.value
d4at <- t.test(actF$act4.t, actM$act4.t, paired=F)
d4at$p.value
d5at <- t.test(actF$act5.t, actM$act5.t, paired=F)
d5at$p.value
d8at <- t.test(actF$act8.t, actM$act8.t, paired=F)
d8at$p.value
senMFt <- t.test(actF$sens, actM$sens, paired=F)
senMFt$p.value


########### COVARIATE ANALYSIS ###################
act <- actWcovs 

sensGen <- lm(act$sens ~ as.factor(act$gen)); summary(sensGen) # gen52, p=0.0303
sensSex <- lm(act$sens ~ act$sex); summary(sensSex) # sex, p=0.0167
sensBox <- lm(act$sens ~ as.factor(act$cpp.box)); summary(sensBox) # box7, p=0.006

d1Gen <- lm(act$act1.t ~ as.factor(act$gen)); summary(d1Gen) # no effect
d1Box <- lm(act$act1.t ~ as.factor(act$cpp.box)); summary(d1Box) # effects of box 3-9,11,12

d2Box <- lm(act$act2.t ~ as.factor(act$cpp.box)); summary(d2Box) # box 5,7,8,10,12
d2Gen <- lm(act$act2.t ~ as.factor(act$gen)); summary(d2Gen) # gen 52, 56

d3Box <- lm(act$act3.t ~ as.factor(act$cpp.box)); summary(d3Box) # box 7,8,10,11,12
d3Gen <- lm(act$act3.t ~ as.factor(act$gen)); summary(d3Gen) # all gen except 52

d4Box <- lm(act$act4.t ~ as.factor(act$cpp.box)); summary(d4Box) # box 7,8,11
d4Gen <- lm(act$act4.t ~ as.factor(act$gen)); summary(d4Gen) # no effect

d5Box <- lm(act$act5.t ~ as.factor(act$cpp.box)); summary(d5Box) # box 7,8
d5Gen <- lm(act$act5.t ~ as.factor(act$gen)); summary(d5Gen) # gen 51, 53-56

d8Box <- lm(act$act8.t ~ as.factor(act$cpp.box)); summary(d8Box) # box 3,4,6-9,11,12
d8Gen <- lm(act$act8.t ~ as.factor(act$gen)); summary(d8Gen) # gen 56


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
        geom_line(lwd=1) +
        geom_point() +
        ggtitle("Day 4 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(300,550))+
        theme_bw()+
        theme(axis.title.x = element_text(colour="black", size=14),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=14),
              plot.title = element_text(colour="black", size=16),
              legend.position = c(0.2, 0.25))
   


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
        geom_line(lwd=1) +
        geom_point() +
        ggtitle("Day 2 Activity")+
        xlab("Time (min)") +
        ylab("Distance (cm)") +
        scale_x_continuous(breaks=c(5,10,15,20,25,30))+
        scale_y_continuous(limits=c(300,500))+
        theme_bw()+
        theme(axis.title.x = element_text(colour="black", size=14),
              axis.text = element_text(colour="black", size=12),
              axis.title.y = element_text(colour="black", size=14),
              plot.title = element_text(colour="black", size=16),
              legend.position = c(0.17, 0.17))

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
         
        








