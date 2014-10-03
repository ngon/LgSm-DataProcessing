
cpp <- read.table(file="jul2014.csv", sep=",", header=TRUE, na.strings="NA")
act <- cbind(data[1:2], data[9:11], data[30:35], data[52:86])

# remove NA entries from activity day 3 data (n=5)
three <- cbind(act[1], act[3], act[19:25], data[17])
three <- three[complete.cases(three),]
## remove cpp data for which there are comerrs on day 8 or 1 (n=32)
cpp <- cbind(data[1], data[4],data[6], data[9],data[16:17], 
             data[37:51], data[87:103], data[108])
cpp.full <- cpp[cpp$comerr8 == 0,]
cpp.t <- cpp.full[cpp.full$comerr1 == 0,]

library(plyr); 
library(reshape); 
library(ggplot2)

# source("colmean.R")
# colMeans <- colmean(data)
# colErrs <- colse(data)
# 
# act1bin <- colMeans[30:35]
# act8bin <- colMeans[80:85]
# act3bin <- colMeans[59:64]
# act5bin <- colMeans[73:78]

# ########### Day 1 and Day 8 activity by bin
# key <- c("Day 1", "Day 8")
# colors <- c("red", "blue")
# png(filename="./Day1vsDay8-ActivityLinePlot.png", height = 400, width=500)
# par(oma = c(.5,.5,0.5,0.5)+0.2, mar=c(3,3,1.5,1)+1.2)
# yy = seq(from=1000, to=3000, by=300)
# plot(act1bin, type = "o", col=colors[1], axes=F, ann=F, lwd=2, lty=1, ylim=range(yy))
# axis(1, at=1:6, lab=c("5", "10", "15", "20", "25", "30"), cex.lab=1.2)
# axis(2, c(800, 1200, 1600, 2000, 2400, 2800), cex.lab=2)
# box()
# lines(act8bin, type="o", lty=1, col=colors[2], lwd=2)
# legend("topright", key, cex=1.5, fill=colors[1:2], 
#        border=colors[1:2], col=colors[1:2], lwd=2, lty=1)
# title(ylab="Distance traveled (cm)", xlab="Time (min)", cex.lab=1.6,
#       main="Activity Day 1 and Day 8")    
# dev.off()

####################
source("SummarySE.R")

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

test1= rename(test1, c("act4.1"="Distance"))
test2=rename(test2, c("act4.2"="Distance"))
test3=rename(test3, c("act4.3"="Distance"))
test4=rename(test4, c("act4.4"="Distance"))
test5=rename(test5, c("act4.5"="Distance"))
test6=rename(test6, c("act4.6"="Distance"))

new= rbind(test1, test2, test3, test4, test5, test6)
new$sex= as.factor(new$sex)

ailD4act= ggplot(data=new, aes(x=Time, y=Distance, colour=sex))+
  geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
  geom_line(lwd=1) +
  geom_point() +
  ggtitle("Day 4 Activity")+
  xlab("Time (min)") +
  ylab("Distance (cm)") +
  scale_x_continuous(breaks=c(5,10,15,20,25,30))+
  scale_y_continuous(limits=c(300,500))+
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20))

############# cpp difference boxplot #################3
boxplot(data$cpp.diff, 
        main="Change in preference for Meth \n paired side (Day8 - Day1)", 
        ylab = "Percent difference", 
        xlab="F50-55")
abline(h=0, col="red", lty=3)
abline(h=0.055, col="blue", lty=3)

##################### sensitization data ###############
bin1 = act[26] - act[12]
bin2 = act[27] - act[13]
bin3 = act[28] - act[14]
bin4 = act[29] - act[15]
bin5 = act[30] - act[16]
bin6 = act[31] - act[17]
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

new= rbind(test1, test2, test3, test4, test5, test6)
new$sex= as.factor(new$sex)

############### sensitization ##############

sens= ggplot(data=new, aes(x=Time, y=Distance, colour=sex))+
  geom_errorbar(aes(ymin=Distance-se, ymax=Distance+se), width=.1) +
  geom_line(lwd=1) +
  geom_point() +
  ggtitle("Locomotor sensitization \n Day4 - Day2 activity")+
  xlab("Time (min)") +
  ylab("Distance (cm)") +
  scale_x_continuous(breaks=c(5,10,15,20,25,30))+
  scale_y_continuous(limits=c(10,100))+
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20))


############# plots of activity and cpp by chamber #########
act1.t= 
  ggplot(data=data, aes(x=cpp.box,y=act1.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act2.t= 
  ggplot(data=data, aes(x=cpp.box,y=act2.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 2 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3.t= 
  ggplot(data=data, aes(x=cpp.box,y=act3.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act4.t= 
  ggplot(data=data, aes(x=cpp.box,y=act4.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 4 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act5.t= 
  ggplot(data=data, aes(x=cpp.box,y=act5.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 5 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act8.t= 
  ggplot(data=data, aes(x=cpp.box,y=act8.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity by Chamber")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

cpp1.t= 
  ggplot(data=data, aes(x=cpp.box,y=cpp1.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 1 Preference")+
  xlab("Locomotor chamber") +
  ylab("Sec spent on Meth-paired side") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

cpp8.t= 
  ggplot(data=data, aes(x=cpp.box,y=cpp8.t, color=as.factor(cpp.box)))+
  geom_boxplot() +
  ggtitle("Day 8 CPP")+
  xlab("Locomotor chamber") +
  ylab("Sec spent on Meth-paired side") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


########### day 1 and 8 activity by chamber stratified by sex ############
females = activity[activity$sex == "F",]
males = activity[activity$sex == "M",]


femAct1= 
  ggplot(data=females, aes(x=box,y=act1, color=as.factor(box)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity in Females")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

malAct1= 
  ggplot(data=males, aes(x=box,y=act1, color=as.factor(box)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity in Males")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

femAct8= 
  ggplot(data=females, aes(x=box,y=act8, color=as.factor(box)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity in Females")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

malAct8= 
  ggplot(data=males, aes(x=box,y=act8, color=as.factor(box)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity in Males")+
  xlab("Locomotor chamber") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


############# plots of activity and cpp by batch #########
## jan - 119, feb - 58, mar - 111, apr - 119, may - 61, jun - 94
## jul - 53, aug - 86, sep - 119, oct - 64, nov - 107, dec - 32

act1batch= 
  ggplot(data=data, aes(x=batch,y=act1.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act2batch= 
  ggplot(data=data, aes(x=batch,y=act2.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 2 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3batch= 
  ggplot(data=data, aes(x=batch,y=act3.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3batch= 
  ggplot(data=data, aes(x=batch,y=act3.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


act5batch= 
  ggplot(data=data, aes(x=batch,y=act5.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 5 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act4batch= 
  ggplot(data=data, aes(x=batch,y=act4.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 4 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act8batch= 
  ggplot(data=data, aes(x=batch,y=act8.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity")+
  xlab("Batch") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


cpp1batch= 
  ggplot(data=data, aes(x=batch,y=cpp1.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 1 Preference")+
  xlab("Batch") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

cpp8batch= 
  ggplot(data=data, aes(x=batch,y=cpp8.t, color=as.factor(batch)))+
  geom_boxplot() +
  ggtitle("Day 8  CPP reference")+
  xlab("Batch") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

############# plots of activity and cpp by month #########
act1Month= 
  ggplot(data=data, aes(x=Month,y=act1.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act2Month= 
  ggplot(data=data, aes(x=Month,y=act2.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 2 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3Month= 
  ggplot(data=data, aes(x=Month,y=act3.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3Month= 
  ggplot(data=data, aes(x=Month,y=act3.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


act5Month= 
  ggplot(data=data, aes(x=Month,y=act5.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 5 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act4Month= 
  ggplot(data=data, aes(x=Month,y=act4.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 4 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act8Month= 
  ggplot(data=data, aes(x=Month,y=act8.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity")+
  xlab("Month") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


cpp1Month= 
  ggplot(data=data, aes(x=Month,y=cpp1.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 1 Preference")+
  xlab("Month") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

cpp8Month= 
  ggplot(data=data, aes(x=Month,y=cpp8.t, color=as.factor(Month)))+
  geom_boxplot() +
  ggtitle("Day 8  CPP reference")+
  xlab("Month") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


################ hists of comerrs #################3

data$comerr2 = as.numeric(data$comerr2)
data$comerr4 = as.numeric(data$comerr4)
comerrTot = data$comerr1 + data$comerr2 + data$comerr3 
            + data$comerr4 + data$comerr5 + data$comerr6 + data$comerr8
data$comerrTot = comerrTot
a <- data[data$comerrTot > 1.0,]

par(mfrow = c(2,2))
hist(as.numeric(a$cpp.box), col="red", main="Comp errors by chamber", xlab="Locomotor chamber")
hist(as.numeric(a$Month), col="orangered", main="Comp errors by month",xlab="Locomotor chamber" )
hist(as.numeric(a$batch), col="orange", main="Comp errors by batch",xlab="Locomotor chamber" )
hist(as.numeric(a$gen), col="yellow", main="Comp errors by generation", xlab="Locomotor chamber")


############# by generation ###########3


############# plots of activity and cpp by gen #########
act1gen= 
  ggplot(data=data, aes(x=gen,y=act1.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 1 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act2gen= 
  ggplot(data=data, aes(x=gen,y=act2.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 2 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3gen= 
  ggplot(data=data, aes(x=gen,y=act3.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act3gen= 
  ggplot(data=data, aes(x=gen,y=act3.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 3 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


act5gen= 
  ggplot(data=data, aes(x=gen,y=act5.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 5 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act4gen= 
  ggplot(data=data, aes(x=gen,y=act4.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 4 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

act8gen= 
  ggplot(data=data, aes(x=gen,y=act8.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 8 Activity")+
  xlab("gen") +
  ylab("Distance traveled (cm)") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


cpp1gen= 
  ggplot(data=data, aes(x=gen,y=cpp1.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 1 Preference")+
  xlab("gen") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")

cpp8gen= 
  ggplot(data=data, aes(x=gen,y=cpp8.t, color=as.factor(gen)))+
  geom_boxplot() +
  ggtitle("Day 8  CPP reference")+
  xlab("gen") +
  ylab("Seconds") +
  theme(axis.title.x = element_text(colour="black", size=18),
        axis.text = element_text(colour="black", size=12),
        axis.title.y = element_text(colour="black", size=18),
        plot.title = element_text(colour="black", size=20),
        legend.position = "none")


