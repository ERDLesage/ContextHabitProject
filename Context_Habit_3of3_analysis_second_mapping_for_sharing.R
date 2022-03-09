#########################################################################################~##
## Context habit project - criterion training analyses
## Script 3: Learning during criterion training: fitting RT and accuracy ####
## 2021 Elise Lesage (elise.lesage@ugent.be) 
## --------------------------------------------------------------------------------------~--
## Notes: 
## - Comparing learning between initial criterion training and second criterion training (group COT only)
## - RT: bootstrapping 2 parameter power function (slope parameter)
## - Accuracy: logistic regression (slope parameter)
##
rm(list = ls()) # clear wm 
##~########################################################################################~##
## SET PATHS AND SET WHETHER YOU WANT OUTPUT ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set your paths 
datadir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/writings/data_for_sharing/output"
outputdir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/writings/data_for_sharing/output"
fitfigdir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/writings/data_for_sharing/fitfigs"


# Switch writing (figures) on (TRUE) or off (FALSE) 
write_it = TRUE

library(doBy)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(afex)
library(phia)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(dplyr)
library(tidyr)
library(Rmisc)
library(lmerTest)
library(RcppRoll)
library(coin)
library(scales)
library(minpack.lm)
library(simpleboot)
library(boot)
library(arm)

# load workspace
setwd(datadir)
load("Context_habit_processed_dataframes.RData")

##~###########################################################################################~##
## Function definition for bootstrap ####
##~###########################################################################################~##

exp_fit <- function(data, indices) {
  sample <- data[indices, ]
  m <- colMeans(sample, na.rm=TRUE)
  x <- 3:55
  d<-data.frame(y=m, Trial=x)
  fo <- m ~ 1/(1 + x^c) #formula
  expon <- nls(fo, data=d, start = list(c = -0.1), nls.control(minFactor=0.0001, warnOnly = TRUE) , alg = "plinear")
  return(coef(expon))
}

##~############################################################################################~##
## MAIN SCRIPT ####
##~############################################################################################~##

## (set colors) ####
grps_cols_order1 = c("#009900", "#FF8000","#0066CC") # CT EOT COT
grps_cols_order2 = c("#FF8000","#0066CC","#009900") # EOT COT CT


## ~ figures: criterion training block length ####
crit_acc <-  summarySE(filter(crit, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "Day_Block"), na.rm=TRUE)
crit_N <- crit_acc[1:4]
crit_N$CritTrials <- crit_N$N
rm(crit_acc)


crit_table <-summarySE(crit_N, measurevar="CritTrials", groupvars =c("Group", "Day_Block"), conf.interval = .95, na.rm=TRUE)
ct1_summ <- summarySE(filter(crit_N, Day_Block == "1_1"), measurevar="CritTrials", conf.interval = .95, na.rm=TRUE)
crit_table$Group2 <- factor(crit_table$Group, levels = c("CT","EOT","COT"))
crit_N$Group2 <- factor(crit_N$Group, levels = c("CT","EOT","COT"))
crit1_bar <- ggplot(filter(crit_table, Day_Block=="1_1"), aes(x=Group2, y=CritTrials)) + scale_y_continuous(limits = c(0,130))+
  geom_bar(aes(x=Group2, y=CritTrials, color=Group2, width=0.4), fill="white", size = 2, stat="identity", position=position_dodge(.4)) +
  geom_errorbar(aes(ymin=CritTrials-se, ymax=CritTrials+se), colour="black", width=0, position=position_dodge(.4), size = 1.5) +
  geom_point(data=filter(crit_N, Day_Block=="1_1"), aes(x=Group2, y= CritTrials, color=Group2),stroke =1, alpha=.5, position=position_jitter(.15), size = 2, stat="identity") +
  xlab("Group") + ylab("# trials needed to reach criterion")+ 
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=grps_cols_order1)+
  #scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
  theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"), axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) #+
#theme(plot.title=element_text(size=28, colour="Black"))
crit1_bar

## (set colors) ####
grps_cols_order1 = c("#009900", "#FF8000","#0066CC") # CT EOT COT
grps_cols_order2 = c("#FF8000","#0066CC","#009900") # EOT COT CT


## ~ figures: criterion training block length ####
crit_acc <-  summarySE(filter(crit, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "Day_Block"), na.rm=TRUE)
crit_N <- crit_acc[1:4]
crit_N$CritTrials <- crit_N$N
rm(crit_acc)
# crit_table <-summarySE(crit_N, measurevar="CritTrials", groupvars =c("Group", "Day_Block"), conf.interval = .95, na.rm=TRUE)
# ct1_summ <- summarySE(filter(crit_N, Day_Block == "1_1"), measurevar="CritTrials", conf.interval = .95, na.rm=TRUE)
# crit1_boxplot <- ggplot(filter(crit_N, Day_Block=="1_1"), aes(x=Group2, y=CritTrials)) + scale_y_continuous(limits = c(0,140))+
#   geom_boxplot(aes(color=Group2), fill="white", width=.5, size = 1.5, position=position_dodge(.4)) +
#   geom_point(aes(x=Group2, y= CritTrials, color=Group2),stroke =1, alpha=.5, position=position_jitter(.25), size = 2, stat="identity") +
#   xlab("Group") + ylab("# trials needed to reach criterion")+ 
#   scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
#   scale_fill_manual(values=grps_cols_order1)+
#   theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
#   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
#   theme(panel.border = element_blank(), axis.line = element_line())+
#   theme(legend.position = "none", legend.title = element_blank()) +
#   theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
#   theme(axis.text.y= element_text(colour="Black"), axis.title.y= element_text(colour="Black")) +
#   theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
#   theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) #+
# #theme(plot.title=element_text(size=28, colour="Black"))
# crit1_boxplot

if (write_it==TRUE){
  ggsave(file="crit1_trials_wdots.svg", plot=crit1_bar, width=4, height=6)
  #ggsave(file="crit1_trials_boxplot.svg", plot=crit1_boxplot, width=4, height=6)
}

# ## bar graph: fewer trials are necessary the second time around for both overtrained groups
# crits_bar <- ggplot(filter(crit_table, Group!="CT"), aes(x=Group, y=CritTrials, group=Day_Block)) + #scale_y_continuous(limits = c(-0.2,.7))+
#   geom_bar(aes(x=Group, y=CritTrials, fill=Day_Block, color=Group, group=Day_Block), width=0.4, size = 2, stat="identity", position=position_dodge(.5)) +
#   geom_errorbar(aes(ymin=CritTrials-se, ymax=CritTrials+se), colour="black", width=0, position=position_dodge(.5), size = 1.5) +
#   xlab("Group") + ylab("Criterion training trials needed Day 1")+ 
#   scale_color_manual(values=grps_cols_order2) + scale_shape_manual(values=c(21,22)) + 
#   scale_fill_manual(values=c("white", "white"))+
#   theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
#   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
#   theme(panel.border = element_blank(), axis.line = element_line())+
#   theme(legend.position = "none", legend.title = element_blank()) +
#   theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
#   theme(axis.text.y= element_text(colour="Black"), axis.title.y= element_text(colour="Black")) +
#   theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
#   theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) #+
# #theme(plot.title=element_text(size=28, colour="Black"))
# crits_bar

## ~ analyses: criterion training block length ####
# at 1_1, no differences between the groups (3 groups)
crit1comp <- lm(CritTrials~Group, data=filter(crit_N, Day_Block=="1_1"))
anova(crit1comp)
# using only 2 groups, effect of Day BUT NO INTERACTION between day and group
critfit <- lmer(CritTrials~Group*Day_Block+(1|Subject), data=filter(crit_N, Group!="CT"))
anova(critfit)
testInteractions(critfit, fixed="Group", pairwise="Day_Block")
testInteractions(critfit, pairwise=c("Day_Block", "Group"))

## ~ figure how many ppl still "in running" at different timepoint (for supplement) ####
# Conclusion: remove first 2 trials (warmup) and trials over 55 (> half of subjects out of pool in at least 1)
crit$Trial <- as.factor(crit$Trial)
crit_data_overtime_ <-summarySEwithin(filter(crit, Group!="CT", Missed==0, OptimalChoice==1), measurevar="RT", withinvars=c("Group","Day_Block","Trial"),conf.interval = .95)
grps_cols_order1 = c("#009900", "#FF8000","#0066CC") # CT EOT COT
stable_overtime_deval <- ggplot(crit_data_overtime_, aes(x=as.numeric(Trial), y=N)) + scale_x_continuous(limits = c(1, 120))+
  geom_line(aes(x=as.numeric(Trial), y=N, color = Group, group=Group), size =1)+
  xlab("Trial") + ylab("Number of datapoints")+
  geom_hline(yintercept=22)+
  geom_vline(xintercept=3)+geom_vline(xintercept=55)+
  facet_grid(.~Day_Block)+
  theme_bw(base_size=16)+ scale_color_manual(values=c("#0066CC","#FF8000"))+
  theme(legend.position = "bottom") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
stable_overtime_deval

if (write_it==TRUE){
  ggsave(file="Suppl_n_trials_crit.svg", plot=stable_overtime_deval, width=10, height=5)
}

## ~ Rolling averages of criterion training on days 1 and 4 (EOT and COT groups) ####
# label condition so the initial criterion training for EOT and COT are together
crit$Condition[crit$Day_Block=="1_1"] <- "Crit1"
crit$Condition[crit$Day_Block=="4_1"&crit$Group=="EOT"] <- "OT"
crit$Condition[crit$Day_Block=="4_1"&crit$Group=="COT"] <- "Crit2"
crit$Condition2[crit$Day_Block=="1_1"&crit$Group=="EOT"] <- "Crit1_EOT"
crit$Condition2[crit$Day_Block=="1_1"&crit$Group=="COT"] <- "Crit1_COT"
crit$Condition2[crit$Day_Block=="4_1"&crit$Group=="EOT"] <- "OT"
crit$Condition2[crit$Day_Block=="4_1"&crit$Group=="COT"] <- "Crit2"

crit_rt_overtime_ <-summarySEwithin(filter(crit, Group!="CT", Missed==0, OptimalChoice==1), measurevar="RT", withinvars=c("Condition","Trial"),conf.interval = .95)
crit_rt_overtime <- filter(crit, Group!="CT", Missed==0, OptimalChoice==1) %>% group_by(Condition, Trial) %>% dplyr::summarize(meanRT=mean(RT, na.rm=TRUE), seRT=(sd(RT, na.rm=TRUE)/sqrt(43))) %>% left_join(crit_rt_overtime_)
crit_acc_overtime_ <-summarySEwithin(filter(crit, Missed==0), measurevar="OptimalChoice", withinvars=c("Condition","Trial"),conf.interval = .95)
crit_acc_overtime <- filter(crit, Group!="CT", Missed==0) %>% group_by(Condition, Trial) %>% dplyr::summarize(meanACC=mean(OptimalChoice, na.rm=TRUE), seACC=(sd(OptimalChoice, na.rm=TRUE)/sqrt(43))) %>% left_join(crit_acc_overtime_)

crit_rt_overtime = filter(crit_rt_overtime, as.numeric(Trial)>2, as.numeric(Trial)<56)
crit_acc_overtime = filter(crit_acc_overtime, as.numeric(Trial)>2, as.numeric(Trial)<56)

# create a running average over 8 trials
ntrials <- 53
rollmean_rtcorrect<-data.matrix(3*ntrials,1)
rollmean_sertcorrect<-data.matrix(3*ntrials,1)
rollmean_acc<-data.matrix((3*ntrials),1)
rollmean_seacc<-data.matrix((3*ntrials),1)

for (q in 0:2){
  rollmean_rtcorrect[(q*ntrials+1):(q*ntrials+ntrials)] <- roll_mean((crit_rt_overtime$meanRT[(q*ntrials+1):(q*ntrials+ntrials)]), n=8, fill=c(NA, NA, NA))
  rollmean_sertcorrect[(q*ntrials+1):(q*ntrials+ntrials)] <- roll_mean((crit_rt_overtime$se[(q*ntrials+1):(q*ntrials+ntrials)]), n=8, fill=c(NA, NA, NA))
  rollmean_acc[(q*ntrials+1):(q*ntrials+ntrials)] <- roll_mean((crit_acc_overtime$meanACC[(q*ntrials+1):(q*ntrials+ntrials)]), n=8, fill=c(NA, NA, NA))
  rollmean_seacc[(q*ntrials+1):(q*ntrials+ntrials)] <- roll_mean((crit_acc_overtime$se[(q*ntrials+1):(q*ntrials+ntrials)]), n=8, fill=c(NA, NA, NA))
}
crit_rt_overtime$rollmean_rtcorrect <- rollmean_rtcorrect
crit_rt_overtime$rollmean_sertcorrect <- rollmean_sertcorrect
crit_acc_overtime$rollmean_acc <- rollmean_acc
crit_acc_overtime$rollmean_seacc <- rollmean_seacc

# Figures rolling averages : learning during criterion training blocks
col = c("grey30", "#0066CC", "#FF8000") # CT EOT COT
roll_rt <- ggplot(filter(crit_rt_overtime, Condition!="qOT_crit"), aes(x=as.numeric(Trial), y=rollmean_rtcorrect, group=Condition, colour = Condition)) + scale_y_continuous(limits = c(0.3, 0.6))+
  geom_ribbon(aes(ymin=rollmean_rtcorrect-rollmean_sertcorrect, ymax=rollmean_rtcorrect+rollmean_sertcorrect, fill=Condition, group=Condition), alpha=0.2, colour=NA) +
  geom_line(aes(color=Condition), size = 1.75, stat="identity") +
  xlab("Trial") + ylab("Response Time (s)\n(rolling average over 8 trials)")+
  theme_classic(base_size=16)+ scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14, colour="Black"), axis.title.x = element_text(size=14, colour="Black")) +
  theme(axis.text.y= element_text(size=14, colour="Black"),axis.title.y= element_text(size=14, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=1,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=14, colour="Black"))  
roll_rt

roll_acc <- ggplot(crit_acc_overtime, aes(x=as.numeric(Trial), y=rollmean_acc, colour = Condition)) + #scale_x_continuous(limits = c(0,58))+
  geom_ribbon(aes(ymin=rollmean_acc-rollmean_seacc, ymax=rollmean_acc+rollmean_seacc, fill=Condition, group=Condition), alpha=0.2, colour=NA) +
  geom_line(aes(x=as.numeric(Trial), y=rollmean_acc, color=Condition), size = 1.75, stat="identity") +
  xlab("Trial") + ylab("Accuracy\n(rolling average over 8 trials)")+#ggtitle("Crit Training accuracy")+
  theme_classic(base_size=16)+ scale_color_manual(values=col)+scale_fill_manual(values=col)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14, colour="Black"), axis.title.x = element_text(size=14, colour="Black")) +
  theme(axis.text.y= element_text(size=14, colour="Black"),axis.title.y= element_text(size=14, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=1,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=14, colour="Black"))  
roll_acc

if (write_it==TRUE){
  ggsave(file="crit_rolling_averages8_rt.svg", plot=roll_rt, width=9, height=6)
  ggsave(file="crit_rolling_averages8_acc.svg", plot=roll_acc, width=9, height=6)
}

## ~ Analysis on the RT's - bootstrapped exponential fit ####
# make wide datasets that are 1 row per person, just the RT
crit <- mutate(crit, Cond = paste0(Group, Day_Block))
cons <- unique(crit$Cond)
rt_wide <- list("COT1"=matrix(), "COT4" = matrix(), "EOT1"=matrix(), "EOT4"=matrix(), "CT1"=matrix())

for (c in 1:length(cons)) {
  dat_cond <- filter(crit, Cond==cons[c])
  subj <- unique(dat_cond$Subject)
  rt_wide[[c]] <- matrix(NA, nrow=length(subj), ncol=53)
  for (n in 1:length(subj)) {
    l = filter(dat_cond, (Subject == subj[n]))
    l<-l[3:55,]
    r <- dplyr::select(l, RT)
    ac = dplyr::select(l, OptimalChoice)
    rt_wide[[c]][n,] <- t(r)
  }
  rm(dat_cond, subj)
}

# simple bootstrap per group
boot_rt <- data.frame("Group"=c("COT1", "COT4", "EOT1", "EOT4", "CT"), "c" = rep(NA, 5), "clower"=rep(NA, 5), "cupper"=rep(NA, 5), "lin"=rep(NA, 5), "linlower"=rep(NA, 5), "linupper"=rep(NA, 5))
all_rt_samples <- list("COT1"=matrix(), "COT4" = matrix(), "EOT1"=matrix(), "EOT4"=matrix(), "CT1"=matrix())


for (c in 1:length(cons)) {
  boot <- boot(data = rt_wide[[c]], statistic = exp_fit, R=10000)
  bci_c <- boot.ci(boot, type="bca", index=1)
  bci_lin <- boot.ci(boot, type="bca", index=2)
  boot_rt$c[c] <- boot$t0[1]
  boot_rt$clower[c] <- bci_c$bca[4]
  boot_rt$cupper[c] <- bci_c$bca[5]
  boot_rt$lin[c] <- boot$t0[2]
  boot_rt$linlower[c] <- bci_lin$bca[4]
  boot_rt$linupper[c] <- bci_lin$bca[5]
  all_rt_samples[[c]] <-boot$t
  rm(boot, bci_c, bci_lin)
}

## ~ Analysis on the accuracy: fit logistic (linear) regression ####
cases <- distinct(crit, Subject, Cond)

# Run a logistic regression per "condition" (day/group) to get estimates and standard error, and with all groups, to estimate
# the differences between the groups
crit$Trial<- as.numeric(crit$Trial)
crit$Cond<-relevel(as.factor(crit$Cond), ref="COT4_1")
overal_log <- glm(OptimalChoice~Trial*Cond, data=filter(crit, Trial<56),family = binomial(link = "logit") )
testInteractions(overal_log, pairwise = "Cond")
summary(overal_log)
# log per group (make dataframe to hold results, and run through the 4 analyses)
acc_log_summary <- data.frame("Cond"<-cons)
acc_log_summary$slope <- NA
acc_log_summary$slope_se <-NA
colnames(acc_log_summary) <- c("Cond", "slope", "slope_se")
acc_log_summary$intercept <- NA
acc_log_summary$intercept_se <-NA

for (c in (1:length(cons))) {
  dat <- filter(crit, Cond==cons[c], Trial<56)
  logfit <- glm(OptimalChoice~Trial, data=dat,family = binomial(link = "logit"))
  acc_log_summary$slope[c] <- coef(logfit)[2]
  acc_log_summary$slope_se[c]<-se.coef(logfit)[2]
  acc_log_summary$intercept[c] <- coef(logfit)[1]
  acc_log_summary$intercept_se[c]<-se.coef(logfit)[1]
  rm(logfit, dat)
}

## ~ Figures for fitted learning (accuracy and RT) ####
cols <-c("Grey", "Grey","#0066CC", "#FF8000")
## . . accuracy ####
acc_log_summary$Cond2 <- factor(acc_log_summary$Cond, levels = c("COT1_1","EOT1_1","COT4_1", "EOT4_1", "CT1_1"))
bar_log_acc_slope <- ggplot(filter(acc_log_summary, Cond!="CT1_1"), aes(x=Cond2, y=slope)) + scale_y_continuous(limits = c(-.005,.045))+
  geom_bar(aes(x=Cond2, y=slope, fill=Cond2), width=.5, stat="identity") +
  geom_errorbar(aes(ymin=slope-slope_se, ymax=slope+slope_se), colour="black", width=0, position=position_dodge(width=0), size =1) +
  xlab("Group") + ylab("Logistic regression\nSlope parameter")+
  theme_bw(base_size=16)+ scale_fill_manual(values=cols)+
  scale_x_discrete(labels=c("COT1_1" = "COT\ninitial\ntraining", "EOT1_1" = "EOT\ninitial\ntraining", "COT4_1"= "COT\nnew\nmapping", "EOT4_1" = "EOT\novertrained\nmapping"))+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_log_acc_slope

bar_log_acc_intercept <- ggplot(filter(acc_log_summary, Cond!="CT1_1"), aes(x=Cond2, y=intercept)) + scale_y_continuous(limits = c(0,3))+
  geom_bar(aes(x=Cond2, y=intercept, fill=Cond2), width=.5, stat="identity") +
  geom_errorbar(aes(ymin=intercept-intercept_se, ymax=intercept+intercept_se), colour="black", width=0, position=position_dodge(width=0), size =1) +
  xlab("Group") + ylab("Logistic regression\nIntercept")+#ggtitle("logistic regression on accuracy")+
  theme_bw(base_size=16)+ scale_fill_manual(values=cols)+
  scale_x_discrete(labels=c("COT1_1" = "COT\ninitial\ntraining", "EOT1_1" = "EOT\ninitial\ntraining", "COT4_1"= "COT\nnew\nmapping", "EOT4_1" = "EOT\novertrained\nmapping"))+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_log_acc_intercept

## . . RT ####
boot_rt$Group2 <- factor(boot_rt$Group, levels = c("COT1","EOT1","COT4", "EOT4", "CT1"))
bar_boot_rt <- ggplot(filter(boot_rt, Group!="CT"), aes(x=Group2, y=c)) + scale_y_continuous(limits = c(-0.05,.35))+
  geom_bar(aes(x=Group2, y=c, fill=Group2), width=.5, stat="identity") +
  geom_errorbar(aes(ymin=clower, ymax=cupper), colour="black", width=0, position=position_dodge(width=0), size =1) +
  xlab("Group") + ylab("Exponential fit\nExponent (speed of decrease)")+#ggtitle("RT Exponential decline")+
  scale_color_manual(values=cols) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=cols)+ 
  scale_x_discrete(labels=c("COT1" = "COT\ninitial\ntraining", "EOT1" = "EOT\ninitial\ntraining", "COT4"= "COT\nnew\nmapping", "EOT4" = "EOT\novertrained\nmapping"))+
  theme_bw(base_size=16)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_boot_rt

# the linear component 
bar_boot_rtlin <- ggplot(filter(boot_rt, Group!="CT"), aes(x=Group2, y=lin)) + scale_y_continuous(limits = c(-.2,1.7))+
  geom_bar(aes(x=Group2, y=lin, fill=Group2), width=.5, stat="identity") +
  geom_errorbar(aes(ymin=linlower, ymax=linupper), colour="black", width=0, position=position_dodge(width=0), size =1) +
  xlab("Group") + ylab("Exponential fit\nLinear component")+#ggtitle("linear component")+
  scale_color_manual(values=cols) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=cols)+ 
  scale_x_discrete(labels=c("COT1" = "COT\ninitial\ntraining", "EOT1" = "EOT\ninitial\ntraining", "COT4"= "COT\nnew\nmapping", "EOT4" = "EOT\novertrained\nmapping"))+
  theme_bw(base_size=16)+ #scale_fill_manual(values=barcols)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_boot_rtlin

if (write_it==TRUE){
  ggsave(file="crit_log_acuracy_slope.svg", plot=bar_log_acc_slope, width=4.5, height=6)
  ggsave(file="crit_log_accuracy_intercept.svg", plot=bar_log_acc_intercept, width=4.5, height=6)
  ggsave(file="crit_boot_rt_slope.svg", plot=bar_boot_rt, width=4.5, height=6)
  ggsave(file="crit_boot_rt_intercept.svg", plot=bar_boot_rtlin, width=4.5, height=6)
}

## . . RT: figures to show significance based on bootstrapping ####

grps_cols_order1 = c("#009900", "#FF8000","#0066CC") # CT EOT COT
samples_cot1 <- data.frame("c" = all_rt_samples[[1]][,1], "lin" = all_rt_samples[[1]][,2])
samples_cot4 <- data.frame("c" = all_rt_samples[[2]][,1], "lin" = all_rt_samples[[2]][,2])
samples_eot1 <- data.frame("c" = all_rt_samples[[3]][,1], "lin" = all_rt_samples[[3]][,2])
samples_eot4 <- data.frame("c" = all_rt_samples[[4]][,1], "lin" = all_rt_samples[[4]][,2])
samples_ct <- data.frame("c" = all_rt_samples[[5]][,1], "lin" = all_rt_samples[[5]][,2])
# histogram to derive "p-value"
hist_cot1 <- ggplot(samples_cot1) +scale_x_continuous(limits = c(-.15,.3))+scale_y_continuous(limits = c(0,700))+
  geom_histogram(aes(x=c), color="grey", fill="gray", bins=100) + 
  geom_vline(aes(xintercept=0.2046), size=1.5, color="#0066CC") +
  xlab("#samples") + xlab("exponent estimate")+
  theme_bw(base_size=20)+ 
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
hist_cot1

hist_eot1 <- ggplot(samples_eot1) +scale_x_continuous(limits = c(-.15,.3))+scale_y_continuous(limits = c(0,700))+
  geom_histogram(aes(x=c), color="grey", fill="gray", bins=100) + 
  geom_vline(aes(xintercept=0.2046), size=1.5, color="#0066CC") +
  xlab("#samples") + xlab("slope estimate")+
  theme_bw(base_size=20)+ 
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
hist_eot1

hist_eot4 <- ggplot(samples_eot4) +scale_x_continuous(limits = c(-.15,.3))+scale_y_continuous(limits = c(0,1000))+
  geom_histogram(aes(x=c), color="#FF8000", fill="#FF8000", bins=100) + 
  geom_vline(aes(xintercept=0.2046), size=1.5, color="#0066CC") +
  xlab("#samples") + xlab("slope estimate")+#ggtitle("Exponential increase")+
  theme_bw(base_size=20)+ 
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
hist_eot4

if (write_it==TRUE){
  lay <- rbind(c(1,2,3))
  svg("bootstrapped_significance.svg", width = 24, height = 6)
  pl <- grid.arrange(hist_cot1, hist_eot1, hist_eot4,layout_matrix = lay) 
  dev.off()
}
