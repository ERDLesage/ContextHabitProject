#########################################################################################~##
## Context habit project
## Script 1 of 3: From raw data to processed data, ready for analysis and figures ####
## 2021 Elise Lesage (elise.lesage@ugent.be) 
## --------------------------------------------------------------------------------------~--
## Note: 
## - After setting a local directory and toggling whether you would like to (over)write, running this script 
##   should write you an R workspace file with the dataframes needed for the analysis and figures script
## - Excluse the mess, this is a (only slightly cleaned) working document :)
##

rm(list = ls()) # clear wm 
##########################################################################################~##
## SET YOUR PATHS AND SET WHETHER YOU WANT OUTPUT
## ---------------------------------------------------------------------------------------~--
# Set your paths 
datadir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/scripts_for_sharing/data"
outputdir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/scripts_for_sharing/analysis_and_output"

# Switch writing on (TRUE) or off (FALSE) 
write_it = TRUE

## FROM HERE ON OUT, JUST RUNNING THE SCRIPT SHOULD WORK
##########################################################################################~##
install.packages("ggplot2", "doBy", "gtable", "grid", "gridExtra", "afex", "phia", "lattice", "latticeExtra") 

library(doBy)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(afex)
library(phia)
library(lattice)
library(latticeExtra)
library(Rmisc)
library(lmerTest)
library(RcppRoll)
library(scales)
library(dplyr)
library(tidyr)

# read in data
setwd(datadir)
RL <- read.table("Context_Habit_data_anonymised.txt", header = TRUE, sep = ",")
modelfits <- read.table("Context_Habit_MBMF_modelfits_anonymised.txt", header=TRUE, sep=',')
wmdata <- read.table("Context_Habit_wmloads_anonymised.txt", header = TRUE, sep = ",")
twostep_vals = read.table("MBMF_task_sample_outcomes.txt", header = TRUE, sep = "\t")

## some basic edits (RL) ####
RL$Congruent = factor(RL$Congruent)
RL$Congruency[RL$Congruent==-1] = "Devalued"
RL$Congruency[RL$Congruent==1] = "Valued"
RL$Num_Days = factor(RL$Num_Days)
RL$Group[RL$Num_Days == 1] <- "COT"
RL$Group[RL$Num_Days == 4] <- "EOT"
RL$Group[RL$Num_Days == 0] <- "CT"
# ensure that accuracy, and (non-WM)missed are coded NA when response is missed on WM trial
RL$OptimalChoice[(RL$WM_Trial==1&RL$WM_correct==666)]<-NA
RL$Missed[(RL$WM_Trial==1&RL$WM_correct==666)]<-NA

## Automaticity ####
simplified_rl <- filter(RL, WM_correct != 0, Missed==0, Group!= "CT", Day_Block!="4_1", Day_Block!="4_3", Day_Block != "4_4")
rl_simplified_summary <-summarySE(simplified_rl, measurevar="OptimalChoice", groupvars = c("Group", "Subject", "Block", "Congruency"), conf.interval = .95)
rl_simplified_summary_rt <-summarySE(simplified_rl, measurevar="RT", groupvars = c("Group", "Subject", "Block", "Congruency"), conf.interval = .95)
firstandlast<- filter(simplified_rl, (Day_Block=="1_3" | Day_Block=="3_6"))

summary_1andlast_acc <-summarySE(firstandlast, measurevar="OptimalChoice", groupvars = c("Subject", "Day_Block"), conf.interval = .95)
summary_1andlast_acc <- dplyr::select(summary_1andlast_acc, -N, -se, -sd, -ci)
summary_1andlast_rt <-summarySE(firstandlast, measurevar="RT", groupvars = c("Subject", "Day_Block"), conf.interval = .95)
summary_1andlast_rt <- dplyr::select(summary_1andlast_rt, -N, -se, -sd, -ci)
summary_1andlast_sdrt <-summarySE(firstandlast, measurevar="RT", groupvars = c("Subject", "Day_Block"), conf.interval = .95)
summary_1andlast_sdrt <- dplyr::select(summary_1andlast_sdrt, -RT, -N, -se, -ci)

# now make it a wide format (for nonparametric test)
ot_summary_wide_acc <- spread(summary_1andlast_acc, Day_Block, OptimalChoice, sep="_")
ot_summary_wide_rt <- spread(summary_1andlast_rt, Day_Block, RT , sep="_")
ot_summary_wide_sdrt <- spread(summary_1andlast_sdrt, Day_Block, sd, sep="_")

## Automaticity for figures 
rl_accuracyg <-summarySEwithin(filter(RL, Group != "CT", WM_correct >0, Missed==0), measurevar="OptimalChoice", idvar = "Subject", withinvars=c("Day_Block", "Congruency"),conf.interval = .95)
rl_accg <- filter(RL, Group != "CT", WM_correct >0, Missed==0) %>% group_by(Day_Block,Congruency) %>% dplyr::summarize(ACC=mean(OptimalChoice, na.rm=TRUE), SE = (sd(OptimalChoice, na.rm=TRUE)/sqrt(85))) %>% left_join(rl_accuracyg)
rm(rl_accuracyg)
rl_reactiontimeg <-summarySEwithin(filter(RL, Group != "CT", OptimalChoice ==1,WM_correct >0, Missed==0), measurevar="RT", idvar = "Subject", withinvars=c("Day_Block", "Congruency"),conf.interval = .95)
rl_RTg <- filter(RL, OptimalChoice ==1, Group != "CT",WM_correct >0, Missed==0) %>% group_by(Day_Block,Congruency) %>% dplyr::summarize(ReactionTime=mean(RT, na.rm=TRUE), SERT = (sd(RT, na.rm=TRUE)/sqrt(85))) %>% left_join(rl_reactiontimeg)
rm(rl_reactiontimeg)
rl_sdrt <- filter(RL, WM_correct != 0, Group != "CT", Missed==0) %>% group_by(Subject,Day_Block,Congruency) %>% dplyr::summarize(SD=sd(RT, na.rm=TRUE)) 
rl_RTSDg <- rl_sdrt %>% group_by(Day_Block,Congruency) %>% dplyr::summarize(ReactionTimeSD=mean(SD, na.rm=TRUE), RTSD_se=(sd(SD, na.rm=TRUE)/sqrt(85))) 
rm(rl_sdrt)


## Devaluation: line-per-trial dataframes for analyses ####
RL_alldeval <- filter(RL, Deval_Trial == 1 )
RL_alldeval$WM_Load <- ifelse(RL_alldeval$WM_Trial==1, "WM Load", "No Load")
RL_deval <- filter(RL, (Deval_Trial == 1 & WM_Trial==0))
RL_deval_wm <- filter(RL, (Deval_Trial == 1 & WM_Trial==1))

## Summaries - line per group/block/condition
rl_accuracy_betweenerrorbars <-summarySE(filter(RL, WM_correct > 0, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "Deval_Trial", "WM_Trial" , "Day_Block", "Congruency"), conf.interval = .95)
rl_accuracy_withinerrorbars <-summarySEwithin(filter(RL, WM_correct >0, Missed==0), measurevar="OptimalChoice", betweenvars = c("Group", "Subject"), withinvars=c("Deval_Trial", "WM_Trial","Day_Block", "Congruency"),conf.interval = .95)
#rl_acc <- filter(RL, Missed==0) %>% group_by(Group, Subject, factor(Deval_Trial), factor(WM_Trial),Day_Block,Congruency) %>% dplyr::summarize(ACC=mean(OptimalChoice, na.rm=TRUE)) %>% left_join(rl_accuracy_withinerrorbars)

rl_reactiontime_betweenerrorbars <-summarySE(filter(RL, OptimalChoice ==1, WM_correct != 0, Missed==0), measurevar="RT", groupvars = c("Group", "Subject", "Deval_Trial", "WM_Trial", "Day_Block", "Congruency"), conf.interval = .95)
rl_reactiontime_withinerrorbars <-summarySEwithin(filter(RL, OptimalChoice ==1,WM_correct >0, Missed==0), measurevar="RT", betweenvars = c("Group", "Subject"), withinvars=c("Deval_Trial", "WM_Trial", "Day_Block", "Congruency"),conf.interval = .95)
rl_RT <- filter(RL, OptimalChoice ==1,WM_correct >0, Missed==0) %>% group_by(Group, Subject, factor(Deval_Trial), factor(WM_Trial), Day_Block,Congruency) %>% dplyr::summarize(meanRT=mean(RT, na.rm=TRUE)) %>% left_join(rl_reactiontime_withinerrorbars, by=c("Group", "Subject","Day_Block", "Congruency"))

rldeval_accuracy_betweenerrorbars <-summarySE(filter(RL_alldeval, WM_correct > 0, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
rldeval_accuracy_withinerrorbars <-summarySEwithin(filter(RL_alldeval, WM_correct >0, Missed==0), measurevar="OptimalChoice", betweenvars = c("Group", "Subject"), withinvars=c("WM_Load", "Congruency"),conf.interval = .95)
rldeval_acc <- filter(RL_alldeval, WM_correct >0, Missed==0) %>% group_by(Group, Subject,WM_Load,Congruency) %>% dplyr::summarize(ACC=mean(OptimalChoice, na.rm=TRUE)) %>% left_join(rldeval_accuracy_withinerrorbars)

rldeval_reactiontime_betweenerrorbars <-summarySE(filter(RL_alldeval, OptimalChoice ==1, WM_correct != 0, Missed==0), measurevar="RT", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
rldeval_reactiontime_withinerrorbars <-summarySEwithin(filter(RL_alldeval, OptimalChoice ==1,WM_correct >0, Missed==0), measurevar="RT", betweenvars = c("Group", "Subject"), withinvars=c("WM_Load", "Congruency"),conf.interval = .95)
rldeval_RT <- filter(RL_alldeval, OptimalChoice ==1,WM_correct >0, Missed==0) %>% group_by(Group, Subject, WM_Load,Congruency) %>% dplyr::summarize(meanRT=mean(RT, na.rm=TRUE), medianRT = median(RT, na.rm=TRUE)) %>% left_join(rldeval_reactiontime_withinerrorbars)

# Criterion training trials - both day 1 and day 4 ####
crit <- filter(RL, Day_Block =="1_1" | Day_Block == "4_1")

## Make line-per-subject (LPS) dataframe ####
rl_acc_deval <- filter(rl_accuracy_betweenerrorbars, (Deval_Trial == 1 & WM_Trial ==0))
rl_acc_deval_wm <- filter(rl_accuracy_betweenerrorbars, (Deval_Trial == 1 & WM_Trial ==1))
deval_wide_acc <- reshape(rl_acc_deval, timevar = "Congruency", drop = c("sd", "ci","se", "N"), v.names = c("OptimalChoice"), idvar = c("Subject"), direction = "wide")
deval_wide_acc$deval_diff_acc <-deval_wide_acc$OptimalChoice.Valued-deval_wide_acc$OptimalChoice.Devalued
deval_wide_acc$deval_mean_acc <- (deval_wide_acc$OptimalChoice.Valued+deval_wide_acc$OptimalChoice.Devalued)/2

rl_rt_deval <- filter(rl_reactiontime_betweenerrorbars, (Deval_Trial == 1 & WM_Trial ==0))
rl_rt_deval_wm <- filter(rl_reactiontime_betweenerrorbars, (Deval_Trial == 1 & WM_Trial ==1))
deval_wide_rt <- reshape(rl_rt_deval, timevar = "Congruency", drop = c("sd", "ci","se", "N"), v.names = c("RT"), idvar = c("Subject"), direction = "wide")
deval_wide_rt$deval_diff_rt <-deval_wide_rt$RT.Valued-deval_wide_rt$RT.Devalued

devalwm_wide_acc <- reshape(rl_acc_deval_wm, timevar = "Congruency", drop = c("RT","sd", "ci","se", "N"), v.names = c("OptimalChoice"), idvar = c("Subject"), direction = "wide")
devalwm_wide_acc$devalwm_diff_acc <-devalwm_wide_acc$OptimalChoice.Valued-devalwm_wide_acc$OptimalChoice.Devalued
devalwm_wide_acc$devalwm_mean_acc <-(devalwm_wide_acc$OptimalChoice.Valued-devalwm_wide_acc$OptimalChoice.Devalued)/2
devalwm_wide_rt <- reshape(rl_rt_deval_wm, timevar = "Congruency", drop = c("sd", "ci","se", "N"), v.names = c("RT"), idvar = c("Subject"), direction = "wide")
devalwm_wide_rt$devalwm_diff_rt <-devalwm_wide_rt$RT.Valued-devalwm_wide_rt$RT.Devalued

deval_effects <- dplyr::select(deval_wide_acc, Group, Subject, deval_diff_acc)
deval_effects$deval_diff_rt <- deval_wide_rt$deval_diff_rt
deval_effects$devalwm_diff_acc <- devalwm_wide_acc$devalwm_diff_acc
deval_effects$devalwm_diff_rt <- devalwm_wide_rt$devalwm_diff_rt
deval_effects$devalwm_mean_acc <-devalwm_wide_acc$devalwm_mean_acc
deval_effects$deval_mean_acc <- deval_wide_acc$deval_mean_acc
LPS__ <- left_join(deval_effects, modelfits, by="Subject")
LPS <- left_join(LPS__, wmdata, by="Subject")
LPS$Group[LPS$Num_Days == 1] <- "COT"
LPS$Group[LPS$Num_Days == 4] <- "EOT"
LPS$Group[LPS$Num_Days == 0] <- "CT"

##################################################################################################~##
## SAVE THE NEEDED DATAFRAMES IN A WORKSPACE ####
setwd(outputdir)

if(write_it==TRUE) {
  save(ot_summary_wide_acc, ot_summary_wide_rt, ot_summary_wide_sdrt, rl_accg, rl_RTg, rl_RTSDg, RL_deval,
     RL_alldeval, RL_deval_wm, RL, LPS, crit, rldeval_acc, twostep_vals, file = "Context_habit_processed_dataframes.RData") 
  }
