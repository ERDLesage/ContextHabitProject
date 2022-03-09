#########################################################################################~##
## Context habit project
## Script 2 of 3: Loads preprocessed dataframes and does analysis and figures ####
## 2021 Elise Lesage (elise.lesage@ugent.be) 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Note: 
## - First set a local directory for data/output and toggle whether you would like to (over)write figures. 
## - Needs the R workspace file created by the preprocessing script. 
## - The script perfoms  analyses and writes figures. (The output from the analyses is not written.)
## - generally all the data needed is either in the preprocessed data (R workspace created by first script)
##   or created locally to the analysis/figure. 
## - The analyses and figures for the criterion training blocks (both initial, and the first block on Day4) 
##   are handled by a third script.

rm(list = ls()) # clear wm 
##~########################################################################################~##
## SET YOUR PATHS, SET WHETHER YOU WANT OUTPUT ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set your paths 
datadir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/scripts_for_sharing/analysis_and_output"
outputdir = "C:/Users/elise.000/OneDrive/Documents/Habit_Project/scripts_for_sharing/analysis_and_output"

# Switch writing (figures) on (TRUE) or off (FALSE) 
write_it = TRUE

## FROM HERE ON OUT, NO EDITS SHOULD BE NEEDED ####
#########################################################################################~###
library(doBy)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(afex)
library(phia)
library(lattice)
library(latticeExtra)
library(dplyr)
library(tidyr)
library(Rmisc)
library(lmerTest)
library(RcppRoll)
library(scales)
library(minpack.lm)
library(svglite)
library(effectsize)

# load workspace
setwd(datadir)
load("Context_habit_processed_dataframes.RData")

# function to calculate pseudo-effect sizes (http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#how-do-i-compute-a-coefficient-of-determination-r2-or-an-analogue-for-glmms)
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

##~###########################################################################################~##
## ANALYSES ####
##~###########################################################################################~##

## ~ overtraining: change in accuracy, rt, sdrt ####
wilcox.test(ot_summary_wide_acc$Day_Block_1_3, ot_summary_wide_acc$Day_Block_3_6, paired = TRUE, alternative = "two.sided")
wilcox.test(ot_summary_wide_rt$Day_Block_1_3, ot_summary_wide_rt$Day_Block_3_6, paired = TRUE, alternative = "two.sided")
wilcox.test(ot_summary_wide_sdrt$Day_Block_1_3, ot_summary_wide_sdrt$Day_Block_3_6, paired = TRUE, alternative = "two.sided")


## ~ devaluation no load ####
# relevel congruency and group
RL_deval$Congruency <- relevel(as.factor(RL_deval$Congruency), ref="Valued")
RL_deval_wm$Congruency <- relevel(as.factor(RL_deval_wm$Congruency), ref="Valued")
RL_alldeval$Congruency <- relevel(as.factor(RL_alldeval$Congruency), ref="Valued")
RL_deval$Group <- relevel(as.factor(RL_deval$Group), ref="CT")
RL_deval_wm$Group <- relevel(as.factor(RL_deval_wm$Group), ref="CT")
RL_alldeval$Group <- relevel(as.factor(RL_alldeval$Group), ref="CT")

# accuracy
Fit5_noWML = glmer(OptimalChoice ~ (1|Subject)+Congruency*Group, data=RL_deval, family=binomial) #
summary(Fit5_noWML)
Anova(Fit5_noWML)
testInteractions(Fit5_noWML, pairwise=c("Congruency","Group"))

# RT
FitRT_noWML = lmer(RT ~ (1|Subject)+Congruency*Group, data=filter(RL_deval, OptimalChoice==1)) #
summary(FitRT_noWML)
Anova(FitRT_noWML)
testInteractions(FitRT_noWML, pairwise=c("Congruency","Group"))
testInteractions(FitRT_noWML, fixed="Group", pairwise=c("Congruency"))
testInteractions(FitRT_noWML, pairwise="Group")
testInteractions(FitRT_noWML, fixed="Group", across=c("Congruency"))

## ~ devaluation effect of load ####
# load does not differ according to group
wm_grp_nodiff <-lm(WM_used~Group, data=LPS)
anova(wm_grp_nodiff)
wm_load_sizes <- summarySE(LPS, measurevar="WM_used", groupvars = c("Group"), na.rm="TRUE")

# performance on the WM task itself
wmperf_grp_nodiff <- lm(WM_response~Group, data = filter(RL_deval_wm, WM_response>0, WM_correct!=666))
anova(wmperf_grp_nodiff)

# fyi, just to see if the number of missed responses differs between grps (it doesn't, and is ~2.3%)
RL_deval_wm$WM_Missed <-0
RL_deval_wm$WM_Missed[RL_deval_wm$WM_correct==666]<-1
WM_missed <- summarySE(RL_deval_wm, measurevar = "WM_Missed", groupvars = c("Group"), na.rm="TRUE")
wmmiss_grp_nodiff <- lm(WM_Missed~Group, data = RL_deval_wm)
anova(wmmiss_grp_nodiff)

# accuracy
Fit5 = glmer(OptimalChoice ~ (1|Subject)+Congruency*Group*WM_Load, data=filter(RL_alldeval, WM_correct!=0), family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) # 
Anova(Fit5)
summary(Fit5)
# same group/congruency differences under load
testInteractions(Fit5, fixed=c("WM_Load"), pairwise=c("Group","Congruency"))
# no effect of load on congruency effect
testInteractions(Fit5, fixed="Group", pairwise=c("WM_Load", "Congruency")) 

# RT
FitRT = lmer(RT ~ (1|Subject)+Congruency*Group*WM_Load, data=filter(RL_alldeval, OptimalChoice==1, WM_correct!=0)) # 
Anova(FitRT)
summary(FitRT)
testInteractions(FitRT,pairwise=c("Group", "Congruency", "WM_Load"))
testInteractions(FitRT,fixed="Congruency", pairwise=c("Group",  "WM_Load"))
testInteractions(FitRT,fixed=c("Group"), pairwise=c( "Congruency","WM_Load"))
testInteractions(FitRT,fixed=c("Group", "Congruency"), pairwise=c("WM_Load"))

testInteractions(FitRT,fixed=c("Group"), pairwise=c("WM_Load", "Congruency"))

testInteractions(FitRT,fixed=c("WM_Load"), pairwise=c("Group", "Congruency"))

## ~ devaluation with load ####
# accuracy
Fit5_WML = glmer(OptimalChoice ~ (1|Subject)+Congruency*Group, data=filter(RL_deval_wm, WM_correct !=0), family=binomial) #
summary(Fit5_WML)
Anova(Fit5_WML)
testInteractions(Fit5_WML, pairwise=c("Congruency","Group"))

# RT
FitRT_WML = lmer(RT ~ (1|Subject)+Congruency*Group, data=filter(RL_deval_wm, OptimalChoice==1, WM_correct!=0)) #
summary(FitRT_WML)
Anova(FitRT_WML)
testInteractions(FitRT_WML, pairwise=c("Congruency","Group"))
testInteractions(FitRT_WML, fixed="Group", pairwise=c("Congruency"))
testInteractions(FitRT_WML, pairwise="Group")

## ~ Prediction w~outome-insensitivity ####
# w does not differ between groups
LPS$Group <- as.factor(LPS$Group)
grp_w <- kruskal.test(w ~ Group, data = LPS)
grp_w

# actual regression
LPS$Group <- relevel(as.factor(LPS$Group), ref="CT")
test <- lm(deval_diff_acc~w*Group, data=LPS)
anova(test)
summary(test)
testInteractions(test, pairwise = "Group")

## non-parametric correlations per group (non-parametric because w is very skewed)
# note that an exact p-value is not possible because of many datapoints w same values
LPS_EOT <- filter(LPS, Group== "EOT")
LPS_COT <- filter(LPS, Group== "COT")
LPS_CT <- filter(LPS, Group== "CT")
cor.test(LPS_EOT$w, LPS_EOT$deval_diff_acc, method = "spearman")
cor.test(LPS_COT$w, LPS_COT$deval_diff_acc, method = "spearman")
cor.test(LPS_CT$w, LPS_CT$deval_diff_acc, method = "spearman")

##~###########################################################################################~##
## FIGURES ####
##~###########################################################################################~##
## (Set outputdir) ####
setwd(outputdir)
## (Set colors) ####
grps_cols_order1 = c("#009900", "#FF8000","#0066CC") # CT EOT COT
grps_cols_order2 = c("#FF8000","#0066CC","#009900") # EOT COT CT

## ~ number of training trials ####
# in context, and on specific stimuli that will be devalued
RLtrain <- filter(RL, Deval_Trial==0)
RLtrain$train_devalstim <- 1
RLtrain$train_devalstim[RLtrain$Group=="COT"&RLtrain$Day!=4] <- 0

train_context_ <- summarySE(filter(RLtrain, Missed==0), measurevar = "OptimalChoice", groupvars = c("Group", "Subject"), na.rm=TRUE)
train_devalstim_ <-summarySE(filter(RLtrain, RLtrain$train_devalstim==1), measurevar = "OptimalChoice", groupvars = c("Group", "Subject"), na.rm=TRUE)
train_context_$Amount <- train_context_$N
train_devalstim_$Amount <- train_devalstim_$N
train_context <- summarySE(train_context_, measurevar = "Amount", groupvars = c("Group"), na.rm=TRUE)
train_devalstim <- summarySE(train_devalstim_, measurevar = "Amount", groupvars = c("Group"), na.rm=TRUE)

train_context$Group2 <- factor(train_context$Group, levels = c("CT","EOT","COT"))
train_devalstim$Group2 <- factor(train_devalstim$Group, levels = c("CT","EOT","COT"))


bar_stimtrain <- ggplot(train_devalstim, aes(x=Group2, y=Amount)) + scale_y_continuous(limits = c(0,3200))+
  geom_bar(aes(x=Group2, y=Amount, colour=Group2, fill=Group2), width=.4, size=2, stat="identity") +
  geom_errorbar(aes(ymin=Amount-se, ymax=Amount+se), colour="black", width=0, position=position_dodge(width=0), size =2) +
  ylab("# trials on S-R-O mapping\nbefore devaluation")+
  theme_bw(base_size=20)+ scale_fill_manual(values=grps_cols_order1)+scale_colour_manual(values=grps_cols_order1)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 18), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_stimtrain

bar_contexttrain <- ggplot(train_context, aes(x=Group2, y=Amount)) + scale_y_continuous(limits = c(0,3200))+
  geom_bar(aes(x=Group2, y=Amount, colour=Group2, fill=Group2), width=.4, size=2, stat="identity") +
  geom_errorbar(aes(ymin=Amount-se, ymax=Amount+se), colour="black", width=0, position=position_dodge(width=0), size =2) +
  ylab("# trials in task context\nbefore devaluation")+#ggtitle("logistic regression on accuracy")+
  theme_bw(base_size=20)+ scale_colour_manual(values=grps_cols_order1)+ scale_fill_manual(values=grps_cols_order1)+
  #scale_x_discrete(labels=c("COT1_1" = "COT\ninitial\ntraining", "EOT1_1" = "EOT\ninitial\ntraining", "COT4_1"= "COT\nnew\nmapping", "EOT4_1" = "EOT\novertrained\nmapping"))+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 18), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_contexttrain

if (write_it==TRUE){
  ggsave(file="n_train_stim.svg", plot=bar_stimtrain, width=4.25, height=5)
  ggsave(file="n_train_context.svg", plot=bar_contexttrain, width=4.25, height=5)
}


## ~ Overtraining: accuracy, RT, sdRT ####
rl_ACCg_line <- ggplot(filter(rl_accg, (Day_Block != "4_1" & Day_Block != "4_3" & Day_Block != "4_4")), aes(x=Day_Block, y=ACC)) + scale_y_continuous(limits = c(0.5, 1))+
  geom_line(aes(x=Day_Block, y=ACC, group=Congruency), size = 1, position=position_dodge(.5), stat="identity") +
  geom_errorbar(aes(ymin=ACC-SE, ymax=ACC+SE, group=Congruency), width=0, size =1, position=position_dodge(.5)) +
  geom_point(aes(shape = Congruency, fill=Congruency, group=Congruency), stroke=1, position=position_dodge(.5), size = 4) +
  xlab("Time") + ylab("accuracy")+ggtitle("Accuracy")+
  scale_color_manual(values=c("black", "white") ) + scale_shape_manual(values=c(21, 21)) + scale_fill_manual(values=c("black", "white"))+
  theme_bw(base_size=20)+ 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14, colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(size=14, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=20, colour="Black"))  
rl_ACCg_line

rl_RTg_line <- ggplot(filter(rl_RTg,  (Day_Block != "4_1" & Day_Block != "4_3" & Day_Block != "4_4")), aes(x=Day_Block, y=ReactionTime)) + scale_y_continuous(limits = c(0.35, 0.55))+
  geom_line(aes(x=Day_Block, y=ReactionTime, group=Congruency), size = 1, position=position_dodge(.5), stat="identity") +
  geom_errorbar(aes(ymin=ReactionTime-SERT, ymax=ReactionTime+SERT, group=Congruency), width=0, size =1, position=position_dodge(.5)) +
  geom_point(aes(shape = Congruency, fill=Congruency, group=Congruency, color=Congruency), stroke=2, position=position_dodge(.5), color = "Black", size = 4) +
  xlab("Overtraining") + ylab("response time (s.)")+ggtitle("Response Time")+
  scale_color_manual(values=c("black", "white") ) + scale_shape_manual(values=c(21, 21)) + scale_fill_manual(values=c("black", "white"))+
  theme_bw(base_size=20)+ 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14, colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(size=14, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=20, colour="Black"))  
rl_RTg_line

rl_RTSDg_line <- ggplot(filter(rl_RTSDg, (Day_Block != "4_1" & Day_Block != "4_3" & Day_Block != "4_4")), aes(x=Day_Block, y=ReactionTimeSD)) + scale_y_continuous(limits = c(0.05, 0.2))+
  geom_line(aes(x=Day_Block, y=ReactionTimeSD, group=Congruency), size = 1, position=position_dodge(.5), stat="identity") +
  geom_errorbar(aes(ymin=ReactionTimeSD-RTSD_se, ymax=ReactionTimeSD+RTSD_se, group=Congruency), width=0, size =1, position=position_dodge(.5)) +
  geom_point(aes(shape = Congruency, fill=Congruency, group=Congruency), stroke=1, position=position_dodge(.5), color = "Black", size = 4) +
  xlab("Overtraining") + ylab("sd(RT; s.)")+ggtitle("Response Time variability")+
  scale_color_manual(values=c("black", "white") ) + scale_shape_manual(values=c(21,21)) + scale_fill_manual(values=c("black", "white"))+
  theme_bw(base_size=20)+ 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=14, colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(size=14, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=20, colour="Black"))  
rl_RTSDg_line

## putting three together for export
if (write_it==TRUE){
  lay <- rbind(c(1,2,3))
  svg("automatizationwide.svg", width = 24, height = 6)
  pl <- grid.arrange(rl_ACCg_line, rl_RTg_line, rl_RTSDg_line,layout_matrix = lay) 
  dev.off()
}

## ~ Devaluation insensitivity - dot ####
rldeval_accuracy_g <-summarySEwithin(filter(RL_alldeval, Missed==0), measurevar="OptimalChoice", withinvars=c("Group","WM_Load", "Congruency"),conf.interval = .95)
rldeval_accg <- filter(RL_alldeval, Missed==0) %>% group_by(Group, WM_Load,Congruency) %>% dplyr::summarize(accuracy=mean(OptimalChoice, na.rm=TRUE), SE=(sd(OptimalChoice, na.rm=TRUE))/sqrt(42)) %>% left_join(rldeval_accuracy_g)
rm(rldeval_accuracy_g)
rldeval_accg$Group2 <- factor(rldeval_accg$Group, levels = c("CT","EOT","COT"))
rldeval_accg$Congruency2 <- factor(rldeval_accg$Congruency, levels = c("Valued","Devalued"))

rl_acc_devalpoint <- ggplot(filter(rldeval_accg, WM_Load == "No Load"), aes(x=Group2, y=accuracy, group=Congruency)) + scale_y_continuous(limits = c(0.4, 1.1))+
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
  geom_point(aes(fill=Group2, shape=Congruency, size=Congruency),color="Black", stroke =1, position=position_dodge(.4), stat="identity") +
  xlab("Group") + ylab("Proportion correct")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_order1)+
  theme_bw(base_size=20)+ scale_size_manual(values=c(6,8))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"), axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) #+
#theme(plot.title=element_text(size=28, colour="Black"))
rl_acc_devalpoint

## ~ Outcome insensitivity - bar plus jitter ####
rldeval_accuracy_betweenerrorbars <-summarySE(filter(RL_alldeval, WM_correct > 0, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
dd <- filter(rldeval_accuracy_betweenerrorbars, WM_Load == "No Load")
d_congruent <- filter(dd, Congruency == "Valued")
d_congruent$error_con <- 1-d_congruent$OptimalChoice
d_congruent <- dplyr::select(d_congruent, Group, Subject, error_con, sd)
d_incongruent <- filter(dd, Congruency == "Devalued")
d_incongruent$error_incon <- 1-d_incongruent$OptimalChoice
d_incongruent <- dplyr::select(d_incongruent, Group, Subject, error_incon, sd)
d <- inner_join(d_congruent, d_incongruent, by=c("Subject", "Group"))
d$habit_errorrate <- d$error_incon - d$error_con
dg <- group_by(d, Group) %>% dplyr::summarize(hab_err=mean(habit_errorrate, na.rm=TRUE))

dg$Group2 <- factor(dg$Group, levels = c("CT","EOT","COT"))
d$Group2 <- factor(d$Group, levels = c("CT","EOT","COT"))

dg_point <- ggplot() + scale_y_continuous(limits = c(-0.4,1))+
  geom_bar(data=dg, aes(x=Group2, y=hab_err, group=Group2, fill=Group2), width=.4, stat="identity") +
  geom_point(data = d, aes(x=Group2, y= habit_errorrate, color=Group2),stroke =1, color="Grey30", position=position_jitter(.25), size = 1, stat="identity") +
  xlab("Group") + ylab("Devaluation-insensitivity")+#ggtitle("Habitual errors")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=grps_cols_order1)+ 
  theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
dg_point
rm(d, dg)

if (write_it==TRUE){
  ggsave(file="deval_accuracy_jitterbar.svg", plot=dg_point, width=4, height=6)
  ggsave(file="accuracy_after_deval2.svg", plot=rl_acc_devalpoint, width=8, height=6)
}

## ~ Individual: Outcome insensitivity - scatter individuals ####
accuracy_acc <- dplyr::select(rldeval_acc, Group, Subject, Congruency, WM_Load, ACC)
accuracy_se <- dplyr::select(rldeval_acc, Group, Subject, Congruency, WM_Load, se)
accuracy_scatter_acc <- spread(accuracy_acc, key=Congruency, value = ACC, sep="_acc_")
accuracy_scatter_se <- spread(accuracy_se, key=Congruency, value = se, sep="_se_")
accuracy_scatter <- left_join(accuracy_scatter_acc, accuracy_scatter_se)

accuracy_scatter$Group2 <- factor(accuracy_scatter$Group, levels = c("CT","EOT","COT"))

a <- ggplot(filter(accuracy_scatter, WM_Load=="No Load"), aes(x=Congruency_acc_Valued, y=Congruency_acc_Devalued)) + 
  geom_point(aes(shape=Group, color=Group2), size = 5, alpha=.7) +
  geom_errorbar(aes(ymin=Congruency_acc_Devalued-Congruency_se_Devalued, ymax=Congruency_acc_Devalued+Congruency_se_Devalued)) +
  geom_errorbarh(aes(xmin=Congruency_acc_Valued-Congruency_se_Valued, xmax=Congruency_acc_Valued+Congruency_se_Valued)) +
  geom_abline() +
  xlab("Valued pair") + ylab("Devalued pair")+ggtitle("Accuracy after devaluation")+
  facet_grid(.~Group2) +
  scale_x_continuous(limits = c(-0.2,1.2)) +
  scale_y_continuous(limits = c(-0.2,1.2)) +
  theme_bw(base_size=20)+ scale_color_manual(values=grps_cols_order1)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=16, colour="Black"), axis.title.x = element_text(size=20, colour="Black")) +
  theme(axis.text.y= element_text(size=20, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=22, colour="Black"))  +
  theme(aspect.ratio=1)
a

# and with the dual task manipulation
awm <- ggplot(filter(accuracy_scatter, WM_Load=="WM Load"), aes(x=Congruency_acc_Valued, y=Congruency_acc_Devalued)) + 
  geom_point(aes(shape=Group, color=Group2), size = 5, alpha=.7) +
  geom_errorbar(aes(ymin=Congruency_acc_Devalued-Congruency_se_Devalued, ymax=Congruency_acc_Devalued+Congruency_se_Devalued)) +
  geom_errorbarh(aes(xmin=Congruency_acc_Valued-Congruency_se_Valued, xmax=Congruency_acc_Valued+Congruency_se_Valued)) +
  geom_abline() +
  xlab("Valued pair") + ylab("Devalued pair")+ggtitle("Accuracy after devaluation - dual task")+
  facet_grid(.~Group2) +
  scale_x_continuous(limits = c(-0.2,1.2)) +
  scale_y_continuous(limits = c(-0.2,1.2)) +
  theme_bw(base_size=20)+ scale_color_manual(values=grps_cols_order1)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=16, colour="Black"), axis.title.x = element_text(size=20, colour="Black")) +
  theme(axis.text.y= element_text(size=20, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=22, colour="Black"))  +
  theme(aspect.ratio=1)
awm

## ~ RT Outcome insensitivity - dot ####
rldeval_rt_g <-summarySEwithin(filter(RL_alldeval, Missed==0, OptimalChoice==1), measurevar="RT", withinvars=c("Group","WM_Load", "Congruency"),conf.interval = .95)
rldeval_rtg <- filter(RL_alldeval, Missed==0, OptimalChoice==1) %>% group_by(Group, WM_Load,Congruency) %>% dplyr::summarize(meanRT=mean(RT, na.rm=TRUE), seRT=(sd(RT, na.rm=TRUE))/sqrt(42)) %>% left_join(rldeval_rt_g)
rm(rldeval_rt_g)
rldeval_rtg$Group2 <- factor(rldeval_rtg$Group, levels = c("CT","EOT","COT"))
rldeval_rtg$Congruency2 <- factor(rldeval_rtg$Congruency, levels = c("Valued","Devalued"))

rl_RT_devalpoint <- ggplot(filter(rldeval_rtg, WM_Load == "No Load"), aes(x=Group2, y=meanRT, group=Congruency)) + scale_y_continuous(limits = c(0.3, 0.8))+
  geom_errorbar(aes(ymin=meanRT-seRT, ymax=meanRT+seRT), colour="black", width=0, position=position_dodge(.4), size = 1) +
  geom_point(aes(fill=Group2, shape=Congruency, size=Congruency),color="Black", stroke =1, position=position_dodge(.4), stat="identity") +
  xlab("Group") + ylab("Response time (s)")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_order1)+
  theme_bw(base_size=20)+ scale_size_manual(values=c(6,8))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"), axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) #+
rl_RT_devalpoint

## ~ RT Outcome insensitivity - bar plus jitter ####
rl_reactiontime_betweenerrorbars <-summarySE(filter(RL_alldeval, OptimalChoice ==1, WM_correct != 0, Missed==0), measurevar="RT", groupvars = c("Group", "Subject", "Deval_Trial", "WM_Trial", "Day_Block", "Congruency"), conf.interval = .95)
rl_rt_deval <- filter(rl_reactiontime_betweenerrorbars, (Deval_Trial == 1 & WM_Trial ==0))
dd <- filter(rl_rt_deval, WM_Trial == 0)
d_congruent <- filter(dd, Congruency == "Valued")
d_congruent$rt_con <- d_congruent$RT
d_congruent <- dplyr::select(d_congruent, Group, Subject, rt_con)
d_incongruent <- filter(dd, Congruency == "Devalued")
d_incongruent$rt_incon <- d_incongruent$RT
d_incongruent <- dplyr::select(d_incongruent, Group, Subject, rt_incon)
d <- inner_join(d_congruent, d_incongruent, by=c("Subject", "Group"))

d$RT_diff <- d$rt_incon - d$rt_con
dg <- group_by(d, Group) %>% dplyr::summarize(rt_diff=mean(RT_diff, na.rm=TRUE))

dg$Group2 <- factor(dg$Group, levels = c("CT","EOT","COT"))
d$Group2 <- factor(d$Group, levels = c("CT","EOT","COT"))

dg_point <- ggplot() + scale_y_continuous(limits = c(-0.2,.4))+
  geom_bar(data=dg, aes(x=Group2, y=rt_diff, group=Group2, fill=Group2, width=0.4), position=position_dodge(1), size = 2, stat="identity") +
  geom_point(data = d, aes(x=Group2, y= RT_diff, color=Group2),stroke =1, color="Grey30", position=position_jitter(.25), size = 1, stat="identity") +
  xlab("Group") + ylab("RT Devalued - RT valued")+
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
dg_point

if (write_it==TRUE){
  ggsave(file="deval_RT_jitterbar_woutlier.svg", plot=dg_point, width=4, height=6)
  ggsave(file="RT_after_deval.svg", plot=rl_RT_devalpoint, width=8, height=6)
}

## ~ Individual: RT Outcome insensitivity - scatter individuals ####
rldeval_reactiontime_betweenerrorbars <-summarySE(filter(RL_alldeval, OptimalChoice ==1, WM_correct != 0, Missed==0), measurevar="RT", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
rldeval_reactiontime_withinerrorbars <-summarySEwithin(filter(RL_alldeval, OptimalChoice ==1,WM_correct >0, Missed==0), measurevar="RT", betweenvars = c("Group", "Subject"), withinvars=c("WM_Load", "Congruency"),conf.interval = .95)
rldeval_RT <- filter(RL_alldeval, OptimalChoice ==1,WM_correct >0, Missed==0) %>% group_by(Group, Subject, WM_Load,Congruency) %>% dplyr::summarize(meanRT=mean(RT, na.rm=TRUE), medianRT = median(RT, na.rm=TRUE)) %>% left_join(rldeval_reactiontime_withinerrorbars)

RT_rt <- dplyr::select(rldeval_RT, Group, Subject, Congruency, WM_Load, RT)
RT_se <- dplyr::select(rldeval_RT, Group, Subject, Congruency, WM_Load, se)
RT_scatter_rt <- spread(RT_rt, key=Congruency, value = RT, sep="_rt_")
RT_scatter_se <- spread(RT_se, key=Congruency, value = se, sep="_se_")
RT_scatter <- left_join(RT_scatter_rt, RT_scatter_se)


RT_scatter$Group2 <- factor(RT_scatter$Group, levels = c("CT","EOT","COT"))
r <- ggplot(filter(RT_scatter, WM_Load=="No Load"), aes(x=Congruency_rt_Valued, y=Congruency_rt_Devalued)) + 
  geom_point(aes(shape=Group, color=Group2), size = 5, alpha=.7) +
  geom_errorbar(aes(ymin=Congruency_rt_Devalued-Congruency_se_Devalued, ymax=Congruency_rt_Devalued+Congruency_se_Devalued)) +
  geom_errorbarh(aes(xmin=Congruency_rt_Valued-Congruency_se_Valued, xmax=Congruency_rt_Valued+Congruency_se_Valued)) +
  geom_abline() +
  facet_grid(WM_Load~Group2) +
  xlab("RT Valued pair") + ylab("RT Devalued pair")+ggtitle("RT after devaluation")+
  scale_x_continuous(limits = c(0.3,.9)) +
  scale_y_continuous(limits = c(0.3,.9)) +
  theme_bw(base_size=20)+ scale_color_manual(values=grps_cols_order1)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=16, colour="Black"), axis.title.x = element_text(size=20, colour="Black")) +
  theme(axis.text.y= element_text(size=20, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=22, colour="Black"))  +
  theme(aspect.ratio=1)
r

# and with the dual task manipulation
RT_scatter$Group2 <- factor(RT_scatter$Group, levels = c("CT","EOT","COT"))
rwm <- ggplot(filter(RT_scatter, WM_Load=="WM Load"), aes(x=Congruency_rt_Valued, y=Congruency_rt_Devalued)) + 
  geom_point(aes(shape=Group, color=Group2), size = 5, alpha=.7) +
  geom_errorbar(aes(ymin=Congruency_rt_Devalued-Congruency_se_Devalued, ymax=Congruency_rt_Devalued+Congruency_se_Devalued)) +
  geom_errorbarh(aes(xmin=Congruency_rt_Valued-Congruency_se_Valued, xmax=Congruency_rt_Valued+Congruency_se_Valued)) +
  geom_abline() +
  facet_grid(.~Group2) +
  xlab("RT Valued pair") + ylab("RT Devalued pair")+ggtitle("RT after devaluation - dual task")+
  scale_x_continuous(limits = c(0.38,1.2)) +
  scale_y_continuous(limits = c(0.38,1.2)) +
  theme_bw(base_size=20)+ scale_color_manual(values=grps_cols_order1)+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x=element_text(size=16, colour="Black"), axis.title.x = element_text(size=20, colour="Black")) +
  theme(axis.text.y= element_text(size=20, colour="Black"),axis.title.y= element_text(size=20, colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(size=22, colour="Black"))  +
  theme(aspect.ratio=1)
rwm

if (write_it==TRUE){
  ggsave(file="Individual_scatter_accuracy.svg", plot=a, width=20, height=7)
  ggsave(file="Individual_scatter_rt.svg", plot=r, width=20, height=7)
  ggsave(file="Individual_scatter_WMaccuracy.svg", plot=awm, width=20, height=7)
  ggsave(file="Individual_scatter_WMrt.svg", plot=rwm, width=20, height=7)
}

## ~ Values for the MBMF task ####
setwd(datadir)
cols <- c("Black", "grey")
lines <- c("solid", "longdash")

twostep_vals_line <- ggplot(twostep_vals, aes(x=Trial, y=Reward, Group=Stimulus)) + scale_y_continuous(limits = c(0, 10), breaks=pretty_breaks())+
  geom_line(aes(x=Trial, y=Reward, linetype = Stimulus), size = 1.5) +
  #geom_point(aes(shape = Congruency, fill=Congruency, group=Congruency), stroke=1, position=position_dodge(.5), color = "Black", size = 4) +
  xlab("Trial number") + ylab("Reward")+#ggtitle("Fluctuating rewards throughout the two-stage task")+
  scale_color_manual(values=cols) + scale_shape_manual(values=c(21,21)) + 
  scale_fill_manual(values=cols)+scale_linetype_manual(values=lines)+
  theme_bw(base_size=18)+ 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(colour="Black"))  
twostep_vals_line

if (write_it==TRUE){
  ggsave(file="mbmf_payoff_bw2.svg", plot=twostep_vals_line, width=9, height=6)
}

## ~ Scatterplot w and outcome-insensitivity ####
LPS$Group2 <- factor(LPS$Group, levels = c("EOT","COT","CT"))
scatter_deval_acc <- ggplot(LPS, aes(x=w, y=deval_diff_acc, color=factor(Group2))) + scale_x_continuous(limits=c(0.4, 1)) + scale_y_continuous(limits=c(-.5, 1)) +
  geom_point(aes(fill=factor(Group2), shape=factor(Group2), color=factor(Group2)), alpha=.8, stroke=.1, size=5) + 
  geom_smooth(method="lm", se=FALSE, fullrange=FALSE, level=0.95, size=1.5) +
  xlab("w (MF/MB balance)") + ylab("Devaluation-insensitivity")+
  theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  scale_color_manual(values=grps_cols_order2) + scale_shape_manual(values=c(21,22,24)) + scale_fill_manual(values=grps_cols_order2)+
  theme(legend.position = "right", legend.title = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(colour="Black"))  
scatter_deval_acc

if (write_it==TRUE){
  ggsave(file="scatter_w_outcome_insens_circleoutline.svg", plot=scatter_deval_acc, width=10, height=6)
}

## ~ Working memory load ####

## . . What load used per group ####
wm_load_sizes <- summarySE(LPS, measurevar="WM_used", groupvars = c("Group"), na.rm="TRUE")
wm_load_sizes$Group2 <- factor(wm_load_sizes$Group, levels = c("CT","EOT","COT"))
LPS$Group2 <- factor(LPS$Group, levels = c("CT","EOT","COT"))

bar_wmload <- ggplot(wm_load_sizes, aes(x=Group2, y=WM_used)) + scale_y_continuous(limits = c(0,8))+
  geom_bar(aes(x=Group2, y=WM_used, colour=Group2), fill="white", width=.4, size=2, , stat="identity") +
  geom_errorbar(aes(ymin=WM_used-se, ymax=WM_used+se), colour="black", width=0, position=position_dodge(width=0), size =2) +
  geom_point(data=LPS, aes(x=Group2, y= WM_used, color=Group2),stroke =1, alpha=.5, position=position_jitter(.15), size = 2, stat="identity") +
  ylab("Dual task: letter string length used")+
  theme_bw(base_size=20)+ scale_colour_manual(values=grps_cols_order1)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_wmload

#boxplot with overlying data (hence make boxplot outlier invisible)
boxplot_wmload <- ggplot(LPS, aes(x=Group2, y=WM_used)) + scale_y_continuous(limits = c(3,8))+
  geom_boxplot(aes(x=Group2, y=WM_used, colour=Group2), outlier.alpha=0,fill="white", width=.4, size=2) +
  geom_point(aes(x=Group2, y= WM_used, color=Group2),stroke =1, alpha=.5, position=position_jitter(.2), size = 2, stat="identity") +
  ylab("Dual task: letter string length used")+
  theme_bw(base_size=20)+ scale_colour_manual(values=grps_cols_order1)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
boxplot_wmload


## . . proportion correct per group ####
p <- filter(RL_deval_wm, WM_response>0, WM_correct!=666) %>% group_by(Group, Subject) %>% dplyr::summarize(MeanWMCorrect=mean(WM_correct, na.rm=TRUE))
WM_performance <- summarySE(p, measurevar = "MeanWMCorrect", groupvars = c("Group"), na.rm="TRUE")
WM_performance$Group2 <- factor(WM_performance$Group, levels = c("CT","EOT","COT"))
p$Group2 <- factor(p$Group, levels = c("CT","EOT","COT"))
bar_wm_perform <- ggplot(WM_performance, aes(x=Group2, y=MeanWMCorrect)) + scale_y_continuous(limits = c(0,1))+
  geom_bar(aes(x=Group2, y=MeanWMCorrect, colour=Group2), fill="white", width=0.4, size=2, stat="identity") +
  geom_errorbar(aes(ymin=MeanWMCorrect-se, ymax=MeanWMCorrect+se), colour="black", width=0, position=position_dodge(width=0), size =2) +
  geom_point(data=p, aes(x=Group2, y= MeanWMCorrect, color=Group2),stroke =1, alpha=.5, position=position_jitter(.15), size = 2, stat="identity") +
  ylab("Dual task: proportion correct")+#ggtitle("logistic regression on accuracy")+
  theme_bw(base_size=20)+ scale_colour_manual(values=grps_cols_order1)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
bar_wm_perform


boxplot_wm_perform <- ggplot(p, aes(x=Group2, y=MeanWMCorrect)) + scale_y_continuous(limits = c(0,1))+
  geom_boxplot(aes(x=Group2, y=MeanWMCorrect, colour=Group2), outlier.alpha=0,fill="white", width=.4, size=2) +
  geom_point(aes(x=Group2, y= MeanWMCorrect, color=Group2),stroke =1, alpha=.5, position=position_jitter(.2), size = 2, stat="identity") +
  ylab("Dual task: proportion correct")+#ggtitle("logistic regression on accuracy")+
  theme_bw(base_size=20)+ scale_colour_manual(values=grps_cols_order1)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black", size = 14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
boxplot_wm_perform

if (write_it==TRUE){
  ggsave(file="wm_bar_dots_loads.svg", plot=bar_wmload, width=3.7, height=6)
  ggsave(file="wm_boxplot_loads.svg", plot=boxplot_wmload, width=3.7, height=6)
  ggsave(file="wm_bar_dots_performance.svg", plot=bar_wm_perform, width=4, height=6)
  ggsave(file="wm_boxplot_performance.svg", plot=boxplot_wm_perform, width=4, height=6)
  
}

## . . WML- Devaluation insensitivity ####
rldeval_accuracy_g <-summarySEwithin(filter(RL_alldeval, Missed==0), measurevar="OptimalChoice", withinvars=c("Group","WM_Load", "Congruency"),conf.interval = .95)
rldeval_accg <- filter(RL_alldeval, Missed==0) %>% group_by(Group, WM_Load,Congruency) %>% dplyr::summarize(accuracy=mean(OptimalChoice, na.rm=TRUE), SE=(sd(OptimalChoice, na.rm=TRUE))/sqrt(42)) %>% left_join(rldeval_accuracy_g)
rm(rldeval_accuracy_g)
rldeval_accg$Group2 <- factor(rldeval_accg$Group, levels = c("CT","EOT","COT"))
rldeval_accg$Congruency2 <- factor(rldeval_accg$Congruency, levels = c("Valued","Devalued"))
grps_cols_lots = c("#0066CC", "#009900", "#FF8000","white", "grey", "pink", "blue", "purple", "black", "cyan") # CT EOT COT

# rlwm_acc_devalpoint <- ggplot(rldeval_accg, aes(x=Group2, y=accuracy, group=Congruency2)) + scale_y_continuous(limits = c(0.5, 1.1))+
#   geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
#   geom_point(aes(fill=Group2, shape=Congruency2, color=Group2), stroke =3, position=position_dodge(.4), size = 6, stat="identity") +
#   geom_point(aes(fill="white", shape=Congruency2, alpha=WM_Load, color=Group2), stroke =0.01, position=position_dodge(.4), size = 6, stat="identity") +
#   xlab("Devaluation") + ylab("Proportion correct")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
#   facet_grid(.~WM_Load) +
#   scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_lots)+ scale_alpha_manual(values=c(0, 1))+
#   scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
#   theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
#   theme(legend.position = "none", legend.title = element_blank()) +
#   theme(panel.border = element_blank(), axis.line = element_line())+
#   theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
#   theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
#   theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
#   theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
#   theme(plot.title=element_text(colour="Black"))
# rlwm_acc_devalpoint

rl_accuracy_withinerrorbars$Group2 <- factor(rl_accuracy_withinerrorbars$Group, levels = c("CT","EOT","COT"))
rl_accuracy_withinerrorbars$Congruency2 <- factor(rl_accuracy_withinerrorbars$Congruency, levels = c("Valued","Devalued"))

alter_rlwm_acc_devalpoint <- ggplot(rldeval_accg, aes(x=Congruency2, y=accuracy, group=Group2)) + #scale_y_continuous(limits = c(0.45, 1.1))+
  #geom_point(data=filter(rl_accuracy_betweenerrorbars, Deval_Trial==1), aes(x=Congruency2, y= OptimalChoice, color=Group2, group=Group2), stroke =1, alpha=.5, position=position_jitter(.2), size = 2, stat="identity") +
  geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
  geom_point(aes(fill=Group2, shape=Congruency2, color=Group2), stroke =3, position=position_dodge(.4), size = 6, stat="identity") +
  geom_point(aes(fill="white", shape=Congruency2, alpha=WM_Load, color=Group2), stroke =0.01, position=position_dodge(.4), size = 5, stat="identity") +
  geom_violin(data=rl_accuracy_betweenerrorbars, aes(x=Congruency2, y=OptimalChoice, group=Group2, color=Group2), position=position_dodge(.4))+
  #geom_violin(data=rl_accuracy_betweenerrorbars, aes(x=Group2, y=OptimalChoice, group=Congruency2, colour=Group2), position=position_dodge(.4))+
  xlab("Devaluation") + ylab("Proportion correct")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
  facet_grid(.~WM_Load) +
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_lots)+ scale_alpha_manual(values=c(0, 1))+
  scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
  theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(colour="Black"))  
alter_rlwm_acc_devalpoint

## b. more condensed versions
# rldeval_accuracy_gb <-summarySEwithin(filter(RL_alldeval, Missed==0), measurevar="OptimalChoice", withinvars=c("Group","WM_Load"),conf.interval = .95)
# rldeval_accgb <- filter(RL_alldeval, Missed==0) %>% group_by(Group, WM_Load) %>% dplyr::summarize(accuracy=mean(OptimalChoice, na.rm=TRUE), SE=(sd(OptimalChoice, na.rm=TRUE))/sqrt(42)) %>% left_join(rldeval_accuracy_gb)
# rm(rldeval_accuracy_gb)
# rldeval_accgb$Group2 <- factor(rldeval_accgb$Group, levels = c("CT","EOT","COT"))
# #rldeval_accgb$Congruency2 <- factor(rldeval_accgb$Congruency, levels = c("Valued","Devalued"))
# 
# rlwm_acc_devalpointb <- ggplot(rldeval_accgb, aes(x=WM_Load, y=accuracy, group=Group2)) + scale_y_continuous(limits = c(0.5, 1.1))+
#   geom_errorbar(aes(ymin=accuracy-SE, ymax=accuracy+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
#   geom_point(aes(fill=Group2, color=Group2), stroke =2, position=position_dodge(.4), size = 6, stat="identity") +
#   geom_point(aes(fill="white", alpha=WM_Load, color=Group2), stroke =0.01, position=position_dodge(.4), size = 6, stat="identity") +
#   xlab("Devaluation") + ylab("Proportion correct")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
#   #facet_grid(.~Group2) +
#   scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_order1)+ scale_alpha_manual(values=c(0, 1))+
#   scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
#   theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
#   theme(legend.position = "none", legend.title = element_blank()) +
#   theme(panel.border = element_blank(), axis.line = element_line())+
#   theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
#   theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
#   theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
#   theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
#   theme(plot.title=element_text(colour="Black"))  
# rlwm_acc_devalpointb


# RT
rldeval_rt_g <-summarySEwithin(filter(RL_alldeval, Missed==0, OptimalChoice==1), measurevar="RT", withinvars=c("Group","WM_Load", "Congruency"), conf.interval = .95)
rldeval_rtg <- filter(RL_alldeval, Missed==0, OptimalChoice==1) %>% group_by(Group, WM_Load,Congruency) %>% dplyr::summarize(SE=(sd(RT))/sqrt(42), RT=mean(RT, na.rm=TRUE)) %>% left_join(rldeval_rt_g)
rm(rldeval_rt_g)
rldeval_rtg$Group2 <- factor(rldeval_rtg$Group, levels = c("CT","EOT","COT"))
rldeval_rtg$Congruency2 <- factor(rldeval_rtg$Congruency, levels = c("Valued","Devalued"))

rlwm_rt_devalpoint <- ggplot(rldeval_rtg, aes(x=Group2, y=RT, group=Congruency2)) + scale_y_continuous(limits = c(0.4, 1.1))+
  geom_errorbar(aes(ymin=RT-SE, ymax=RT+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
  geom_point(aes(fill=Group2, shape=Congruency2, color=Group2), stroke =3, position=position_dodge(.4), size = 6, stat="identity") +
  geom_point(aes(fill="white", shape=Congruency2, alpha=WM_Load, color=Group2), stroke =0.1, position=position_dodge(.4), size = 6, stat="identity") +
  xlab("Devaluation") + ylab("Response time (s)")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
  facet_grid(.~WM_Load) +
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_lots)+ scale_alpha_manual(values=c(0, 1))+
  scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
  theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
rlwm_rt_devalpoint

alter_rlwm_rt_devalpoint <- ggplot(rldeval_rtg, aes(x=Congruency2, y=RT, group=Group2)) + scale_y_continuous(limits = c(0.3, 1.1))+
  geom_errorbar(aes(ymin=RT-SE, ymax=RT+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
  geom_point(aes(fill=Group2, shape=Congruency2, color=Group2), stroke =3, position=position_dodge(.4), size = 5, stat="identity") +
  geom_point(aes(fill="white", shape=Congruency2, alpha=WM_Load, color=Group2), stroke =0.1, position=position_dodge(.4), size = 5, stat="identity") +
  xlab("Devaluation") + ylab("Response time (s)")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
  facet_grid(.~WM_Load) +
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_lots)+ scale_alpha_manual(values=c(0, 1))+
  scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
  theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
  theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
  theme(plot.title=element_text(colour="Black"))  
alter_rlwm_rt_devalpoint

# rldeval_rt_overg <-summarySEwithin(filter(RL_alldeval, Missed==0, OptimalChoice==1), measurevar="RT", withinvars=c("WM_Load", "Congruency"), conf.interval = .95)
# rldeval_rtoverg <- filter(RL_alldeval, Missed==0, OptimalChoice==1) %>% group_by(WM_Load,Congruency) %>% dplyr::summarize(SE=(sd(RT))/sqrt(42), RT=mean(RT, na.rm=TRUE)) %>% left_join(rldeval_rt_overg)
# rm(rldeval_rt_overg)
# #rldeval_rtoverg$Group2 <- factor(rldeval_rtoverg$Group, levels = c("CT","EOT","COT"))
# rldeval_rtoverg$Congruency2 <- factor(rldeval_rtoverg$Congruency, levels = c("Valued","Devalued"))
# 
# alter_rlwm_rt_devalpoint <- ggplot(rldeval_rtoverg, aes(x=WM_Load, y=RT, group=Congruency2)) + scale_y_continuous(limits = c(0.4, 1.1))+
#   geom_errorbar(aes(ymin=RT-SE, ymax=RT+SE), colour="black", width=0, position=position_dodge(.4), size = 1) +
#   geom_point(aes(fill=Congruency2, shape=Congruency2, color=Congruency2), stroke =1, position=position_dodge(.4), size = 6, stat="identity") +
#   #geom_point(aes(fill="white", shape=Congruency2, alpha=WM_Load, color=Group2), stroke =0.1, position=position_dodge(.4), size = 6, stat="identity") +
#   xlab("Devaluation") + ylab("Response time (s)")+#ggtitle("Accuracy: WM Load effects on devaluation-insensitivity")+
#   #facet_grid(.~Group2) +
#   scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + scale_fill_manual(values=grps_cols_order1)+ scale_alpha_manual(values=c(0, 1))+
#   scale_x_discrete(labels=c("Criterion train" = "CT", "Overtrain on other pairs" = "OT other", "Overtrain on same pairs" = "OT same"))+
#   theme_bw(base_size = 20)  + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
#   theme(legend.position = "bottom", legend.title = element_blank()) +
#   theme(panel.border = element_blank(), axis.line = element_line())+
#   theme(axis.text.x=element_text(colour="Black", size=14), axis.title.x = element_blank()) +
#   theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
#   theme(axis.line.y = element_line(size=0.75,colour="Black"),axis.line.x = element_line(size=0.75,colour="Black")) +
#   theme(axis.ticks.y = element_line(size=1,colour="Black"),axis.ticks.x = element_line(size=1,colour="Black")) +
#   theme(plot.title=element_text(colour="Black"))  
# alter_rlwm_rt_devalpoint


if (write_it==TRUE){
  ggsave(file="alt_wm_rt.svg", plot=alter_rlwm_rt_devalpoint, width=10, height=6)
  ggsave(file="alt_wm_acc.svg", plot=alter_rlwm_acc_devalpoint, width=10, height=6)
}

## . . WML Outcome insensitivity - bar plus jitter ####
rldeval_accuracy_betweenerrorbars <-summarySE(filter(RL_alldeval, WM_correct = 1, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
dd <- filter(rldeval_accuracy_betweenerrorbars, WM_Load == "WM Load")
d_congruent <- filter(dd, Congruency == "Valued")
d_congruent$error_con <- 1-d_congruent$OptimalChoice
d_congruent <- dplyr::select(d_congruent, Group, Subject, error_con, sd)
d_incongruent <- filter(dd, Congruency == "Devalued")
d_incongruent$error_incon <- 1-d_incongruent$OptimalChoice
d_incongruent <- dplyr::select(d_incongruent, Group, Subject, error_incon, sd)
d <- inner_join(d_congruent, d_incongruent, by=c("Subject", "Group"))
d$habit_errorrate <- d$error_incon - d$error_con
dg <- group_by(d, Group) %>% dplyr::summarize(hab_err=mean(habit_errorrate, na.rm=TRUE))

dg$Group2 <- factor(dg$Group, levels = c("CT","EOT","COT"))
d$Group2 <- factor(d$Group, levels = c("CT","EOT","COT"))

dg_point <- ggplot() + #scale_y_continuous(limits = c(-0.4,1))+
  geom_bar(data=dg, aes(x=Group2, y=hab_err, group=Group2, fill=Group2), width=.4, stat="identity") +
  geom_point(data = d, aes(x=Group2, y= habit_errorrate, color=Group2),stroke =1, color="Grey30", position=position_jitter(.15), size = 1, stat="identity") +
  xlab("Group") + ylab("Accuracy difference of differences")+#ggtitle("Habitual errors")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=grps_cols_order1)+ 
  theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
dg_point
rm(d, dg)

## . . effect of WM load, bars with jitter overlaid ####
dd <- summarySE(filter(RL_alldeval, Missed==0), measurevar="OptimalChoice", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
d_congruent_NWM <- filter(dd, Congruency == "Valued", WM_Load == "No Load")
d_congruent_NWM$error_con_NWM <- 1-d_congruent_NWM$OptimalChoice
d_congruent_NWM <- dplyr::select(d_congruent_NWM, Group, Subject, error_con_NWM, sd)
d_incongruent_NWM <- filter(dd, Congruency == "Devalued", WM_Load == "No Load")
d_incongruent_NWM$error_incon_NWM <- 1-d_incongruent_NWM$OptimalChoice
d_incongruent_NWM <- dplyr::select(d_incongruent_NWM, Group, Subject, error_incon_NWM, sd)
d <- inner_join(d_congruent_NWM, d_incongruent_NWM, by=c("Subject", "Group"))
d$habit_errorrate_NWM <- d$error_incon_NWM - d$error_con_NWM

d_congruent_WM <- filter(dd, Congruency == "Valued", WM_Load == "WM Load")
d_congruent_WM$error_con_WM <- 1-d_congruent_WM$OptimalChoice
d_congruent_WM <- dplyr::select(d_congruent_WM, Group, Subject, error_con_WM, sd)
d_incongruent_WM <- filter(dd, Congruency == "Devalued", WM_Load == "WM Load")
d_incongruent_WM$error_incon_WM <- 1-d_incongruent_WM$OptimalChoice
d_incongruent_WM <- dplyr::select(d_incongruent_WM, Group, Subject, error_incon_WM, sd)
d2 <- inner_join(d_congruent_WM, d_incongruent_WM, by=c("Subject", "Group"))
d$error_incon_WM <- d2$error_incon_WM
d$error_con_WM <- d2$error_con_WM
d$habit_errorrate_WM <- d$error_incon_WM - d$error_con_WM
d$haberr_WM_NWM <- d$habit_errorrate_WM - d$habit_errorrate_NWM


dg <- group_by(d, Group) %>% dplyr::summarize(hab_err_WM=mean(haberr_WM_NWM, na.rm=TRUE))

dg$Group2 <- factor(dg$Group, levels = c("CT","EOT","COT"))
d$Group2 <- factor(d$Group, levels = c("CT","EOT","COT"))
dg_point_acc <- ggplot() + scale_y_continuous(limits = c(-0.4,1))+
  geom_bar(data=dg, aes(x=Group2, y=hab_err_WM, group=Group2, fill=Group2), width=.4, stat="identity") +
  geom_point(data = d, aes(x=Group2, y= haberr_WM_NWM, color=Group2),stroke =1, color="Grey50", position=position_jitter(.15), size = 1, stat="identity") +
  xlab("Group") + ylab("di")+#ggtitle("Habitual errors")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=grps_cols_order1)+ 
  theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
dg_point_acc
rm(d, dg, d2, dd)

## same for the RT
dd <- summarySE(filter(RL_alldeval, OptimalChoice ==1, Missed==0), measurevar="RT", groupvars = c("Group", "Subject", "WM_Load", "Congruency"), conf.interval = .95)
d_congruent_NWM <- filter(dd, Congruency == "Valued", WM_Load == "No Load")
d_congruent_NWM$rt_con_NWM <- d_congruent_NWM$RT
d_congruent_NWM <- dplyr::select(d_congruent_NWM, Group, Subject, rt_con_NWM, sd)
d_incongruent_NWM <- filter(dd, Congruency == "Devalued", WM_Load == "No Load")
d_incongruent_NWM$rt_incon_NWM <- d_incongruent_NWM$RT
d_incongruent_NWM <- dplyr::select(d_incongruent_NWM, Group, Subject, rt_incon_NWM, sd)
d <- inner_join(d_congruent_NWM, d_incongruent_NWM, by=c("Subject", "Group"))
d$habit_RTdiff_NWM <- d$rt_incon_NWM - d$rt_con_NWM

d_congruent_WM <- filter(dd, Congruency == "Valued", WM_Load == "WM Load")
d_congruent_WM$rt_con_WM <- d_congruent_WM$RT
d_congruent_WM <- dplyr::select(d_congruent_WM, Group, Subject, rt_con_WM, sd)
d_incongruent_WM <- filter(dd, Congruency == "Devalued", WM_Load == "WM Load")
d_incongruent_WM$rt_incon_WM <- d_incongruent_WM$RT
d_incongruent_WM <- dplyr::select(d_incongruent_WM, Group, Subject, rt_incon_WM, sd)
d2 <- inner_join(d_congruent_WM, d_incongruent_WM, by=c("Subject", "Group"))
d <- full_join(d, d2, by=c("Subject", "Group"))
d$habit_RTdiff_WM <- d$rt_incon_WM - d$rt_con_WM
d$habitRTdiff_WM_NWM <- d$habit_RTdiff_WM - d$habit_RTdiff_NWM


dg <- group_by(d, Group) %>% dplyr::summarize(hab_rtdiff_WM=mean(habitRTdiff_WM_NWM, na.rm=TRUE))

dg$Group2 <- factor(dg$Group, levels = c("CT","EOT","COT"))
d$Group2 <- factor(d$Group, levels = c("CT","EOT","COT"))
dg_point_rt <- ggplot() + #scale_y_continuous(limits = c(-0.4,1))+
  geom_bar(data=dg, aes(x=Group2, y=hab_rtdiff_WM, group=Group2, fill=Group2), width=.4, stat="identity") +
  geom_point(data = d, aes(x=Group2, y= habitRTdiff_WM_NWM, color=Group2),stroke =1, color="Grey50", position=position_jitter(.15), size = 1, stat="identity") +
  xlab("Group") + ylab("rt")+#ggtitle("Habitual errors")+
  scale_color_manual(values=grps_cols_order1) + scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=grps_cols_order1)+ 
  theme_bw(base_size=20)+ #scale_fill_manual(values=barcols)+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(axis.text.x=element_text(colour="Black"), axis.title.x = element_blank()) +
  theme(axis.text.y= element_text(colour="Black"),axis.title.y= element_text(colour="Black")) +
  theme(axis.line.y = element_line(colour="Black"),axis.line.x = element_line(colour="Black")) +
  theme(axis.ticks.y = element_line(colour="Black"),axis.ticks.x = element_line(colour="Black")) +
  theme(plot.title=element_text(colour="Black"))
dg_point_rt

if (write_it==TRUE){
  ggsave(file="WM_interactions_plusdots_accuracy.svg", plot=dg_point_acc, width=4, height=6)
  ggsave(file="WM_interactions_plusdots_rt.svg", plot=dg_point_rt, width=4, height=6)
}