###################################################################.
# Does trauma predict gPPI?
# gPPI right amyg fcMRI with: bl hipp, bl ACC, bl PCC, bl PHG
###################################################################.

rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", 
       "lmerTest", "effects", "MASS", "mgcv", "splines", "sjPlot", "gvlma", "mediation", "ggstatsplot",
       "tidyverse", "broom")


hipp_rh_amyg <- read.table("data/roi-extraction-files/gPPI/gPPI_Final_Hippocampus_bl_Harvard-Oxford_thr50_binarized_gPPI_Amygdala_rh_Harvard-Oxford_thr50_binarized_FearLearning_p05.feat_CSPvCSM_zstat1.txt", header=FALSE)       
hipp_rh_amyg <- hipp_rh_amyg %>%
  dplyr::select(subject = V1, gPPI = V2) # gPPI = CSP vs. CSM
hipp_rh_amyg$subject <- factor(hipp_rh_amyg$subject)

acc_rh_amyg <- read.table("data/roi-extraction-files/gPPI/gPPI_Final_ACC_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_gPPI_Amygdala_rh_Harvard-Oxford_thr50_binarized_FearLearning_p05.feat_CSPvCSM_zstat1.txt", header=FALSE)       
acc_rh_amyg <- acc_rh_amyg %>%
  dplyr::select(subject = V1, gPPI = V2) # gPPI = CSP vs. CSM
acc_rh_amyg$subject <- factor(acc_rh_amyg$subject)

pcc_rh_amyg <- read.table("data/roi-extraction-files/gPPI/gPPI_Final_PrecuneusPCC_VentralBlob_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_gPPI_Amygdala_rh_Harvard-Oxford_thr50_binarized_FearLearning_p05.feat_CSPvCSM_zstat1.txt", header=FALSE)       
pcc_rh_amyg <- pcc_rh_amyg %>%
  dplyr::select(subject = V1, gPPI = V2) # gPPI = CSP vs. CSM
pcc_rh_amyg$subject <- factor(pcc_rh_amyg$subject)

vvs_rh_amyg <- read.table("data/roi-extraction-files/gPPI/gPPI_Final_Parahipp_post_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_gPPI_Amygdala_rh_Harvard-Oxford_thr50_binarized_FearLearning_p05.feat_CSPvCSM_zstat1.txt", header=FALSE)       
vvs_rh_amyg <- vvs_rh_amyg %>%
  dplyr::select(subject = V1, gPPI = V2) # gPPI = CSP vs. CSM
vvs_rh_amyg$subject <- factor(vvs_rh_amyg$subject)

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "mt_fearcond_data"
fear_data <- data 

hipp_rh_amyg <- hipp_rh_amyg %>%
  left_join(fear_data, by = "subject")

acc_rh_amyg <- acc_rh_amyg %>%
  left_join(fear_data, by = "subject")

pcc_rh_amyg <- pcc_rh_amyg %>%
  left_join(fear_data, by = "subject")

vvs_rh_amyg <- vvs_rh_amyg %>%
  left_join(fear_data, by = "subject")


###################################################################. 
# TRAUMA TO GPPI --------------------------------------------------
###################################################################.

#------------------------------------------------------------------.
# VISUALIZING OUTCOME VARS 
hist(hipp_rh_amyg$gPPI, breaks = 20, xlab = "gPPI", main = "Histogram hipp_rh_amyg") 
hist(acc_rh_amyg$gPPI, breaks = 20, xlab = "gPPI", main = "Histogram acc_rh_amyg") 
hist(pcc_rh_amyg$gPPI, breaks = 20, xlab = "gPPI", main = "Histogram pcc_rh_amyg") 
hist(vvs_rh_amyg$gPPI, breaks = 20, xlab = "gPPI", main = "Histogram vvs_rh_amyg") 
# normal, use normal lms not negbin

#------------------------------------------------------------------.
# REGRESSION - TRAUMA GROUP BINARY OUTCOME VAR 
# Does trauma predict connectivity? 
# Already know this; it's extracted from the between-group map...
# do not use these for anything, just for conceptual reminder

#fithipp <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg)
#summary(fithipp) # yes p < 0.001

#fitacc <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg)
#summary(fitacc) # yes p = 0.04

#fitpcc <- lm(gPPI ~ trauma + white + inc_needs, data = pcc_rh_amyg)
#summary(fitpcc) # yes p = 0.002

#fitvvs <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg)
#summary(fitvvs) # yes p < 0.001




