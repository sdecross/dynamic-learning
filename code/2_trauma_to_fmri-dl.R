############################################################.
# SOURCE & LOAD & MANIPULATE DATA
############################################################.
rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", "lmerTest", "effects", "tidyverse")

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "data"
fear_data <- data %>%
  dplyr::select(-trauma) # get rid of "trauma" col; duplicate col causes trouble in merge of dfs

load("data/full-data/fMRI_dl_func_roilist.RData") # object "fMRI_dl_func_roilist"

roilist_modeling <- list() # initialize list 

for (roi in seq_along(fMRI_dl_func_roilist)){
  roilist_modeling[[roi]] <- fMRI_dl_func_roilist[[roi]] %>%
    dplyr::select(subject, trauma, CSP_1, CSP_2, CSP_3, CSP_4, CSM_1, CSM_2, CSM_3, CSM_4) %>%
    gather(key = "stim", value = "response", 3:10) %>%
    mutate(block = if_else(substr(stim, 4, 5) == "_1", 1,
                           if_else(substr(stim, 4, 5) == "_2", 2,
                                   if_else(substr(stim, 4, 5) == "_3", 3, 4)))) %>%
    mutate(stimsign = if_else(substr(stim, 1, 3) == "CSP", "CS+", "CS-")) %>%
    left_join(fear_data, by = "subject") %>%
    dplyr::select(subject, trauma, stim, response, block, stimsign, white, white_dm, inc_needs, inc_needs_dm, num_threat, abuse_dv_dur) 
# can use either white/inc_needs or white_dm/inc_needs_dm ; same results and stats for the analyses in this code
  roilist_modeling[[roi]]$block <- as.integer(roilist_modeling[[roi]]$block) # do not want block as factor; want to model linear trajectory over blocks;
                                                                             # if block is factor, it will bounce around between blocks
  roilist_modeling[[roi]]$stimsign <- factor(roilist_modeling[[roi]]$stimsign, levels = c("CS+", "CS-")) # factor stimsign; CS+ is reference
}
names(roilist_modeling) <- c("ACC_bl", "amyg_lh", "amyg_rh", "frontorbcort_bl", "frontpole_bl", "hipp_lh", "hipp_rh", "insula_bl", 
                             "parahipppost_bl", "precunpcc_bl", "subcallosal_bl", "thalamus_bl")

ACC_bl <- roilist_modeling$ACC_bl
amyg_lh <- roilist_modeling$amyg_lh
amyg_rh <- roilist_modeling$amyg_rh
frontorbcort_bl <- roilist_modeling$frontorbcort_bl
frontpole_bl <- roilist_modeling$frontpole_bl
hipp_lh <- roilist_modeling$hipp_lh
hipp_rh <- roilist_modeling$hipp_rh
insula_bl <- roilist_modeling$insula_bl
parahipppost_bl <- roilist_modeling$parahipppost_bl
precunpcc_bl <- roilist_modeling$precunpcc_bl
subcallosal_bl <- roilist_modeling$subcallosal_bl
thalamus_bl <- roilist_modeling$thalamus_bl


############################################################.
# MIXED MODELS - TRAUMA (GROUP; BINARY 0/1)              ####
############################################################.

# Does TRAUMA (GROUP) predict pattern of brain response across blocks?

# Test models: 
#   1. random intercept (RI)
#   2. RI + random slope of block (RS)
# AIC criteria used to pick the model of best fit
# Notate any findings of interest, aka effects of trauma

# ACC bl -----------------------------------------------------------------------------------
m4_ACC <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = ACC_bl)
m5_ACC <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=ACC_bl)
anova(m4_ACC, m5_ACC)
AIC(m4_ACC, m5_ACC)
anova(m5_ACC, type = "II") # m5; no

# amyg lh ----------------------------------------------------------------------------------
m4_amyg_lh <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = amyg_lh)
m5_amyg_lh <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=amyg_lh)
anova(m4_amyg_lh, m5_amyg_lh)
AIC(m4_amyg_lh, m5_amyg_lh)
exp((3638.922-3640.963)/2) # 36% chance m5 is better than m4; go with m4
anova(m4_amyg_lh, type = "II") # m4; no

# amyg rh ----------------------------------------------------------------------------------
m4_amyg_rh <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = amyg_rh)
m5_amyg_rh <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=amyg_rh)
anova(m4_amyg_rh, m5_amyg_rh)
AIC(m4_amyg_rh, m5_amyg_rh)
exp((3542.254-3544.254)/2) # 36% chance m5 is better
#Anova(m4_amyg_rh) # m4; stim*block*trauma p = 0.005
anova(m4_amyg_rh, type = "II")

# frontorbcort bl --------------------------------------------------------------------------
m4_frontorbcort <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = frontorbcort_bl)
m5_frontorbcort <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=frontorbcort_bl)
anova(m4_frontorbcort, m5_frontorbcort)
AIC(m4_frontorbcort, m5_frontorbcort)
exp((3550.232-3551.565)/2) # 51% chance m5 is better
anova(m5_frontorbcort, type = "II") # m5; stim*block*trauma p = 0.066

# frontpole bl -----------------------------------------------------------------------------
m4_frontpole <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = frontpole_bl)
m5_frontpole <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=frontpole_bl)
anova(m4_frontpole, m5_frontpole)
AIC(m4_frontpole, m5_frontpole)
anova(m5_frontpole, type = "II") # m5; stim*block*trauma p = 0.034

# hipp lh ----------------------------------------------------------------------------------
m4_hipp_lh <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = hipp_lh)
m5_hipp_lh <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=hipp_lh)
anova(m4_hipp_lh, m5_hipp_lh)
AIC(m4_hipp_lh, m5_hipp_lh)
anova(m5_hipp_lh, type = "II") # m5; no

# hipp rh ----------------------------------------------------------------------------------
m4_hipp_rh <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = hipp_rh)
m5_hipp_rh <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=hipp_rh)
anova(m4_hipp_rh, m5_hipp_rh)
AIC(m4_hipp_rh, m5_hipp_rh)
exp((3541.463-3542.241)/2) # 67% chance m4 better fit
anova(m4_hipp_rh, type = "II") # m4; stim*block*trauma p = 0.021

# insula bl --------------------------------------------------------------------------------
m4_insula <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = insula_bl)
m5_insula <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=insula_bl)
anova(m4_insula, m5_insula)
AIC(m4_insula, m5_insula)
exp((3479.480-3480.246)/2) # 68% chance m4 better fit
anova(m4_insula, type = "II") # m4; stim*block*trauma p = 0.086

# parahipppost bl --------------------------------------------------------------------------
m4_parahipppost <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = parahipppost_bl)
m5_parahipppost <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=parahipppost_bl)
anova(m4_parahipppost, m5_parahipppost)
AIC(m4_parahipppost, m5_parahipppost)
anova(m5_parahipppost, type = "II")  # m5; stim*block*trauma p = 0.090

# precunpcc bl -----------------------------------------------------------------------------
m4_precunpcc <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = precunpcc_bl)
m5_precunpcc <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=precunpcc_bl)
anova(m4_precunpcc, m5_precunpcc)
AIC(m4_precunpcc, m5_precunpcc)
anova(m5_precunpcc, type = "II") # m5; no

# subcallosal bl ---------------------------------------------------------------------------
m4_subcallosal <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = subcallosal_bl)
m5_subcallosal <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=subcallosal_bl)
anova(m4_subcallosal, m5_subcallosal)
AIC(m4_subcallosal, m5_subcallosal)
exp((4229.695-4230.225)/2) # 77% chance m4 is better 
anova(m4_subcallosal, type = "II") # no

# thalamus bl ------------------------------------------------------------------------------
m4_thalamus <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = thalamus_bl)
m5_thalamus <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=thalamus_bl)
anova(m4_thalamus, m5_thalamus)
AIC(m4_thalamus, m5_thalamus)
exp((3521.270-3523.914)/2) # 26% chance m4 is better
anova(m5_thalamus, type = "II") # m5; no






############################################################
# SUMMARY
############################################################

#------------------------------------------------------------------------------
# TRAUMA (GROUP)
#                              stim*block*trauma        model       stim*trauma
#
#   Right amygdala                 p = 0.005*          m4 (RI)
#   Bl insula                      p = 0.086           m4 (RI)
#   Right hippocampus              p = 0.021*          m4 (RI)
#   Bl front pole                  p = 0.034*          m5 (RI+RS)
#   BL front orb cort              p = 0.066           m5 (RI+RS)
#   Bl parahipppost                p = 0.090           m5 (RI+RS)


# Whether you demean white and inc_needs or not, same results and stats here in this script (did not demean here)



