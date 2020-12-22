############################################################.
# SOURCE & LOAD & MANIPULATE DATA
############################################################.
rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", "lmerTest", "effects", 
       "tidyverse", "plyr", "ggdark")

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "data"
fear_data <- data %>%
  dplyr::select(-trauma) # get rid of "trauma" col; duplicate col causes trouble in merge of dfs

load("data/full-data/fMRI_dl_func_roilist_US.RData") # object "fMRI_dl_func_roilist_US"

roilist_modeling <- list() # initialize list 

for (roi in seq_along(fMRI_dl_func_roilist_US)){
  roilist_modeling[[roi]] <- fMRI_dl_func_roilist_US[[roi]] %>%
    dplyr::select(subject, trauma, CSPR_1, CSPR_2, CSPR_3, CSPR_4, CSPN_1, CSPN_2, CSPN_3, CSPN_4) %>%
    gather(key = "stim", value = "response", 3:10) %>%
    mutate(block = if_else(substr(stim, 5, 6) == "_1", 1,
                           if_else(substr(stim, 5, 6) == "_2", 2,
                                   if_else(substr(stim, 5, 6) == "_3", 3, 4)))) %>%
    mutate(stimsign = if_else(substr(stim, 1, 4) == "CSPR", "CS+R", "CS+")) %>%
    left_join(fear_data, by = "subject") %>%
    dplyr::select(subject, trauma, stim, response, block, stimsign, white, white_dm, inc_needs, inc_needs_dm, num_threat, abuse_dv_dur) 
# can use either white/inc_needs or white_dm/inc_needs_dm ; same results and stats for the analyses in this code
  roilist_modeling[[roi]]$block <- as.integer(roilist_modeling[[roi]]$block) # do not want block as factor; want to model linear trajectory over blocks;
                                                                             # if block is factor, it will bounce around between blocks
  roilist_modeling[[roi]]$stimsign <- factor(roilist_modeling[[roi]]$stimsign, levels = c("CS+R", "CS+")) # factor stimsign; CS+R is reference
}
names(roilist_modeling) <- c("amyg_lh", "amyg_rh")


amyg_lh_US <- roilist_modeling$amyg_lh
amyg_rh_US <- roilist_modeling$amyg_rh

amyg_rh_US_1 <- amyg_rh_US %>%
  dplyr::filter(stim == "CSPR_1")

amyg_rh_US_2 <- amyg_rh_US %>%
  dplyr::filter(stim == "CSPR_2")

amyg_rh_US_3 <- amyg_rh_US %>%
  dplyr::filter(stim == "CSPR_3")

amyg_rh_US_4 <- amyg_rh_US %>%
  dplyr::filter(stim == "CSPR_4")


# interleaved histograms

mu <- ddply(amyg_rh_US_1, "trauma", summarise, grp.mean=mean(response))
ggplot(amyg_rh_US_1, aes(x=response, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("amyg rh - 1")


mu <- ddply(amyg_rh_US_2, "trauma", summarise, grp.mean=mean(response))
ggplot(amyg_rh_US_2, aes(x=response, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("amyg rh - 2")


mu <- ddply(amyg_rh_US_3, "trauma", summarise, grp.mean=mean(response))
ggplot(amyg_rh_US_3, aes(x=response, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("amyg rh - 3")

mu <- ddply(amyg_rh_US_4, "trauma", summarise, grp.mean=mean(response))
ggplot(amyg_rh_US_4, aes(x=response, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("amyg rh - 4")









############################################################.
# MIXED MODELS - TRAUMA (GROUP; BINARY 0/1)             ####
############################################################.

# Does TRAUMA (GROUP) predict pattern of brain response across blocks?

# Test models: 
#   1. random intercept (RI)
#   2. RI + random slope of block (RS)
# AIC criteria used to pick the model of best fit
# Notate any findings of interest, aka effects of trauma

# amyg lh ----------------------------------------------------------------------------------
m4_amyg_lh_US <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = amyg_lh_US)
m5_amyg_lh_US <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=amyg_lh_US)
anova(m4_amyg_lh_US, m5_amyg_lh_US)
AIC(m4_amyg_lh_US, m5_amyg_lh_US)
anova(m5_amyg_lh_US, type = "II") # m5; NO INTERACTIONS WITH TRAUMA

# amyg rh ----------------------------------------------------------------------------------
m4_amyg_rh_US <- lmer(response ~ 1 + (1|subject) + stimsign*block*trauma + white + inc_needs, data = amyg_rh_US)
m5_amyg_rh_US <- lmer(response ~ 1 + (1 + block|subject) + stimsign*block*trauma + white + inc_needs, data=amyg_rh_US)
anova(m4_amyg_rh_US, m5_amyg_rh_US)
AIC(m4_amyg_rh_US, m5_amyg_rh_US)
anova(m5_amyg_rh_US, type = "II") # m5; TREND FOR 3-WAY INTERACTION stimsign x block x trauma p = 0.08


#Type II Analysis of Variance Table with Satterthwaite's method
#                      Sum Sq   Mean Sq   NumDF  DenDF    F value   Pr(>F) 

#stimsign:block:trauma 2.7699   2.7699     1      877.98  3.0436    0.081404 . 




# Whether you demean white and inc_needs or not, same results and stats here in this script (did not demean here)



