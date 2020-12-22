############################################################.
# SOURCE & LOAD & MANIPULATE DATA
############################################################.
rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", 
       "lmerTest", "effects", "MASS", "mgcv", "splines", "sjPlot", "plyr","tidyverse", "ggdark", "sjPlot", "sjmisc")

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "data"
fear_data <- data %>%
  dplyr::select(-trauma) # get rid of "trauma" col; duplicate col causes trouble in merge of dfs

load("data/full-data/fMRI_dl_func_roilist.RData") # object "fMRI_dl_func_roilist"

roilist_regressions <- list() # initialize list

for (roi in seq_along(fMRI_dl_func_roilist)){
  roilist_regressions[[roi]] <- fMRI_dl_func_roilist[[roi]] %>%
    left_join(fear_data, by = "subject") %>%
    dplyr::select(subject, trauma, white, white_dm, inc_needs, inc_needs_dm, 
                  CSP_1, CSP_2, CSP_3, CSP_4, 
                  CSM_1, CSM_2, CSM_3, CSM_4,
                  CSPM_1, CSPM_2, CSPM_3, CSPM_4,
                  cdi, scared_panic, scared_gad, external, ptsd,
                  cdi_t2, scared_panic_t2, scared_gad_t2, external_t2, ptsd_t2) %>%
    mutate(CSMP_4 = CSM_4 - CSP_4) %>% # adding the MP so stats are consistent with safety plots
    mutate(CSMP_3 = CSM_3 - CSP_3) %>%
    mutate(CSMP_2 = CSM_2 - CSP_2) %>%
    mutate(CSMP_1 = CSM_1 - CSP_1) %>%
    mutate(CSPM_slope_linear = NA) %>%
    mutate(CSMP_slope_linear = NA)
}

names(roilist_regressions) <- c("ACC_bl", "amyg_lh", "amyg_rh", "frontorbcort_bl", "frontpole_bl", "hipp_lh", "hipp_rh", "insula_bl", "parahipppost_bl", 
                                "precunpcc_bl", "subcallosal_bl", "thalamus_bl")

amyg_rh <- roilist_regressions$amyg_rh
insula <- roilist_regressions$insula_bl
hipp_rh <- roilist_regressions$hipp_rh
frontpole <- roilist_regressions$frontpole_bl
frontorbcort <- roilist_regressions$frontorbcort_bl
parahipppost <- roilist_regressions$parahipppost_bl

# calculate sophisticated slopes
get_slope_block <- c(1, 2, 3, 4) # want to be number not factor so it fits a linear slope

# amyg_rh
for (i in 1:147){
  get_slope_bold_PM <- c(amyg_rh[i, "CSPM_1"], amyg_rh[i, "CSPM_2"], amyg_rh[i, "CSPM_3"], amyg_rh[i, "CSPM_4"])
  get_slope_bold_MP <- c(amyg_rh[i, "CSMP_1"], amyg_rh[i, "CSMP_2"], amyg_rh[i, "CSMP_3"], amyg_rh[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  amyg_rh[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  amyg_rh[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}

# insula
for (i in 1:147){
  get_slope_bold_PM <- c(insula[i, "CSPM_1"], insula[i, "CSPM_2"], insula[i, "CSPM_3"], insula[i, "CSPM_4"])
  get_slope_bold_MP <- c(insula[i, "CSMP_1"], insula[i, "CSMP_2"], insula[i, "CSMP_3"], insula[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  insula[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  insula[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}

# hipp_rh
for (i in 1:147){
  get_slope_bold_PM <- c(hipp_rh[i, "CSPM_1"], hipp_rh[i, "CSPM_2"], hipp_rh[i, "CSPM_3"], hipp_rh[i, "CSPM_4"])
  get_slope_bold_MP <- c(hipp_rh[i, "CSMP_1"], hipp_rh[i, "CSMP_2"], hipp_rh[i, "CSMP_3"], hipp_rh[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  hipp_rh[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  hipp_rh[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}

# frontpole
for (i in 1:147){
  get_slope_bold_PM <- c(frontpole[i, "CSPM_1"], frontpole[i, "CSPM_2"], frontpole[i, "CSPM_3"], frontpole[i, "CSPM_4"])
  get_slope_bold_MP <- c(frontpole[i, "CSMP_1"], frontpole[i, "CSMP_2"], frontpole[i, "CSMP_3"], frontpole[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  frontpole[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  frontpole[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}

# frontorbcort
for (i in 1:147){
  get_slope_bold_PM <- c(frontorbcort[i, "CSPM_1"], frontorbcort[i, "CSPM_2"], frontorbcort[i, "CSPM_3"], frontorbcort[i, "CSPM_4"])
  get_slope_bold_MP <- c(frontorbcort[i, "CSMP_1"], frontorbcort[i, "CSMP_2"], frontorbcort[i, "CSMP_3"], frontorbcort[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  frontorbcort[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  frontorbcort[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}

# parahipppost
for (i in 1:147){
  get_slope_bold_PM <- c(parahipppost[i, "CSPM_1"], parahipppost[i, "CSPM_2"], parahipppost[i, "CSPM_3"], parahipppost[i, "CSPM_4"])
  get_slope_bold_MP <- c(parahipppost[i, "CSMP_1"], parahipppost[i, "CSMP_2"], parahipppost[i, "CSMP_3"], parahipppost[i, "CSMP_4"])
  get_slope_df <- data.frame(get_slope_block, get_slope_bold_PM, get_slope_bold_MP)
  get_slope_model_PM <- lm(get_slope_bold_PM ~ get_slope_block)
  parahipppost[i, "CSPM_slope_linear"] <- get_slope_model_PM$coefficients["get_slope_block"]
  get_slope_model_MP <- lm(get_slope_bold_MP ~ get_slope_block)
  parahipppost[i, "CSMP_slope_linear"] <- get_slope_model_MP$coefficients["get_slope_block"]
}



# checking whether average slopes of each regions reflect the values they should

# salience regions: control should be more negative than trauma
amyg_rh %>% group_by(trauma) %>% summarise(mean_slope = mean(CSPM_slope_linear), sd_slope = sd(CSPM_slope_linear))

# interleaved histograms
mu <- ddply(amyg_rh, "trauma", summarise, grp.mean=mean(CSPM_slope_linear))
ggplot(amyg_rh, aes(x=CSPM_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("amyg rh") #+
  #scale_color_brewer(palette="Dark2")


insula %>% group_by(trauma) %>% summarise(mean_slope = mean(CSPM_slope_linear), sd_slope = sd(CSPM_slope_linear))
# interleaved histograms
mu <- ddply(insula, "trauma", summarise, grp.mean=mean(CSPM_slope_linear))
ggplot(insula, aes(x=CSPM_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("insula")#+
  #scale_color_brewer(palette="Dark2")

# default regions: control should be more positive than trauma
hipp_rh %>% group_by(trauma) %>% summarise(mean_slope = mean(CSMP_slope_linear), sd_slope = sd(CSMP_slope_linear))
# interleaved histograms
mu <- ddply(hipp_rh, "trauma", summarise, grp.mean=mean(CSMP_slope_linear))
ggplot(hipp_rh, aes(x=CSMP_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("hipp_rh")#+
  #scale_color_brewer(palette="Dark2")

frontpole %>% group_by(trauma) %>% summarise(mean_slope = mean(CSMP_slope_linear), sd_slope = sd(CSMP_slope_linear))
# interleaved histograms
mu <- ddply(frontpole, "trauma", summarise, grp.mean=mean(CSMP_slope_linear))
ggplot(frontpole, aes(x=CSMP_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("frontpole")#+
  #scale_color_brewer(palette="Dark2")

frontorbcort %>% group_by(trauma) %>% summarise(mean_slope = mean(CSMP_slope_linear), sd_slope = sd(CSMP_slope_linear))
# interleaved histograms
mu <- ddply(frontorbcort, "trauma", summarise, grp.mean=mean(CSMP_slope_linear))
ggplot(frontorbcort, aes(x=CSMP_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top")+
  ggtitle("frontorbcort")#+
  #scale_color_brewer(palette="Dark2")

parahipppost %>% group_by(trauma) %>% summarise(mean_slope = mean(CSMP_slope_linear), sd_slope = sd(CSMP_slope_linear))
# interleaved histograms
mu <- ddply(parahipppost, "trauma", summarise, grp.mean=mean(CSMP_slope_linear))
ggplot(parahipppost, aes(x=CSMP_slope_linear, color=trauma)) +
  dark_theme_classic() +
  geom_histogram(fill="black", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=trauma),
             linetype="dashed")+
  theme(legend.position="top") +
  ggtitle("parahipppost")#+
  #scale_color_brewer(palette="Dark2")


#####################################################################################################.
# DO NOT WANT TO DEMEAN SLOPES; WILL MAKE SLOPE = 0
# MEAN "AVERAGE" SLOPE INSTEAD OF "FLAT" SLOPE/INFLECTION
# BETWEEN POSITIVE AND NEGATIVE SLOPES

# Only regress the regions that showed differences between groups on symptoms (or were trending)
# Use: rh amyg, bl insula, rh hipp, bl frontpole, bl frontorbcort, bl parahipppost
# (Slopes calculated above)

# IN FULL GROUP WILL RUN THESE MODELS:
# with lm and glm.nb to see which fits better
#
# x = CSPM slope of: amyg rh, insula bl 
#     CSMP slope of: hipp rh, front pole bl, frontorbcort bl, parahippost bl
# y = cdi, scared_panic, scared_gad, external, ptsd
# covarying for white and inc_needs (not demeaned; will manually adjust if plotting just regressions)

############################################################.
# VISUALIZE OUTCOME VARS
############################################################.

# my DVS are all metric variables with lower bounds at 0 with no upper bound
# expect negbin versions to fit better for all models, except possibly external 
hist(amyg_rh$cdi, breaks = 20, xlab = "CDI counts", main = "Histogram CDI") # right-skewed
hist(amyg_rh$scared_panic, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC") # right-skewed
hist(amyg_rh$scared_gad, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD") # right-skewed
hist(amyg_rh$external, breaks = 25, xlab = "External counts", main = "Histogram external") # possibly normal
hist(amyg_rh$ptsd, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD sev") # right-skewed

############################################################.
# REGRESSIONS WITH FMRI SLOPE ####
############################################################.

# amyg rh ---------------------------------------------------------------------------------------

# cdi
scatterplot(cdi ~ amyg_rh$CSPM_slope_linear, data = amyg_rh, id = list(n = 5)) 
amyg_cdi_nor <- lm(cdi ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
amyg_cdi_nb <- glm.nb(cdi ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
AIC(amyg_cdi_nor, amyg_cdi_nb)
summary(amyg_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ amyg_rh$CSPM_slope_linear, data = amyg_rh, id = list(n = 5)) 
amyg_scared_panic_nor <- lm(scared_panic ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
amyg_scared_panic_nb <- glm.nb(scared_panic ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_panic_nor, amyg_scared_panic_nb)
summary(amyg_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ amyg_rh$CSPM_slope_linear, data = amyg_rh, id = list(n = 5)) 
amyg_scared_gad_nor <- lm(scared_gad ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
amyg_scared_gad_nb <- glm.nb(scared_gad ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_gad_nor, amyg_scared_gad_nb)
summary(amyg_scared_gad_nb) # no

# external
scatterplot(external ~ amyg_rh$CSPM_slope_linear, data = amyg_rh, id = list(n = 5)) 
amyg_external_nor <- lm(external ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
amyg_external_nb <- glm.nb(external ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
AIC(amyg_external_nor, amyg_external_nb)
exp((1074.083-1075.194)/2) # 57% chance normal is better
summary(amyg_external_nor) # no

# ptsd
scatterplot(ptsd ~ amyg_rh$CSPM_slope_linear, data = amyg_rh, id = list(n = 5)) 
amyg_ptsd_nor <- lm(ptsd ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
amyg_ptsd_nb <- glm.nb(ptsd ~ CSPM_slope_linear + white + inc_needs, data = amyg_rh)
AIC(amyg_ptsd_nor, amyg_ptsd_nb)
summary(amyg_ptsd_nb) # 0.080

# insula bl -------------------------------------------------------------------------------------

# cdi
scatterplot(cdi ~ insula$CSPM_slope_linear, data = insula, id = list(n = 5)) 
ins_cdi_nor <- lm(cdi ~ CSPM_slope_linear + white + inc_needs, data = insula)
ins_cdi_nb <- glm.nb(cdi ~ CSPM_slope_linear + white + inc_needs, data = insula)
AIC(ins_cdi_nor, ins_cdi_nb)
summary(ins_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ insula$CSPM_slope_linear, data = insula, id = list(n = 5)) 
ins_scared_panic_nor <- lm(scared_panic ~ CSPM_slope_linear + white + inc_needs, data = insula)
ins_scared_panic_nb <- glm.nb(scared_panic ~ CSPM_slope_linear + white + inc_needs, data = insula)
AIC(ins_scared_panic_nor, ins_scared_panic_nb)
summary(ins_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ insula$CSPM_slope_linear, data = insula, id = list(n = 5)) 
ins_scared_gad_nor <- lm(scared_gad ~ CSPM_slope_linear + white + inc_needs, data = insula)
ins_scared_gad_nb <- glm.nb(scared_gad ~ CSPM_slope_linear + white + inc_needs, data = insula)
AIC(ins_scared_gad_nor, ins_scared_gad_nb)
summary(ins_scared_gad_nb) # no

# external
scatterplot(external ~ insula$CSPM_slope_linear, data = insula, id = list(n = 5)) 
ins_external_nor <- lm(external ~ CSPM_slope_linear + white + inc_needs, data = insula)
ins_external_nb <- glm.nb(external ~ CSPM_slope_linear + white + inc_needs, data = insula)
AIC(ins_external_nor, ins_external_nb)
exp((1072.884-1074.031)/2) # 56% chance normal is better
summary(ins_external_nor) # no

# ptsd
scatterplot(ptsd ~ insula$CSPM_slope_linear, data = insula, id = list(n = 5)) 
ins_ptsd_nor <- lm(ptsd ~ CSPM_slope_linear + white + inc_needs, data = insula)
ins_ptsd_nb <- glm.nb(ptsd ~ CSPM_slope_linear + white + inc_needs, data = insula)
AIC(ins_ptsd_nor, ins_ptsd_nb)
summary(ins_ptsd_nb) # no

# hipp rh ---------------------------------------------------------------------------------------
# cdi
scatterplot(cdi ~ hipp_rh$CSMP_slope_linear, data = hipp_rh, id = list(n = 5)) 
hipp_cdi_nor <- lm(cdi ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
hipp_cdi_nb <- glm.nb(cdi ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
AIC(hipp_cdi_nor, hipp_cdi_nb)
summary(hipp_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ hipp_rh$CSMP_slope_linear, data = hipp_rh, id = list(n = 5)) 
hipp_scared_panic_nor <- lm(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
hipp_scared_panic_nb <- glm.nb(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
AIC(hipp_scared_panic_nor, hipp_scared_panic_nb)
summary(hipp_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ hipp_rh$CSMP_slope_linear, data = hipp_rh, id = list(n = 5)) 
hipp_scared_gad_nor <- lm(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
hipp_scared_gad_nb <- glm.nb(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
AIC(hipp_scared_gad_nor, hipp_scared_gad_nb)
summary(hipp_scared_gad_nb) # no

# external
scatterplot(external ~ hipp_rh$CSMP_slope_linear, data = hipp_rh, id = list(n = 5)) 
hipp_external_nor <- lm(external ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
hipp_external_nb <- glm.nb(external ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
AIC(hipp_external_nor, hipp_external_nb)
exp((1074.019-1075.128)/2) # 57% chance normal is better
summary(hipp_external_nor) # no

# ptsd
scatterplot(ptsd ~ hipp_rh$CSMP_slope_linear, data = hipp_rh, id = list(n = 5)) 
hipp_ptsd_nor <- lm(ptsd ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
hipp_ptsd_nb <- glm.nb(ptsd ~ CSMP_slope_linear + white + inc_needs, data = hipp_rh)
AIC(hipp_ptsd_nor, hipp_ptsd_nb)
summary(hipp_ptsd_nb) # no

# front pole bl ---------------------------------------------------------------------------------
# cdi
scatterplot(cdi ~ frontpole$CSMP_slope_linear, data = frontpole, id = list(n = 5)) 
frontpl_cdi_nor <- lm(cdi ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
frontpl_cdi_nb <- glm.nb(cdi ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
AIC(frontpl_cdi_nor, frontpl_cdi_nb)
summary(frontpl_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ frontpole$CSMP_slope_linear, data = frontpole, id = list(n = 5)) 
frontpl_scared_panic_nor <- lm(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
frontpl_scared_panic_nb <- glm.nb(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
AIC(frontpl_scared_panic_nor, frontpl_scared_panic_nb)
summary(frontpl_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ frontpole$CSMP_slope_linear, data = frontpole, id = list(n = 5)) 
frontpl_scared_gad_nor <- lm(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
frontpl_scared_gad_nb <- glm.nb(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
AIC(frontpl_scared_gad_nor, frontpl_scared_gad_nb)
summary(frontpl_scared_gad_nb) # no

# external
scatterplot(external ~ frontpole$CSMP_slope_linear, data = frontpole, id = list(n = 5)) 
frontpl_external_nor <- lm(external ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
frontpl_external_nb <- glm.nb(external ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
AIC(frontpl_external_nor, frontpl_external_nb)
exp((1073.695-1074.884)/2) # 55% chance normal is better
summary(frontpl_external_nor) # no

# ptsd
scatterplot(ptsd ~ frontpole$CSMP_slope_linear, data = frontpole, id = list(n = 5)) 
frontpl_ptsd_nor <- lm(ptsd ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
frontpl_ptsd_nb <- glm.nb(ptsd ~ CSMP_slope_linear + white + inc_needs, data = frontpole)
AIC(frontpl_ptsd_nor, frontpl_ptsd_nb)
summary(frontpl_ptsd_nb) # no

# front orb cort bl -----------------------------------------------------------------------------
# cdi
scatterplot(cdi ~ frontorbcort$CSMP_slope_linear, data = frontorbcort, id = list(n = 5)) 
frontoc_cdi_nor <- lm(cdi ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
frontoc_cdi_nb <- glm.nb(cdi ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
AIC(frontoc_cdi_nor, frontoc_cdi_nb)
summary(frontoc_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ frontorbcort$CSMP_slope_linear, data = frontorbcort, id = list(n = 5)) 
frontoc_scared_panic_nor <- lm(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
frontoc_scared_panic_nb <- glm.nb(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
AIC(frontoc_scared_panic_nor, frontoc_scared_panic_nb)
summary(frontoc_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ frontorbcort$CSMP_slope_linear, data = frontorbcort, id = list(n = 5)) 
frontoc_scared_gad_nor <- lm(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
frontoc_scared_gad_nb <- glm.nb(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
AIC(frontoc_scared_gad_nor, frontoc_scared_gad_nb)
summary(frontoc_scared_gad_nb) # no

# external
scatterplot(external ~ frontorbcort$CSMP_slope_linear, data = frontorbcort, id = list(n = 5)) 
frontoc_external_nor <- lm(external ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
frontoc_external_nb <- glm.nb(external ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
AIC(frontoc_external_nor, frontoc_external_nb)
exp((1075.044-1076.134)/2) # 58% chance normal is better
summary(frontoc_external_nor) # no

# ptsd
scatterplot(ptsd ~ frontorbcort$CSMP_slope_linear, data = frontorbcort, id = list(n = 5)) 
frontoc_ptsd_nor <- lm(ptsd ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
frontoc_ptsd_nb <- glm.nb(ptsd ~ CSMP_slope_linear + white + inc_needs, data = frontorbcort)
AIC(frontoc_ptsd_nor, frontoc_ptsd_nb)
summary(frontoc_ptsd_nb) # no

# parahipppost bl -------------------------------------------------------------------------------
# cdi
scatterplot(cdi ~ parahipppost$CSMP_slope_linear, data = parahipppost, id = list(n = 5)) 
parahipp_cdi_nor <- lm(cdi ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
parahipp_cdi_nb <- glm.nb(cdi ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
AIC(parahipp_cdi_nor, parahipp_cdi_nb)
summary(parahipp_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ parahipppost$CSMP_slope_linear, data = parahipppost, id = list(n = 5)) 
parahipp_scared_panic_nor <- lm(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
parahipp_scared_panic_nb <- glm.nb(scared_panic ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
AIC(parahipp_scared_panic_nor, parahipp_scared_panic_nb)
summary(parahipp_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ parahipppost$CSMP_slope_linear, data = parahipppost, id = list(n = 5)) 
parahipp_scared_gad_nor <- lm(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
parahipp_scared_gad_nb <- glm.nb(scared_gad ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
AIC(parahipp_scared_gad_nor, parahipp_scared_gad_nb)
summary(parahipp_scared_gad_nb) # no

# external
scatterplot(external ~ parahipppost$CSMP_slope_linear, data = parahipppost, id = list(n = 5)) 
parahipp_external_nor <- lm(external ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
parahipp_external_nb <- glm.nb(external ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
AIC(parahipp_external_nor, parahipp_external_nb)
exp((1071.626-1072.906)/2) # 53% chance normal is better
summary(parahipp_external_nor) # p = 0.076

# ptsd
scatterplot(ptsd ~ parahipppost$CSMP_slope_linear, data = parahipppost, id = list(n = 5)) 
parahipp_ptsd_nor <- lm(ptsd ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
parahipp_ptsd_nb <- glm.nb(ptsd ~ CSMP_slope_linear + white + inc_needs, data = parahipppost)
AIC(parahipp_ptsd_nor, parahipp_ptsd_nb)
summary(parahipp_ptsd_nb) # no




############################################################.
# REGRESSIONS WITH BLOCK 1 DIFFERENCE; SALIENCE ONLY ####
############################################################.

# amyg rh ---------------------------------------------------------------------------------------

# cdi
scatterplot(cdi ~ amyg_rh$CSPM_1, data = amyg_rh, id = list(n = 5)) 
amyg_cdi_nor <- lm(cdi ~ CSPM_1 + white + inc_needs, data = amyg_rh)
amyg_cdi_nb <- glm.nb(cdi ~ CSPM_1 + white + inc_needs, data = amyg_rh)
AIC(amyg_cdi_nor, amyg_cdi_nb)
summary(amyg_cdi_nb) # nb; no

# scared_panic
scatterplot(scared_panic ~ amyg_rh$CSPM_1, data = amyg_rh, id = list(n = 5)) 
amyg_scared_panic_nor <- lm(scared_panic ~ CSPM_1 + white + inc_needs, data = amyg_rh)
amyg_scared_panic_nb <- glm.nb(scared_panic ~ CSPM_1 + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_panic_nor, amyg_scared_panic_nb)
summary(amyg_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ amyg_rh$CSPM_1, data = amyg_rh, id = list(n = 5)) 
amyg_scared_gad_nor <- lm(scared_gad ~ CSPM_1 + white + inc_needs, data = amyg_rh)
amyg_scared_gad_nb <- glm.nb(scared_gad ~ CSPM_1 + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_gad_nor, amyg_scared_gad_nb)
summary(amyg_scared_gad_nb) # no

# external
scatterplot(external ~ amyg_rh$CSPM_1, data = amyg_rh, id = list(n = 5)) 
amyg_external_nor <- lm(external ~ CSPM_1 + white + inc_needs, data = amyg_rh)
amyg_external_nb <- glm.nb(external ~ CSPM_1 + white + inc_needs, data = amyg_rh)
AIC(amyg_external_nor, amyg_external_nb)
exp((1074.639-1075.777)/2) # 56% chance normal is better
summary(amyg_external_nor) # no

# ptsd
scatterplot(ptsd ~ amyg_rh$CSPM_1, data = amyg_rh, id = list(n = 5)) 
amyg_ptsd_nor <- lm(ptsd ~ CSPM_1 + white + inc_needs, data = amyg_rh)
amyg_ptsd_nb <- glm.nb(ptsd ~ CSPM_1 + white + inc_needs, data = amyg_rh)
AIC(amyg_ptsd_nor, amyg_ptsd_nb)
summary(amyg_ptsd_nb) # no

# insula bl -------------------------------------------------------------------------------------

# cdi
scatterplot(cdi ~ insula$CSPM_1, data = insula, id = list(n = 5)) 
ins_cdi_nor <- lm(cdi ~ CSPM_1 + white + inc_needs, data = insula)
ins_cdi_nb <- glm.nb(cdi ~ CSPM_1 + white + inc_needs, data = insula)
AIC(ins_cdi_nor, ins_cdi_nb)
summary(ins_cdi_nb) # no

# scared_panic
scatterplot(scared_panic ~ insula$CSPM_1, data = insula, id = list(n = 5)) 
ins_scared_panic_nor <- lm(scared_panic ~ CSPM_1 + white + inc_needs, data = insula)
ins_scared_panic_nb <- glm.nb(scared_panic ~ CSPM_1 + white + inc_needs, data = insula)
AIC(ins_scared_panic_nor, ins_scared_panic_nb)
summary(ins_scared_panic_nb) # no

# scared_gad
scatterplot(scared_gad ~ insula$CSPM_1, data = insula, id = list(n = 5)) 
ins_scared_gad_nor <- lm(scared_gad ~ CSPM_1 + white + inc_needs, data = insula)
ins_scared_gad_nb <- glm.nb(scared_gad ~ CSPM_1 + white + inc_needs, data = insula)
AIC(ins_scared_gad_nor, ins_scared_gad_nb)
summary(ins_scared_gad_nb) # no

# external
scatterplot(external ~ insula$CSPM_1, data = insula, id = list(n = 5)) 
ins_external_nor <- lm(external ~ CSPM_1 + white + inc_needs, data = insula)
ins_external_nb <- glm.nb(external ~ CSPM_1 + white + inc_needs, data = insula)
AIC(ins_external_nor, ins_external_nb)
exp((1073.285-1074.495)/2) # 54% chance normal is better
summary(ins_external_nor) # no

# ptsd
scatterplot(ptsd ~ insula$CSPM_1, data = insula, id = list(n = 5)) 
ins_ptsd_nor <- lm(ptsd ~ CSPM_1 + white + inc_needs, data = insula)
ins_ptsd_nb <- glm.nb(ptsd ~ CSPM_1 + white + inc_needs, data = insula)
AIC(ins_ptsd_nor, ins_ptsd_nb)
summary(ins_ptsd_nb) # no


############################################################.
# VISUALIZE OUTCOME VARS: CHANGE IN SXS; T2 CONTROL FOR T1
############################################################.

# my DVS are all metric variables with lower bounds at 0 with no upper bound
# expect negbin versions to fit better for all models, except possibly external 
hist(amyg_rh$cdi_t2, breaks = 20, xlab = "CDI counts", main = "Histogram CDI") # right-skewed
hist(amyg_rh$scared_panic_t2, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC") # right-skewed
hist(amyg_rh$scared_gad_t2, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD") # right-skewed
hist(amyg_rh$external_t2, breaks = 25, xlab = "External counts", main = "Histogram external") # possibly normal
hist(amyg_rh$ptsd_t2, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD sev") # right-skewed


############################################################.
# T2 REGRESSIONS WITH FMRI SLOPE:                      ####
# PREDICT T2, CONTROLLING FOR T1
############################################################.

# amyg rh ---------------------------------------------------------------------------------------

# cdi
amyg_cdi_nor_t2 <- lm(cdi_t2 ~ CSPM_slope_linear + cdi + white + inc_needs, data = amyg_rh)
amyg_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSPM_slope_linear + cdi + white + inc_needs, data = amyg_rh)
AIC(amyg_cdi_nor_t2, amyg_cdi_nb_t2)
summary(amyg_cdi_nb_t2) # no

# scared_panic
amyg_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSPM_slope_linear + scared_panic + white + inc_needs, data = amyg_test)
amyg_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSPM_slope_linear + scared_panic + white + inc_needs, data = amyg_test)
AIC(amyg_scared_panic_nor_t2, amyg_scared_panic_nb_t2)
summary(amyg_scared_panic_nb_t2) # no

# scared_gad
amyg_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSPM_slope_linear + scared_gad + white + inc_needs, data = amyg_rh)
amyg_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSPM_slope_linear + scared_gad + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_gad_nor_t2, amyg_scared_gad_nb_t2)
summary(amyg_scared_gad_nb_t2) # no

# external
amyg_external_nor_t2 <- lm(external_t2 ~ CSPM_slope_linear + external + white + inc_needs, data = amyg_rh)
#amyg_external_nb_t2 <- glm.nb(external_t2 ~ CSPM_slope_linear + external + white + inc_needs, data = amyg_rh) # alternation limit reached
amyg_external_nb_t2 <- gam(external_t2 ~ CSPM_slope_linear + external + white + inc_needs, data = amyg_rh, family = nb)
AIC(amyg_external_nor_t2, amyg_external_nb_t2)
exp((846.8404-848.7251)/2) # 32% chance normal fits best
summary(amyg_external_nb_t2) # 0.002

# ptsd
amyg_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSPM_slope_linear + ptsd + white + inc_needs, data = amyg_rh)
#amyg_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSPM_slope_linear + ptsd + white + inc_needs, data = amyg_rh) # alternation limit reached
amyg_ptsd_nb_t2 <- gam(ptsd_t2 ~ CSPM_slope_linear + ptsd + white + inc_needs, data = amyg_rh, family = nb) # alternation limit reached
AIC(amyg_ptsd_nor_t2, amyg_ptsd_nb_t2)
summary(amyg_ptsd_nb_t2) # no

# insula bl -------------------------------------------------------------------------------------

# cdi
ins_cdi_nor_t2 <- lm(cdi_t2 ~ CSPM_slope_linear + cdi + white + inc_needs, data = insula)
ins_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSPM_slope_linear + cdi + white + inc_needs, data = insula)
AIC(ins_cdi_nor_t2, ins_cdi_nb_t2)
summary(ins_cdi_nb_t2) # no

# scared_panic
ins_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSPM_slope_linear + scared_panic + white + inc_needs, data = insula)
ins_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSPM_slope_linear + scared_panic + white + inc_needs, data = insula)
AIC(ins_scared_panic_nor_t2, ins_scared_panic_nb_t2)
summary(ins_scared_panic_nb_t2) # no

# scared_gad
ins_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSPM_slope_linear + scared_gad + white + inc_needs, data = insula)
ins_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSPM_slope_linear + scared_gad + white + inc_needs, data = insula)
AIC(ins_scared_gad_nor_t2, ins_scared_gad_nb_t2)
summary(ins_scared_gad_nb_t2) # no

# external
ins_external_nor_t2 <- lm(external_t2 ~ CSPM_slope_linear + external + white + inc_needs, data = insula)
ins_external_nb_t2 <- glm.nb(external_t2 ~ CSPM_slope_linear + external + white + inc_needs, data = insula)
AIC(ins_external_nor_t2, ins_external_nb_t2)
exp((851.9331-854.3100)/2) # 30% chance normal fits better
summary(ins_external_nb_t2) # 0.026

# ptsd
ins_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSPM_slope_linear + ptsd + white + inc_needs, data = insula)
ins_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSPM_slope_linear + ptsd + white + inc_needs, data = insula)
AIC(ins_ptsd_nor_t2, ins_ptsd_nb_t2)
summary(ins_ptsd_nb_t2) # no

# hipp rh ---------------------------------------------------------------------------------------
# cdi
hipp_cdi_nor_t2 <- lm(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = hipp_rh)
hipp_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = hipp_rh)
AIC(hipp_cdi_nor_t2, hipp_cdi_nb_t2)
summary(hipp_cdi_nb_t2) # no

# scared_panic
hipp_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = hipp_rh)
hipp_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = hipp_rh)
AIC(hipp_scared_panic_nor_t2, hipp_scared_panic_nb_t2)
summary(hipp_scared_panic_nb_t2) # no

# scared_gad
hipp_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = hipp_rh)
hipp_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = hipp_rh)
AIC(hipp_scared_gad_nor_t2, hipp_scared_gad_nb_t2)
summary(hipp_scared_gad_nb_t2) # no

# external
hipp_external_nor_t2 <- lm(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = hipp_rh)
hipp_external_nb_t2 <- glm.nb(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = hipp_rh)
AIC(hipp_external_nor_t2, hipp_external_nb_t2)
summary(hipp_external_nb_t2) # no

# ptsd
hipp_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = hipp_rh)
hipp_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = hipp_rh)
AIC(hipp_ptsd_nor_t2, hipp_ptsd_nb_t2)
summary(hipp_ptsd_nb_t2) # no

# frontpole bl ----------------------------------------------------------------------------------
# cdi
frontpl_cdi_nor_t2 <- lm(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = frontpole)
frontpl_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = frontpole)
AIC(frontpl_cdi_nor_t2, frontpl_cdi_nb_t2)
summary(frontpl_cdi_nb_t2) # no

# scared_panic
frontpl_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = frontpole)
frontpl_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = frontpole)
AIC(frontpl_scared_panic_nor_t2, frontpl_scared_panic_nb_t2)
summary(frontpl_scared_panic_nb_t2) # no

# scared_gad
frontpl_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = frontpole)
frontpl_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = frontpole)
AIC(frontpl_scared_gad_nor_t2, frontpl_scared_gad_nb_t2)
summary(frontpl_scared_gad_nb_t2) # no

# external
frontpl_external_nor_t2 <- lm(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = frontpole)
frontpl_external_nb_t2 <- glm.nb(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = frontpole)
AIC(frontpl_external_nor_t2, frontpl_external_nb_t2)
exp((856.6284-858.9073)/2) # 31% chance normal fits better
summary(frontpl_external_nb_t2) # no

# ptsd
frontpl_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = frontpole)
frontpl_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = frontpole)
AIC(frontpl_ptsd_nor_t2, frontpl_ptsd_nb_t2)
summary(frontpl_ptsd_nb_t2) # no

# frontorbcort bl -------------------------------------------------------------------------------
# cdi
frontoc_cdi_nor_t2 <- lm(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = frontorbcort)
frontoc_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = frontorbcort)
AIC(frontoc_cdi_nor_t2, frontoc_cdi_nb_t2)
summary(frontoc_cdi_nb_t2) # no

# scared_panic
frontoc_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = frontorbcort)
frontoc_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = frontorbcort)
AIC(frontoc_scared_panic_nor_t2, frontoc_scared_panic_nb_t2)
summary(frontoc_scared_panic_nb_t2) # no

# scared_gad
frontoc_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = frontorbcort)
frontoc_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = frontorbcort)
AIC(frontoc_scared_gad_nor_t2, frontoc_scared_gad_nb_t2)
summary(frontoc_scared_gad_nb_t2) # no

# external
frontoc_external_nor_t2 <- lm(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = frontorbcort)
frontoc_external_nb_t2 <- glm.nb(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = frontorbcort)
AIC(frontoc_external_nor_t2, frontoc_external_nb_t2)
summary(frontoc_external_nb_t2) # no

# ptsd
frontoc_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = frontorbcort)
frontoc_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = frontorbcort)
AIC(frontoc_ptsd_nor_t2, frontoc_ptsd_nb_t2)
summary(frontoc_ptsd_nb_t2) # no

# parahipppost bl -------------------------------------------------------------------------------
# cdi
parahipp_cdi_nor_t2 <- lm(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = parahipppost)
parahipp_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSMP_slope_linear + cdi + white + inc_needs, data = parahipppost)
AIC(parahipp_cdi_nor_t2, parahipp_cdi_nb_t2)
summary(parahipp_cdi_nb_t2) # 0.019

# scared_panic
parahipp_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = parahipppost)
parahipp_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSMP_slope_linear + scared_panic + white + inc_needs, data = parahipppost)
AIC(parahipp_scared_panic_nor_t2, parahipp_scared_panic_nb_t2)
summary(parahipp_scared_panic_nb_t2) # no

# scared_gad
parahipp_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = parahipppost)
parahipp_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSMP_slope_linear + scared_gad + white + inc_needs, data = parahipppost)
AIC(parahipp_scared_gad_nor_t2, parahipp_scared_gad_nb_t2)
summary(parahipp_scared_gad_nb_t2) # no

# external
parahipp_external_nor_t2 <- lm(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = parahipppost)
parahipp_external_nb_t2 <- glm.nb(external_t2 ~ CSMP_slope_linear + external + white + inc_needs, data = parahipppost)
AIC(parahipp_external_nor_t2, parahipp_external_nb_t2)
summary(parahipp_external_nb_t2) # no

# ptsd
parahipp_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = parahipppost)
parahipp_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSMP_slope_linear + ptsd + white + inc_needs, data = parahipppost)
AIC(parahipp_ptsd_nor_t2, parahipp_ptsd_nb_t2)
summary(parahipp_ptsd_nb_t2) # no





############################################################.
# T2 REGRESSIONS WITH BLOCK 1 DIFFERENCE:               ####
# PREDICT T2, CONTROLLING FOR T1; SALIENCE ONLY
############################################################.

# amyg rh ---------------------------------------------------------------------------------------

# cdi
amyg_cdi_nor_t2 <- lm(cdi_t2 ~ CSPM_1 + cdi + white + inc_needs, data = amyg_rh)
amyg_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSPM_1 + cdi + white + inc_needs, data = amyg_rh)
AIC(amyg_cdi_nor_t2, amyg_cdi_nb_t2)
summary(amyg_cdi_nor_t2) # no

# scared_panic
amyg_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSPM_1 + scared_panic + white + inc_needs, data = amyg_rh)
amyg_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSPM_1 + scared_panic + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_panic_nor_t2, amyg_scared_panic_nb_t2)
summary(amyg_scared_panic_nb_t2) # no

# scared_gad
amyg_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSPM_1 + scared_gad + white + inc_needs, data = amyg_rh)
amyg_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSPM_1 + scared_gad + white + inc_needs, data = amyg_rh)
AIC(amyg_scared_gad_nor_t2, amyg_scared_gad_nb_t2)
summary(amyg_scared_gad_nb_t2) # no

# external
amyg_external_nor_t2 <- lm(external_t2 ~ CSPM_1 + external + white + inc_needs, data = amyg_rh)
amyg_external_nb_t2 <- glm.nb(external_t2 ~ CSPM_1 + external + white + inc_needs, data = amyg_rh)
AIC(amyg_external_nor_t2, amyg_external_nb_t2)
summary(amyg_external_nb_t2) # no

# ptsd
amyg_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSPM_1 + ptsd + white + inc_needs, data = amyg_rh)
amyg_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSPM_1 + ptsd + white + inc_needs, data = amyg_rh)
AIC(amyg_ptsd_nor_t2, amyg_ptsd_nb_t2)
summary(amyg_ptsd_nb_t2) # no

# insula bl -------------------------------------------------------------------------------------

# cdi
ins_cdi_nor_t2 <- lm(cdi_t2 ~ CSPM_1 + cdi + white + inc_needs, data = insula)
ins_cdi_nb_t2 <- glm.nb(cdi_t2 ~ CSPM_1 + cdi + white + inc_needs, data = insula)
AIC(ins_cdi_nor_t2, ins_cdi_nb_t2)
summary(ins_cdi_nb_t2) # no

# scared_panic
ins_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ CSPM_1 + scared_panic + white + inc_needs, data = insula)
ins_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ CSPM_1 + scared_panic + white + inc_needs, data = insula)
AIC(ins_scared_panic_nor_t2, ins_scared_panic_nb_t2)
summary(ins_scared_panic_nb_t2) # no

# scared_gad
ins_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ CSPM_1 + scared_gad + white + inc_needs, data = insula)
ins_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ CSPM_1 + scared_gad + white + inc_needs, data = insula)
AIC(ins_scared_gad_nor_t2, ins_scared_gad_nb_t2)
summary(ins_scared_gad_nb_t2) # no

# external
ins_external_nor_t2 <- lm(external_t2 ~ CSPM_1 + external + white + inc_needs, data = insula)
ins_external_nb_t2 <- glm.nb(external_t2 ~ CSPM_1 + external + white + inc_needs, data = insula)
AIC(ins_external_nor_t2, ins_external_nb_t2)
summary(ins_external_nb_t2) # no 

# ptsd
ins_ptsd_nor_t2 <- lm(ptsd_t2 ~ CSPM_1 + ptsd + white + inc_needs, data = insula)
ins_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ CSPM_1 + ptsd + white + inc_needs, data = insula)
AIC(ins_ptsd_nor_t2, ins_ptsd_nb_t2)
summary(ins_ptsd_nb_t2) # 0.036


############################################################.
# ADDED VARIABLE PLOTS OF SIGNIFICANT RESULTS           ####
############################################################.
# Added variable plots (aka partial regression plots, adjusted variable plots, individual coefficient plots), 
# are “results” plots. They are plots showing the estimated relationship between the response and an explanatory 
# variable after accounting for the other variables in the model. Use these to plot the findings, 
# holding covariates white and inc_needs constant (plotting for the median income, white participant).
# https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
# If you plot like this it doesn't account for the influence of race or inc_needs on the model:
# ggplot(brain_data, aes(x = brain_data, y = symptom)) + geom_point() + theme_classic() +
# labs(x = "Brain data", y = "Symptom") + geom_smooth(method = MASS::glm.nb, se = TRUE)


# RE-RUN CORRECT MODELS FIRST; NAMES OVERWRITE IN EACH ANALYSIS SECTION

# PLOT SIGNIFICANT MODELS ONLY

#--------------------------------------------------------------------------------------------.
# Blunted habituation/safety slopes predicting baseline symptoms (NONE)

#--------------------------------------------------------------------------------------------.
# Baseline block1 discrimination predicting baseline symptoms (NONE)

#--------------------------------------------------------------------------------------------.
# Baseline block1 discrimination predicting change in symptoms, controlling for T1 symptoms:

# bl insula and ptsd 
newdata <- data.frame(CSPM_1 = insula$CSPM_1, white = factor(1), inc_needs = median(insula$inc_needs), ptsd = mean(insula$ptsd))
newdata$p.hat <- predict(ins_ptsd_nb_t2, newdata, type = "response")
newdata <- cbind(newdata, predict(ins_ptsd_nb_t2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  ptsd_t2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(CSPM_1, ptsd_t2)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = insula, colour = 'black') +
  labs(x = "Bilateral insula block 1 CS+ vs. CS-", y = "PTSD at 2-year follow-up")

#ggsave("plots/symptom-regressions/fMRI-calc-slope/discrim_insula-bl_T2_ptsd.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

#--------------------------------------------------------------------------------------------.
# Baseline habituation/safety slopes predicting change in symptoms, controlling for T1 symptoms:

# rh amyg and external
newdata <- data.frame(CSPM_slope_linear = amyg_rh$CSPM_slope_linear, white = factor(1), inc_needs = median(amyg_rh$inc_needs), external = mean(amyg_rh$external))
newdata$p.hat <- predict(amyg_external_nb_t2, newdata, type = "response")
newdata <- cbind(newdata, predict(amyg_external_nb_t2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  external_t2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(CSPM_slope_linear, external_t2)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = amyg_rh, colour = 'black') +
  labs(x = "Right amygdala habituation slope", y = "Externalizing at 2-year follow-up")

#ggsave("plots/symptom-regressions/fMRI-calc-slope/slope_amyg-rh_T2_external.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# insula and external
newdata <- data.frame(CSPM_slope_linear = insula$CSPM_slope_linear, white = factor(1), inc_needs = median(insula$inc_needs), external = mean(insula$external))
newdata$p.hat <- predict(ins_external_nb_t2, newdata, type = "response")
newdata <- cbind(newdata, predict(ins_external_nb_t2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  external_t2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(CSPM_slope_linear, external_t2)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = insula, colour = 'black') +
  labs(x = "Bilateral insula habituation slope", y = "Externalizing at 2-year follow-up")

#ggsave("plots/symptom-regressions/fMRI-calc-slope/slope_insula-bl_T2_external.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# parahipppost and cdi
newdata <- data.frame(CSMP_slope_linear = parahipppost$CSMP_slope_linear, white = factor(1), inc_needs = median(parahipppost$inc_needs), cdi = mean(parahipppost$cdi))
newdata$p.hat <- predict(parahipp_cdi_nb_t2, newdata, type = "response")
newdata <- cbind(newdata, predict(parahipp_cdi_nb_t2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  cdi_t2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(CSMP_slope_linear, cdi_t2)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = parahipppost, colour = 'black') +
  labs(x = "Bilateral posterior parahippocampal gyrus cortex safety-signaling slope", y = "Depression at 2-year follow-up")

#ggsave("plots/symptom-regressions/fMRI-calc-slope/slope_parahipppost-bl_T2_cdi.pdf", width = 4.5, height = 4, units = "in", dpi = 600)



############################################################.
# SUMMARY (FINDINGS ONLY; REMOVED OTHER ROWS)
############################################################.
#---------------------------------------------------------------------------------------------------------.
# SLOPE LINEAR              cdi            scared_panic      scared_gad      external        ptsd
# amyg rh                                                                                    0.080         << TREND (up slope, up sxs)                                                                                     
# parahipppost bl                                                            0.076                         << TREND (up slope, down sxs)
#---------------------------------------------------------------------------------------------------------.
# BLOCK1 DIFF               cdi            scared_panic      scared_gad      external        ptsd
#---------------------------------------------------------------------------------------------------------.
# CHANGE IN SXS; SLOPE      cdi            scared_panic      scared_gad      external        ptsd
# amyg rh                                                                    0.002                         (up slope, up sxs)
# insula bl                                                                  0.026                         (up slope, up sxs)
# parahipppost bl           0.019                                                                          (up slope, up sxs)?
#---------------------------------------------------------------------------------------------------------.
# CHANGE IN SXS; BLOCK1 DIFFcdi            scared_panic      scared_gad      external        ptsd
# insula bl                                                                                  0.036         (up discrim, down sxs)
#---------------------------------------------------------------------------------------------------------.



############################################################.
# SUMMARY (TRENDS REMOVED)
############################################################.
#---------------------------------------------------------------------------------------------------------.
# SLOPE LINEAR (NONE)       cdi            scared_panic      scared_gad      external        ptsd
#---------------------------------------------------------------------------------------------------------.
# BLOCK1 DIFF (NONE)        cdi            scared_panic      scared_gad      external        ptsd
#---------------------------------------------------------------------------------------------------------.
# CHANGE IN SXS; SLOPE      cdi            scared_panic      scared_gad      external        ptsd
# amyg rh                                                                    0.002                         (up slope, up sxs)
# insula bl                                                                  0.026                         (up slope, up sxs)
# parahipppost bl           0.019                                                                          (up slope, up sxs)?
#---------------------------------------------------------------------------------------------------------.
# CHANGE IN SXS; BLOCK1 DIFFcdi            scared_panic      scared_gad      external        ptsd
# insula bl                                                                                  0.036         (up discrim, down sxs)
#---------------------------------------------------------------------------------------------------------.










