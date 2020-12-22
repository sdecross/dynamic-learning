############################################################
# SOURCE & LOAD & MANIPULATE DATA
############################################################
rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", 
       "lmerTest", "effects", "MASS", "mgcv", "splines", "sjPlot", "tidyverse")

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

insula <- roilist_regressions$insula_bl


# Does trauma predict block 1 discrimination?

hist(insula$CSPM_1, breaks = 20, xlab = "BOLD", main = "Histogram block 1 discrim") # normal

fit <- lm(CSPM_1 ~ trauma + white + inc_needs, data = insula)
summary(fit) # no

# Doing this to determine if a mediation model should be tested for activation metric of block 1 discrimination in salience regions.
# Already ran arm b tests at this point to see whether block 1 discrimination in salience regions predicted
# psychopathology; insula was only region where relevant findings for discrim and therefore relevant to test a mediation model.
# So that's why insula is only salience region here; this analysis was run second.




