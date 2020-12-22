##############################.
# SOURCE & LOAD DATA
##############################.

rm(list=ls())

source("code/tools/fMRI-plot-functions.R")

if (!require("pacman")) {install.packages("pacman"); require("pacman")}
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm",
       "coin", "gridExtra", "tidyverse")

subjectgroups <- read.table("data/behavioral/FearLearning_SubjectGroups_n147.txt", header=TRUE)
subjectgroups[ , "trauma"] <- factor(subjectgroups[ , "trauma"]) # convert to factor

orig_ACC_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_ACC_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_amyg_lh <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Amygdala_lh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_amyg_rh <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Amygdala_rh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_frontorbcort_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_FrontOrbCort_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_frontpole_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_FrontPole_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_hipp_lh <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Hippocampus_lh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_hipp_rh <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Hippocampus_rh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_insula_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Insula_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_parahipppost_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Parahipp_post_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_precunpcc_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_PrecuneusPCC_VentralBlob_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_subcallosal_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Subcallosal_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")
orig_thalamus_bl <- read.csv("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DynamicLearning_Thalamus_bl_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_in_FearLearning_space.txt")

# naming impt variables, making named list to refer to data later (when multiple ROIs)
roilist <- list(orig_ACC_bl, orig_amyg_lh, orig_amyg_rh, orig_frontorbcort_bl, orig_frontpole_bl, orig_hipp_lh, orig_hipp_rh, orig_insula_bl,
                orig_parahipppost_bl, orig_precunpcc_bl, orig_subcallosal_bl, orig_thalamus_bl)
names(roilist) <- c("ACC_bl", "amyg_lh", "amyg_rh", "frontorbcort_bl", "frontpole_bl", "hipp_lh", "hipp_rh", "insula_bl", "parahipppost_bl", 
                    "precunpcc_bl", "subcallosal_bl", "thalamus_bl")
trauma_names <- c(`0` = "Control", `1` = "Trauma")
block_names <- c("1", "2", "3", "4")

# re-name vars to merge with roi lists
subjectgroups <- subjectgroups %>%
  dplyr::select(Subject = subject, Trauma = trauma)

# manipulating loaded data to obtain variables of interest
for (roi in seq_along(roilist)){
  roilist[[roi]] <- roilist[[roi]] %>%
    left_join(subjectgroups, by = "Subject") %>% # add in trauma/grouping col, matching to subject ID
    mutate(Subject = as.factor(Subject)) %>% # change subject IDs to factor
    dplyr::select(subject = Subject, trauma = Trauma, 
                  CSP_1 = CSplusNR_1_zstat5, CSP_2 = CSplusNR_2_zstat6, CSP_3 = CSplusNR_3_zstat7, CSP_4 = CSplusNR_4_zstat8, 
                  CSM_1 = CSminus_1_zstat1, CSM_2 = CSminus_2_zstat2, CSM_3 = CSminus_3_zstat3, CSM_4 = CSminus_4_zstat4) %>%
    mutate(CSPM_1 = CSP_1 - CSM_1) %>% # adding CSPM differential cols
    mutate(CSPM_2 = CSP_2 - CSM_2) %>%
    mutate(CSPM_3 = CSP_3 - CSM_3) %>%
    mutate(CSPM_4 = CSP_4 - CSM_4) %>%
    mutate(CSP_WR = (CSP_1 + CSP_2 + CSP_3 + CSP_4)/4) %>% # adding cols for the avg of CS+ and CS- over the whole run
    mutate(CSM_WR = (CSM_1 + CSM_2 + CSM_3 + CSM_4)/4) %>%
    mutate(CSPM_WR = CSP_WR - CSM_WR) # adding CSPM differential for avg over the whole run
}



##############################.
# PLOTS                   ####
##############################.

# need titles for graphs that go in same order as list elements
titles <- c("Bilateral ACC", "Left amygdala", "Right amygdala", "Bilateral frontal orbital cortex", "Bilateral frontal pole", "Left hippocampus", 
            "Right hippocampus", "Bilateral insula", "Bilateral posterior parahippocampal gyrus", "Bilateral Precuneus/PCC", "Bilateral subcallosal cortex",
            "Bilateral thalamus")

# plot prep
diff_ctl_list <- lapply(roilist, diff_ctl)
diff_tr_list <- lapply(roilist, diff_tr)
diff_ctl_summary <- lapply(diff_ctl_list, diff_ctl_sum)
diff_tr_summary <- lapply(diff_tr_list, diff_tr_sum)
diff_summary <- map2(diff_ctl_summary, diff_tr_summary, combine_summary_diff) # map2 simultaneously iterates a function across two parameters in parallel
diff_final <- Map(list, diff_summary, titles) # glue the summary dataframe list and the titles into a bigger master list

stim_ctl_list <- lapply(roilist, stim_ctl)
stim_tr_list <- lapply(roilist, stim_tr)
stim_ctl_summary <- lapply(stim_ctl_list, stim_ctl_sum)
stim_tr_summary <- lapply(stim_tr_list, stim_tr_sum)
stim_pre1_summary <- map2(stim_ctl_summary, stim_tr_summary, combine_summary_stim)
stim_pre2_summary <- lapply(stim_pre1_summary, addstimsign)
stim_summary <- lapply(stim_pre2_summary, addblockcol)
stim_final <- Map(list, stim_summary,titles)

# plot creation and saving
plots_diff <- lapply(diff_final, diff_plot)
plots_diff 
lapply(names(plots_diff),
       function(x) ggsave(filename=paste("plots/fMRI/diff_plots/", x, "_diff.pdf", sep=""), plot=plots_diff[[x]], width=6.5, height=4.5, units="in", dpi=600))

plots_stim <- lapply(stim_final, stim_plot)
plots_stim 
lapply(names(plots_stim), 
       function(x) ggsave(filename=paste("plots/fMRI/stim_plots/", x, "_stim.pdf", sep=""), plot=plots_stim[[x]], width=6.5, height=4.5, units="in", dpi=600))

# -------------------------------------------------------------------------.
# create version of plots emphasizing safety-signaling; upwards trajectory

# pull orig files for regions to make safety plots of into new safety list
roilist_safety <- list(orig_frontorbcort_bl, orig_frontpole_bl, orig_hipp_lh, orig_hipp_rh, orig_parahipppost_bl, 
                       orig_precunpcc_bl, orig_subcallosal_bl) 
names(roilist_safety) <- c("frontorbcort_bl", "frontpole_bl", "hipp_lh", "hipp_rh", "parahipppost_bl", "precunpcc_bl", "subcallosal_bl")

# manipulating loaded data to obtain variables of interest; MP instead of PM differentials needed
for (roi in seq_along(roilist_safety)){
  roilist_safety[[roi]] <- roilist_safety[[roi]] %>%
    left_join(subjectgroups, by = "Subject") %>% # add in trauma/grouping col, matching to subject ID
    mutate(Subject = as.factor(Subject)) %>% # change subject IDs to factor
    dplyr::select(subject = Subject, trauma = Trauma, 
                  CSP_1 = CSplusNR_1_zstat5, CSP_2 = CSplusNR_2_zstat6, CSP_3 = CSplusNR_3_zstat7, CSP_4 = CSplusNR_4_zstat8, 
                  CSM_1 = CSminus_1_zstat1, CSM_2 = CSminus_2_zstat2, CSM_3 = CSminus_3_zstat3, CSM_4 = CSminus_4_zstat4) %>%
    mutate(CSMP_1 = CSM_1 - CSP_1) %>% # adding CSMP differential cols
    mutate(CSMP_2 = CSM_2 - CSP_2) %>%
    mutate(CSMP_3 = CSM_3 - CSP_3) %>%
    mutate(CSMP_4 = CSM_4 - CSP_4) 
}
names(roilist_safety) <- c("frontorbcort_bl", "frontpole_bl", "hipp_lh", "hipp_rh", "parahipppost_bl", "precunpcc_bl", "subcallosal_bl")

# need titles for graphs that go in same order as list elements
titles_safety <- c("Bilateral frontal orbital cortex", "Bilateral frontal pole", "Left hippocampus", "Right hippocampus", 
                   "Bilateral posterior parahippocampal gyrus", "Bilateral PCC", "Bilateral subcallosal cortex")

# plot prep for safety-signaling
diff_ctl_list_safety <- lapply(roilist_safety, diff_ctl_safety)
diff_tr_list_safety <- lapply(roilist_safety, diff_tr_safety)
diff_ctl_summary_safety <- lapply(diff_ctl_list_safety, diff_ctl_sum)
diff_tr_summary_safety <- lapply(diff_tr_list_safety, diff_tr_sum)
diff_summary_safety <- map2(diff_ctl_summary_safety, diff_tr_summary_safety, combine_summary_diff) # map2 simultaneously iterates a function across two parameters in parallel
diff_final_safety <- Map(list, diff_summary_safety, titles_safety) # glue the summary dataframe list and the titles into a bigger master list

# plot creation and safety
plots_diff_safety <- lapply(diff_final_safety, diff_plot_safety)
plots_diff_safety 
lapply(names(plots_diff_safety),
       function(x) ggsave(filename=paste("plots/fMRI/diff_plots/", x, "_diff_safety.pdf", sep=""), plot=plots_diff_safety[[x]], width=6.5, height=4.5, units="in", dpi=600))


#############################################.
# REFINE AND SAVE DATA FOR LATER
#############################################.

# name trauma factor levels
for (roi in seq_along(roilist)){
  levels(roilist[[roi]]$trauma) <- c("Control", "Trauma")
}

# rename object to save
fMRI_dl_func_roilist <- roilist
fMRI_dl_func_roilist_safety <- roilist_safety

# save the whole roilist list 
#save(fMRI_dl_func_roilist, file = "data/full-data/fMRI_dl_func_roilist.RData") # differentials CSPM
#save(fMRI_dl_func_roilist_safety, file = "data/full-data/fMRI_dl_func_roilist_safety.RData") # opposite direction differentials (CSMP)

