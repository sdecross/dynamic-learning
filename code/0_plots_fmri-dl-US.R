##############################
# SOURCE & LOAD DATA
##############################

rm(list=ls())

source("code/tools/fMRI-plot-functions-US.R")

if (!require("pacman")) {install.packages("pacman"); require("pacman")}
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm",
       "coin", "gridExtra", "tidyverse")

subjectgroups <- read.table("data/behavioral/FearLearning_SubjectGroups_n147.txt", header=TRUE)
subjectgroups[ , "trauma"] <- factor(subjectgroups[ , "trauma"]) # convert to factor

orig_amyg_lh <- read.table("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DL-US_Final_Amygdala_lh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_FearLearning_US.feat.txt", header = T)
orig_amyg_rh <- read.table("data/roi-extraction-files/fMRI/dl-func-rois/FearLearning_DL-US_Final_Amygdala_rh_H-O_thr20_b_x_FearRevPM_Alln147_p05_CSMvCSP_FearLearning_US.feat.txt", header = T)

# naming impt variables, making named list to refer to data later (when multiple ROIs)
roilist <- list(orig_amyg_lh, orig_amyg_rh)
names(roilist) <- c("amyg_lh", "amyg_rh")
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
                  CSPR_1 = CSplusR_1_zstat1, CSPR_2 = CSplusR_2_zstat2, CSPR_3 = CSplusR_3_zstat3, CSPR_4 = CSplusR_4_zstat4,
                  CSPN_1 = CSplusNR_1_zstat5, CSPN_2 = CSplusNR_2_zstat6, CSPN_3 = CSplusNR_3_zstat7, CSPN_4 = CSplusNR_4_zstat8) %>% 
    mutate(CSPRvPN_1 = CSPR_1 - CSPN_1) %>% # adding CSPR vs. CSP differential cols
    mutate(CSPRvPN_2 = CSPR_2 - CSPN_2) %>%
    mutate(CSPRvPN_3 = CSPR_3 - CSPN_3) %>%
    mutate(CSPRvPN_4 = CSPR_4 - CSPN_4) 
}

##############################
# PLOTS
##############################

# need titles for graphs that go in same order as list elements
titles <- c("Left amygdala", "Right amygdala")

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
       function(x) ggsave(filename=paste("plots/fMRI/US_plots/diff_plots/", x, "_diff.pdf", sep=""), plot=plots_diff[[x]], width=6.5, height=4.5, units="in", dpi=600))

plots_stim <- lapply(stim_final, stim_plot)
plots_stim 
lapply(names(plots_stim), 
       function(x) ggsave(filename=paste("plots/fMRI/US_plots/stim_plots/", x, "_stim.pdf", sep=""), plot=plots_stim[[x]], width=6.5, height=4.5, units="in", dpi=600))



#############################################
# REFINE AND SAVE DATA FOR LATER
#############################################
# name trauma factor levels
for (roi in seq_along(roilist)){
  levels(roilist[[roi]]$trauma) <- c("Control", "Trauma")
}

# rename object to save
fMRI_dl_func_roilist_US <- roilist

# save the whole roilist list 
#save(fMRI_dl_func_roilist_US, file = "data/full-data/fMRI_dl_func_roilist_US.RData") 



