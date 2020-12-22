##############################################################
# INFO
##############################################################

# RE: JEANETTE MUMFORD 
# gPPI preferable; but cannot average over the whole block, must be each TR
# Truncate the TRs to account for hemodynamic response; shift to nearest TR

# In order to make different plots exploring the different findings from the whole-brain maps,
# edit the patterns searched for in loading the files for roi1 and roi2, and edit the last
# few lines of code in this script to change the titles of the x and y axes in the plots.
# Load the seed region as roi1, which will end up on the x axis.

# Modeled after http://www.brainmapping.org/NITP/2011Data/ppi_fsl_lab_2011.pdf
# https://fmri-training-course.psych.lsa.umich.edu/wp-content/uploads/2018/08/mumford_thursday_ppi.pdf
# Question of interest: do the slopes between red line and blue line differ?
# Whole brain maps are picking up on whether: difference in slopes for CTL is DIFFERENT than difference in slopes for TR.


##############################################################
# PREP & LOAD
##############################################################
rm(list=ls())
if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
# purrr is the part of tidyverse with pmap (lets function iterate over multiple inputs simultaneously)
source("code/tools/gPPI-plot-functions.R")  # functions live here


filenames_roi1_ctl <- list.files(path = "data/roi-extraction-files/gPPI/visualization/control/", pattern="*_rh_amyg.txt", full.names = TRUE) # make seed region roi1
filenames_roi2_ctl <- list.files(path = "data/roi-extraction-files/gPPI/visualization/control/", pattern="_acc.txt", full.names = TRUE)
filenames_task_ctl <- list.files(path = "data/roi-extraction-files/gPPI/visualization/control/", pattern="*_task.txt", full.names = TRUE)

filenames_roi1_tr <- list.files(path = "data/roi-extraction-files/gPPI/visualization/trauma/", pattern="*_rh_amyg.txt", full.names = TRUE) # make seed region roi1
filenames_roi2_tr <- list.files(path = "data/roi-extraction-files/gPPI/visualization/trauma/", pattern="*_acc.txt", full.names = TRUE)
filenames_task_tr <- list.files(path = "data/roi-extraction-files/gPPI/visualization/trauma/", pattern="*_task.txt", full.names = TRUE)

roi1_ctl <- lapply(filenames_roi1_ctl, read.table) 
roi2_ctl <- lapply(filenames_roi2_ctl, read.table) 
task_ctl <- lapply(filenames_task_ctl, read.table) 

roi1_tr <- lapply(filenames_roi1_tr, read.table) 
roi2_tr <- lapply(filenames_roi2_tr, read.table) 
task_tr <- lapply(filenames_task_tr, read.table) 


##############################################################
# DATA WRANGLE
##############################################################

# Manipulate roi1 and roi2 ----------------------------------------------------------------

# Names cols, adds TR col, demeans timeseries
roi1_ctl <- lapply(roi1_ctl, SetUpROI)
roi2_ctl <- lapply(roi2_ctl, SetUpROI)
roi1_tr <- lapply(roi1_tr, SetUpROI)
roi2_tr <- lapply(roi2_tr, SetUpROI)

# Manipulate task --------------------------------------------------------------------------

# Names cols, adds cols for stim type and block number, sorts rows
task_ctl <- lapply(task_ctl, SetUpTask) 
task_tr <- lapply(task_tr, SetUpTask) 

# Adds col listing the TR that begins each respective block
task_ctl <- lapply(task_ctl, StartBlocks)
task_tr <- lapply(task_tr, StartBlocks)

# Each block has 10 TRs (2 sec each) - pull TRs 4-10 (leave out first 6 sec to account for HRF)
# For each subject individually, get list of roi1 CS+ values and CS- values, and roi2 CS+ values and CS- values
data_ctl <- pmap(list(x = task_ctl, y = roi1_ctl, z = roi2_ctl), ActivityDuringBlocks)
data_tr <- pmap(list(x = task_tr, y = roi1_tr, z = roi2_tr), ActivityDuringBlocks)


##############################################################
# AVERAGING OVER CONTROLS/TRAUMA
##############################################################

# Make skeleton of final dataframes
ctl_avg <- data_ctl[[1]]
ctl_avg[ , "roi1"] <- NA
ctl_avg[ , "roi2"] <- NA
tr_avg <- data_tr[[1]]
tr_avg[ , "roi1"] <- NA
tr_avg[ , "roi2"] <- NA



# CONTROL, ROI1 ---------------------------------------
i <- 1 # counter for the 56 rows of data
j <- 1 # counter for the 70 control subjects
vec <- rep(NA, 70) # make vector of 70 empty NAs for that row for each subject

for (i in 1:56){ # repeat over 56 rows in order 
  
  for (j in 1:70){ # take the value from each subject and stick in temp vector "vec"
    vec[j] <- data_ctl[[j]][i, "roi1"]    
  }

  ctl_avg[i, "roi1"] <- mean(vec) # take mean of temp vector and save it in final df
  vec <- rep(NA, 70) # re-make vector of 70 empty NAs
  
  i <- i + 1 # move on to next row out of 56 rows (turns out this line is redundant; do not need)
  
}

# CONTROL, ROI2 ---------------------------------------
i <- 1 # counter for the 56 rows of data
j <- 1 # counter for the 70 control subjects
vec <- rep(NA, 70) # make vector of 70 empty NAs for that row for each subject

for (i in 1:56){ # repeat over 56 rows in order
  
  for (j in 1:70){ # take the value from each subject and stick in temp vector "vec"
    vec[j] <- data_ctl[[j]][i, "roi2"]    
  }
  
  ctl_avg[i, "roi2"] <- mean(vec) # take mean of temp vector and save it in final df
  vec <- rep(NA, 70) # re-make vector of 70 empty NAs
  
  i <- i + 1 # move on to next row out of 56 rows (turns out this line is redundant; do not need)
  
}


# TRAUMA, ROI1 ---------------------------------------
i <- 1 # counter for the 56 rows of data
j <- 1 # counter for the 77 trauma subjects
vec <- rep(NA, 77) # make vector of 77 empty NAs for that row for each subject

for (i in 1:56){ # repeat over 56 rows in order
  
  for (j in 1:77){ # take the value from each subject and stick in temp vector "vec"
    vec[j] <- data_tr[[j]][i, "roi1"]    
  }
  
  tr_avg[i, "roi1"] <- mean(vec) # take mean of temp vector and save it in final df
  vec <- rep(NA, 77) # re-make vector of 77 empty NAs
  
  i <- i + 1 # move on to next row out of 56 rows (turns out this line is redundant; do not need)
  
}


# TRAUMA, ROI2 ---------------------------------------
i <- 1 # counter for the 56 rows of data
j <- 1 # counter for the 77 trauma subjects
vec <- rep(NA, 77) # make vector of 77 empty NAs for that row for each subject

for (i in 1:56){ # repeat over 56 rows in order
  
  for (j in 1:77){ # take the value from each subject and stick in temp vector "vec"
    vec[j] <- data_tr[[j]][i, "roi2"]    
  }
  
  tr_avg[i, "roi2"] <- mean(vec) # take mean of temp vector and save it in final df
  vec <- rep(NA, 77) # re-make vector of 77 empty NAs
  
  i <- i + 1 # move on to next row out of 56 rows (turns out this line is redundant; do not need)
  
}


##############################################################
# PLOTS
##############################################################

# Combining into one dataframe to put both groups in one plot
ctl_avg <- ctl_avg %>%
  mutate(trauma = 0)
tr_avg <- tr_avg %>%
  mutate(trauma = 1)
data_combined <- bind_rows(ctl_avg, tr_avg)
trauma_names <- c(`0` = "Control", `1` = "Trauma")

# Creating base plot
plotPPI_base <- ggplot(data_combined, aes(x = roi1, y = roi2, color = stim)) +
  theme_classic(base_size = 12) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", fill = NA, fullrange = TRUE) + # fullrange makes the line extend across the whole plot
  theme(panel.border = element_rect(fill=NA)) + # completes the border into a box shape
  theme(plot.title = element_text(hjust = 0.5)) + # centers title
  scale_color_manual(name = "Stimuli", values = c("CS-" = "blue", "CS+" = "red")) +
  guides(color = guide_legend(override.aes = list(size = 2, shape = 15, linetype = 0), reverse = TRUE)) + # legend squares with no line, red on top 
  theme(legend.title = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  facet_grid(~trauma, labeller = as_labeller(trauma_names))

# Edit this last section here to customize for each graph
plotPPI <- plotPPI_base +
  xlab("Right amygdala BOLD response") +                  # roi1 (seed region)
  ylab("Bilateral ACC BOLD response") +           # roi2
  ggtitle("Generalized psychophysiological interaction") #+
  #coord_cartesian(ylim=c(-25,25))                           # controls the zoom

# Creation/previewing and saving (filename is "seed-finding.pdf")
plotPPI
ggsave(filename = "FINAL_rh-amyg_bl-acc.pdf", plot = plotPPI, width = 6.5, height = 4.5, units = "in", dpi = 600)




