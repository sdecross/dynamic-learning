###################################################################.
# fMRI MEDIATION ----------------- --------------------------------
###################################################################.

############################################################.
# SOURCE & LOAD & MANIPULATE DATA
############################################################.
rm(list=ls())

require(pacman)
p_load("pander", "ggplot2", "ggthemes", "data.table", "purrr", "lmPerm", "car", "coin", "gridExtra", 
       "lmerTest", "effects", "MASS", "mgcv", "splines", "sjPlot", "gvlma", "mediation", "tidyverse")

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "data"

fear_data <- data %>%
  dplyr::select(-trauma) # get rid of "trauma" col; duplicate col causes trouble in merge of dfs with fMRI data

fear_data_2 <- data # doesn't get rid of "trauma" col; PPI data needs this included


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

############################################################.
# fMRI baseline mediation models                        ####
############################################################.

# for models where a and b branches are significant or trending

# trauma --- rh amyg slope --- ptsd
set.seed(999)
xm <- lm(CSPM_slope_linear ~ trauma + white + inc_needs, data = amyg_rh)  
xmy <- lm(ptsd ~ CSPM_slope_linear + trauma + white + inc_needs, data = amyg_rh) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "CSPM_slope_linear", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no; ACME 0.92
plot(m1)


# trauma --- vvs slope --- external
set.seed(999)
xm <- lm(CSMP_slope_linear ~ trauma + white + inc_needs, data = parahipppost)  
xmy <- lm(external ~ CSMP_slope_linear + trauma + white + inc_needs, data = parahipppost) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "CSMP_slope_linear", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no; ACME 0.44
plot(m1)

############################################################.
# fMRI T2 mediation models                              ####
############################################################.

# trauma --- rh amyg slope --- external_t2
amyg_rh_ext_t2_omit <- amyg_rh %>% #first eliminating NAs 
  filter(!is.na(external_t2))

set.seed(999)
xm <- lm(CSPM_slope_linear ~ trauma + white + inc_needs, data = amyg_rh_ext_t2_omit)  
summary(xm) # traumaTrauma  0.47742    0.13273   3.597 0.000473 ***
xmy <- lm(external_t2 ~ CSPM_slope_linear + trauma + external + white + inc_needs, data = amyg_rh_ext_t2_omit) ####################### added T1
summary(xmy) # CSPM_slope_linear   3.6319     1.3243   2.742  0.00708 ** 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "CSPM_slope_linear", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # yes; ACME p004
plot(m1)

#                 Estimate    95% CI Lower 95% CI Upper p-value   
# ACME             1.7339       0.6087         3.34     0.0042 **
# ADE              2.8127      -1.7182         8.07     0.2360   
# Total Effect     4.5466       0.0214         9.76     0.0542 . 
# Prop. Mediated   0.3814       0.3988        75.08      0.0584 .


# trauma --- insula slope --- external_t2
insula_ext_t2_omit <- insula %>% #first eliminating NAs 
  filter(!is.na(external_t2))

set.seed(999)
xm <- lm(CSPM_slope_linear ~ trauma + white + inc_needs, data = insula_ext_t2_omit)  
summary(xm) #traumaTrauma  0.34237    0.11710   2.924  0.00416 **
xmy <- lm(external_t2 ~ CSPM_slope_linear + trauma + external + white + inc_needs, data = insula_ext_t2_omit) ######################## added T1
summary(xmy) #CSPM_slope_linear  2.74995    1.53112   1.796  0.07511 .  
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "CSPM_slope_linear", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # yes; ACME p05
plot(m1)
#               Estimate 95% CI Lower 95% CI Upper p-value  
#ACME             0.9415       0.0837         2.59   0.053 .
#ADE              3.7756      -0.8991         9.07   0.125  
#Total Effect     4.7171       0.2247         9.96   0.047 *
#Prop. Mediated   0.1996      -0.0377         1.85   0.097 .


##########################################################.
# fMRI mediation model summary ---------------------------
##########################################################.

# Baseline fMRI mediation models ------------------------------------ NONE

# T2 fMRI mediation models ----------------------------------------- BELOW 

# trauma --- rh amyg slope --- external_t2       yes p004
#                 Estimate    95% CI Lower 95% CI Upper p-value   
# ACME             1.7339       0.6087         3.34     0.0042 **
# ADE              2.8127      -1.7182         8.07     0.2360   
# Total Effect     4.5466       0.0214         9.76     0.0542 . 
# Prop. Mediated   0.3814       0.3988        75.08      0.0584 .

# trauma --- insula slope --- external_t2        yes p05
#               Estimate 95% CI Lower 95% CI Upper p-value  
#ACME             0.9415       0.0837         2.59   0.053 .
#ADE              3.7756      -0.8991         9.07   0.125  
#Total Effect     4.7171       0.2247         9.96   0.047 *
#Prop. Mediated   0.1996      -0.0377         1.85   0.097 .






###################################################################.
# gPPI MEDIATIONS
# gPPI right amyg fcMRI with: bl hipp, bl ACC, bl PCC, bl VVS
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
fear_data <- mt_fearcond_data 

hipp_rh_amyg <- hipp_rh_amyg %>%
  left_join(fear_data, by = "subject")

acc_rh_amyg <- acc_rh_amyg %>%
  left_join(fear_data, by = "subject")

pcc_rh_amyg <- pcc_rh_amyg %>%
  left_join(fear_data, by = "subject")

vvs_rh_amyg <- vvs_rh_amyg %>%
  left_join(fear_data, by = "subject")



###################################################################.
# gPPI MEDIATION ----------------- --------------------------------
###################################################################.

# for models where a and b branches are significant or trending

# gPPI baseline mediations -------------------------------------------------------

#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- depression
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg)  
xmy <- lm(cdi ~ gPPI + trauma + white + inc_needs, data = hipp_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # trend; p = 0.082
plot(m1)
#                Estimate    95% CI Lower     95% CI Upper    p-value    
#ACME            0.82972     -0.03007         2.06            0.0820   

#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- panic
hipp_rh_amyg_panic_omit <- hipp_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_panic))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg_panic_omit)  
summary(xm)  #trauma1     -0.44369    0.10201  -4.349  2.6e-05 ***
xmy <- lm(scared_panic ~ gPPI + trauma + white + inc_needs, data = hipp_rh_amyg_panic_omit)
summary(xmy) #gPPI         -1.7276     0.8232  -2.099 0.037637 *  
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # yes; ACME p = 0.0026; prop. mediated estimate 0.2641
plot(m1)

#                  Estimate     95% CI Lower   95% CI Upper p-value   
# ACME             0.7665       0.2532         1.54         0.0026 **
# ADE              2.1359      -0.0601         4.48         0.0610 . 
# Total Effect     2.9025       0.7341         5.27         0.0110 * 
# Prop. Mediated   0.2641       0.1273         3.02         0.0136 * 

#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- gad
hipp_rh_amyg_gad_omit <- hipp_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_gad))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg_gad_omit) 
summary(xm) #trauma1     -0.44369    0.10201  -4.349  2.6e-05 ***
xmy <- lm(scared_gad ~ gPPI + trauma + white + inc_needs, data = hipp_rh_amyg_gad_omit) 
summary(xmy) #gPPI         -1.5865     0.7846  -2.022   0.0451 *
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # ACME p = 0.037
plot(m1)

#                  Estimate     95% CI Lower   95% CI Upper p-value  
# ACME             0.7039       0.0827         1.50         0.037 *
# ADE              1.2959      -0.7916         3.55         0.236  
# Total Effect     1.9999      -0.0320         4.05         0.056 .
# Prop. Mediated   0.3520     -10.6217         0.73         0.091 .

#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- external
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg)  
xmy <- lm(external ~ gPPI + trauma + white + inc_needs, data = hipp_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no; ACME p = 0.8
plot(m1)

#                 Estimate      95% CI Lower  95% CI Upper p-value    
# ACME             0.1687      -0.9401         1.93         0.8    
# ADE             13.0706       8.6344        17.00         <2e-16 ***
# Total Effect    13.2393       9.4737        16.77         <2e-16 ***
# Prop. Mediated   0.0127      -0.0646         0.17         0.8    

#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- ptsd
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg)  
xmy <- lm(ptsd ~ gPPI + trauma + white + inc_needs, data = hipp_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)

# Causal Mediation Analysis 
# Nonparametric Bootstrap Confidence Intervals with the BCa Method
# Simulations: 10000 
#                 Estimate      95% CI Lower 95% CI Upper p-value    
# ACME             -0.727       -2.979         0.69       0.36    
# ADE              27.668       22.315        33.39       <2e-16 ***
# Total Effect     26.941       21.849        32.13       <2e-16 ***
# Prop. Mediated   -0.027       -0.111         0.03       0.36    


#----------------------------------------------------------------------------.
# trauma --- acc_rh_amyg gppi --- depression
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg)  
xmy <- lm(cdi ~ gPPI + trauma + white + inc_needs, data = acc_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME             0.2279      -0.1376         0.97  0.3150    
#ADE              5.9935       2.9364         9.57  0.0002 ***
#  Total Effect     6.2215       3.2425         9.59  <2e-16 ***
#  Prop. Mediated   0.0366      -0.0169         0.21  0.3150    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147

#----------------------------------------------------------------------------.
# trauma --- acc_rh_amyg gppi --- panic
acc_rh_amyg_panic_omit <- acc_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_panic))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg_panic_omit)  
xmy <- lm(scared_panic ~ gPPI + trauma + white + inc_needs, data = acc_rh_amyg_panic_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # trend,  ACME p = 0.088
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME            0.27961      0.00561         0.90   0.088 .
#ADE             2.62287      0.44227         4.99   0.026 *
#  Total Effect    2.90248      0.73407         5.27   0.011 *
#  Prop. Mediated  0.09633      0.01696         5.27   0.098 .
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 145 

#----------------------------------------------------------------------------.
# trauma --- acc_rh_amyg gppi --- GAD
acc_rh_amyg_gad_omit <- acc_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_gad))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg_gad_omit)  
xmy <- lm(scared_gad ~ gPPI + trauma + white + inc_needs, data = acc_rh_amyg_gad_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # trend, ACME p = 0.093
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME            0.277621     0.000838         0.82   0.093 .
#ADE             1.722236    -0.319766         3.78   0.098 .
#Total Effect    1.999857    -0.032015         4.05   0.056 .
#Prop. Mediated  0.138820    -0.029907         2.04   0.143  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 145

#----------------------------------------------------------------------------.
# trauma --- acc_rh_amyg gppi --- external
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg)  
xmy <- lm(external ~ gPPI + trauma + white + inc_needs, data = acc_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#ACME            0.25597     -0.09665         1.35    0.32    
#ADE            12.98336      9.13150        16.51  <2e-16 ***
#  Total Effect   13.23933      9.47367        16.77  <2e-16 ***
#  Prop. Mediated  0.01933     -0.00714         0.11    0.32    



#----------------------------------------------------------------------------.
# trauma --- pcc_rh_amyg gppi --- external
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = pcc_rh_amyg)  
xmy <- lm(external ~ gPPI + trauma + white + inc_needs, data = pcc_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME             0.4872      -0.1869         1.60    0.22    
#ADE             12.7521       8.8861        16.40  <2e-16 ***
#Total Effect    13.2393       9.4737        16.77  <2e-16 ***
#Prop. Mediated   0.0368      -0.0127         0.14    0.22    
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147 

#----------------------------------------------------------------------------.
# trauma --- pcc_rh_amyg gppi --- ptsd
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = pcc_rh_amyg)  
xmy <- lm(ptsd ~ gPPI + trauma + white + inc_needs, data = pcc_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME             0.5036      -0.4332         2.11    0.35    
#ADE             26.4373      21.1682        31.59  <2e-16 ***
#Total Effect    26.9410      21.8486        32.13  <2e-16 ***
#Prop. Mediated   0.0187      -0.0169         0.08    0.35    
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147



#----------------------------------------------------------------------------.
# trauma --- vvs_rh_amyg gppi --- depression
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg)  
xmy <- lm(cdi ~ gPPI + trauma + white + inc_needs, data = vvs_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME            0.85266     -0.05297         2.28   0.104    
#ADE             5.36880      2.47893         8.68   0.001 ***
#  Total Effect    6.22146      3.24245         9.59  <2e-16 ***
#  Prop. Mediated  0.13705     -0.00951         0.38   0.104    
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147

#----------------------------------------------------------------------------.
# trauma --- vvs_rh_amyg gppi --- panic
vvs_rh_amyg_panic_omit <- vvs_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_panic))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg_panic_omit)
summary(xm) #trauma1     -0.56071    0.14509  -3.865 0.000169 ***
xmy <- lm(scared_panic ~ gPPI + trauma + white + inc_needs, data = vvs_rh_amyg_panic_omit) 
summary(xmy) #gPPI         -1.3766     0.5762  -2.389 0.018217 *  
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # yes, ACME p = 0.0042
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value   
#ACME             0.7719       0.2745         1.49  0.0042 **
#ADE              2.1306      -0.0237         4.46  0.0570 . 
#Total Effect     2.9025       0.7341         5.27  0.0110 * 
#Prop. Mediated   0.2659       0.1400         5.57  0.0152 * 
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 145 


#----------------------------------------------------------------------------.
# trauma --- vvs_rh_amyg gppi --- gad
vvs_rh_amyg_gad_omit <- vvs_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_gad))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg_gad_omit)  
summary(xm) #trauma1     -0.56071    0.14509  -3.865 0.000169 ***
xmy <- lm(scared_gad ~ gPPI + trauma + white + inc_needs, data = vvs_rh_amyg_gad_omit) 
summary(xmy) #gPPI         -1.3064     0.5487  -2.381   0.0186 *
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # yes, ACME = 0.014
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME              0.733        0.199         1.49   0.014 *
#ADE               1.267       -0.721         3.42   0.221  
#Total Effect      2.000       -0.032         4.05   0.056 .
#Prop. Mediated    0.366      -31.645         0.49   0.069 .
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 145 


#----------------------------------------------------------------------------.
# trauma --- vvs_rh_amyg gppi --- external
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg)  
xmy <- lm(external ~ gPPI + trauma + white + inc_needs, data = vvs_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME             0.1992      -0.6473         1.72    0.71    
#ADE             13.0401       8.9809        16.60  <2e-16 ***
#Total Effect    13.2393       9.4737        16.77  <2e-16 ***
#Prop. Mediated   0.0150      -0.0511         0.14    0.71    
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147 

#----------------------------------------------------------------------------.
# trauma --- vvs_rh_amyg gppi --- ptsd
set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg)  
xmy <- lm(ptsd ~ gPPI + trauma + white + inc_needs, data = vvs_rh_amyg) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value    
#ACME             0.3594      -1.1075         2.26    0.63    
#ADE             26.5815      21.0704        32.18  <2e-16 ***
#Total Effect    26.9410      21.8486        32.13  <2e-16 ***
#Prop. Mediated   0.0133      -0.0393         0.09    0.63    
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 147 


# gPPI T2 mediations -------------------------------------------------------


#----------------------------------------------------------------------------.
# trauma --- hipp_rh_amyg gppi --- ext_t2 (B BRANCH IS TREND)
hipp_rh_amyg_ext_t2_omit <- hipp_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(external_t2))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = hipp_rh_amyg_ext_t2_omit)  
xmy <- lm(external_t2 ~ gPPI + trauma + external + white + inc_needs, data = hipp_rh_amyg_ext_t2_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)

#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME              0.864       -0.588         3.01    0.28  
#ADE               3.758       -0.968         9.31    0.13  
#Total Effect      4.622        0.160         9.87    0.05 *
#  Prop. Mediated    0.187       -0.240         2.11    0.32  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 121 



# trauma --- acc_rh_amyg gppi --- ptsd_t2 
acc_rh_amyg_ptsd_t2_omit <- acc_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(ptsd_t2))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = acc_rh_amyg_ptsd_t2_omit)  
xmy <- lm(ptsd_t2 ~ gPPI + trauma + ptsd + white + inc_needs, data = acc_rh_amyg_ptsd_t2_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)

#Estimate 95% CI Lower 95% CI Upper p-value
#ACME            -0.1809      -1.7286          0.3    0.59
#ADE              4.0361      -3.7340         11.6    0.28
#Total Effect     3.8552      -3.8993         11.3    0.31
#Prop. Mediated  -0.0469       0.0147         82.5    0.70
#Sample Size Used: 121 


# trauma --- pcc_rh_amyg gppi --- ext_t2 (B BRANCH IS TREND)
pcc_rh_amyg_ext_t2_omit <- pcc_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(external_t2))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = pcc_rh_amyg_ext_t2_omit)  
xmy <- lm(external_t2 ~ gPPI + trauma + external + white + inc_needs, data = pcc_rh_amyg_ext_t2_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME           -0.68723     -2.06431         0.17   0.177  
#ADE             5.13440      0.29994        10.34   0.040 *
#  Total Effect    4.44717     -0.00205         9.50   0.055 .
#Prop. Mediated -0.15453     -1.10451         0.16   0.215  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 121 


# trauma --- vvs_rh_amyg gppi --- panic_t2 (B BRANCH IS TREND)
vvs_rh_amyg_panic_t2_omit <- vvs_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(scared_panic_t2)) %>%
  filter(!is.na(scared_panic))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg_panic_t2_omit)  
xmy <- lm(scared_panic_t2 ~ gPPI + trauma + scared_panic + white + inc_needs, data = vvs_rh_amyg_panic_t2_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value
#ACME              0.1589      -0.4269         0.75    0.56
#ADE              -0.1698      -1.6396         1.18    0.79
#Total Effect     -0.0109      -1.5819         1.33    0.96
#Prop. Mediated  -14.5848    -107.4107         1.10    0.87
#Sample Size Used: 118 

# trauma --- vvs_rh_amyg gppi --- ext_t2 (B BRANCH IS TREND)
vvs_rh_amyg_ext_t2_omit <- vvs_rh_amyg %>% #first eliminating NAs 
  filter(!is.na(external_t2))

set.seed(999)
xm <- lm(gPPI ~ trauma + white + inc_needs, data = vvs_rh_amyg_ext_t2_omit)  
xmy <- lm(external_t2 ~ gPPI + trauma + external + white + inc_needs, data = vvs_rh_amyg_ext_t2_omit) 
m1 <- mediate(xm, xmy, treat = "trauma", mediator = "gPPI", 
              sims = 10000, boot = TRUE, boot.ci.type = "bca")    
summary(m1) # no
plot(m1)
#Estimate 95% CI Lower 95% CI Upper p-value  
#ACME             0.4497      -0.6412         2.48   0.503  
#ADE              4.1804      -0.5525         9.68   0.094 .
#Total Effect     4.6301       0.1689         9.90   0.048 *
#  Prop. Mediated   0.0971      -0.1291         3.15   0.525  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Sample Size Used: 121 




##########################################################.
# gPPI mediation model summary ---------------------------
##########################################################.

######## BASELINE MED MODELS

# trauma --- rh amyg-hipp gppi --- depression   (ACME TREND 0.08)
# trauma --- rh amyg-hipp gppi --- panic        (ACME 0.003)
# trauma --- rh amyg-hipp gppi --- GAD          (ACME 0.04)

# trauma --- rh amyg-acc gppi --- panic         (ACME TREND 0.09)
# trauma --- rh amyg-acc gppi --- GAD           (ACME TREND 0.09) 

# trauma --- rh amyg-vvs gppi --- panic         (ACME 0.004)
# trauma --- rh amyg-vvs gppi --- GAD           (ACME 0.01)


######## T2 MED MODELS
# none significant



