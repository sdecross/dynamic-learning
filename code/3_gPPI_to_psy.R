###################################################################.
# gPPI to psy
# gPPI amyg rh fcMRI with: bl hipp, bl ACC, bl PCC, bl VVS
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

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "mdata"
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
# GPPI TO PSY -----------------------------------------------------
###################################################################.

#------------------------------------------------------------------.
# VISUALIZING OUTCOME VAR 

# baseline symptoms
# except for externalizing symptoms, DVS are all metric variables with lower bounds at 0 with no upper bound 
# >> neg bin will prob fit betterfor all but externalizing
hist(hipp_rh_amyg$cdi, breaks = 20, xlab = "CDI counts", main = "Histogram CDI") # right-skewed
hist(hipp_rh_amyg$scared_panic, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC") # right-skewed
hist(hipp_rh_amyg$scared_gad, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD") # right-skewed
hist(hipp_rh_amyg$external, breaks = 25, xlab = "External counts", main = "Histogram external") # normal
hist(hipp_rh_amyg$ptsd, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD sev") # right-skewed

# T2 symptoms
# my DVS are all metric variables with lower bounds at 0 with no upper bound
# expect negbin versions to fit better for all models, except possibly external 
hist(hipp_rh_amyg$cdi_t2, breaks = 20, xlab = "CDI counts", main = "Histogram CDI") # right-skewed
hist(hipp_rh_amyg$scared_panic_t2, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC") # right-skewed
hist(hipp_rh_amyg$scared_gad_t2, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD") # right-skewed
hist(hipp_rh_amyg$external_t2, breaks = 25, xlab = "External counts", main = "Histogram external") # possibly normal
hist(hipp_rh_amyg$ptsd_t2, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD sev") # right-skewed

# initial look at gPPI data
ggstatsplot::ggbetweenstats(data = hipp_rh_amyg, x = trauma, y = gPPI, messages = FALSE, 
                            xlab = "Group", ylab = "gPPI interaction extracted zstat", title = "hipp_rh_amyg gPPI CS+ vs. CS-") 

ggstatsplot::ggbetweenstats(data = acc_rh_amyg, x = trauma, y = gPPI, messages = FALSE, 
                            xlab = "Group", ylab = "gPPI interaction extracted zstat", title = "acc_rh_amyg gPPI CS+ vs. CS-") 

ggstatsplot::ggbetweenstats(data = pcc_rh_amyg, x = trauma, y = gPPI, messages = FALSE, 
                            xlab = "Group", ylab = "gPPI interaction extracted zstat", title = "pcc_rh_amyg gPPI CS+ vs. CS-") 

ggstatsplot::ggbetweenstats(data = vvs_rh_amyg, x = trauma, y = gPPI, messages = FALSE, 
                            xlab = "Group", ylab = "gPPI interaction extracted zstat", title = "vvs_rh_amyg gPPI CS+ vs. CS-") 




#------------------------------------------------------------------.
# regressions with baseline symptoms ------------------

# hipp_amyg-rh 

# cdi
scatterplot(cdi ~ hipp_rh_amyg$gPPI, data = hipp_rh_amyg, id = list(n = 5)) 
m1 <- lm(cdi ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(cdi ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.00159, -0.44585

# scared_panic
scatterplot(scared_panic ~ hipp_rh_amyg$gPPI, data = hipp_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_panic ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(scared_panic ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.00097, b= -0.51407

# scared_gad
scatterplot(scared_gad ~ hipp_rh_amyg$gPPI, data = hipp_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_gad ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(scared_gad ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.0103, b = -0.42962

# external
scatterplot(external ~ hipp_rh_amyg$gPPI, data = hipp_rh_amyg, id = list(n = 5)) 
m1 <- lm(external ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(external ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit 
exp((1068.595-1069.625)/2) # 59% normal is better
summary(m1) # p = 0.01187, b= -3.8625

# pstd
scatterplot(ptsd ~ hipp_rh_amyg$gPPI, data = hipp_rh_amyg, id = list(n = 5)) 
m1 <- lm(ptsd ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(ptsd ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.017074 b= -0.62304


# acc_amyg-rh 

# cdi
scatterplot(cdi ~ acc_rh_amyg$gPPI, data = acc_rh_amyg, id = list(n = 5)) 
m1 <- lm(cdi ~ gPPI + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(cdi ~ gPPI + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.107

# scared_panic
scatterplot(scared_panic ~ acc_rh_amyg$gPPI, data = acc_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_panic ~ gPPI + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(scared_panic ~ gPPI + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.021, b = 0.28699

# scared_gad
scatterplot(scared_gad ~ acc_rh_amyg$gPPI, data = acc_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_gad ~ gPPI + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(scared_gad ~ gPPI + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.0696

# external
scatterplot(external ~ acc_rh_amyg$gPPI, data = acc_rh_amyg, id = list(n = 5)) 
m1 <- lm(external ~ gPPI + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(external ~ gPPI + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit 
exp((1071.308-1072.57)/2) # 53% normal is better
summary(m1) # p = 0.062

# pstd
scatterplot(ptsd ~ acc_rh_amyg$gPPI, data = acc_rh_amyg, id = list(n = 5)) 
m1 <- lm(ptsd ~ gPPI + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(ptsd ~ gPPI + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no


# pcc_amyg-rh 

# cdi
scatterplot(cdi ~ pcc_rh_amyg$gPPI, data = pcc_rh_amyg, id = list(n = 5)) 
m1 <- lm(cdi ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(cdi ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_panic
scatterplot(scared_panic ~ pcc_rh_amyg$gPPI, data = pcc_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_panic ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(scared_panic ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_gad
scatterplot(scared_gad ~ pcc_rh_amyg$gPPI, data = pcc_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_gad ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(scared_gad ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# external
scatterplot(external ~ pcc_rh_amyg$gPPI, data = pcc_rh_amyg, id = list(n = 5)) 
m1 <- lm(external ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(external ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit 
exp((1068.014-1069.155)/2) # 69% normal is better
summary(m1) # p = 0.009, b = -2.6575

# pstd
scatterplot(ptsd ~ pcc_rh_amyg$gPPI, data = pcc_rh_amyg, id = list(n = 5)) 
m1 <- lm(ptsd ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(ptsd ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.066

# vvs_amyg-rh 

# cdi
scatterplot(cdi ~ vvs_rh_amyg$gPPI, data = vvs_rh_amyg, id = list(n = 5)) 
m1 <- lm(cdi ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(cdi ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.003, b = -0.29475

# scared_panic
scatterplot(scared_panic ~ vvs_rh_amyg$gPPI, data = vvs_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_panic ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(scared_panic ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.000553, b = -0.38092

# scared_gad
scatterplot(scared_gad ~ vvs_rh_amyg$gPPI, data = vvs_rh_amyg, id = list(n = 5)) 
m1 <- lm(scared_gad ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(scared_gad ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.0037, b = -0.34411

# external
scatterplot(external ~ vvs_rh_amyg$gPPI, data = vvs_rh_amyg, id = list(n = 5)) 
m1 <- lm(external ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(external ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit 
exp((1069.698-1070.67)/2) # 61% normal is better
summary(m1) # p = 0.02113, b = -2.5072

# pstd
scatterplot(ptsd ~ vvs_rh_amyg$gPPI, data = vvs_rh_amyg, id = list(n = 5)) 
m1 <- lm(ptsd ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(ptsd ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.034798, b = -0.38893

#------------------------------------------------------------------.
# regressions with T2 symptoms ------------------------

# hipp_amyg-rh -----------------------------------------------------.

# cdi
m1 <- lm(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_panic
m1 <- lm(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_gad
m1 <- lm(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

#external
m1 <- lm(external_t2 ~ gPPI + external + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(external_t2 ~ gPPI + external + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
exp((853.7594-856.0746)/2) # 31% chance normal is better, go with negbin
summary(m2) # 0.0793

# pstd
m1 <- lm(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = hipp_rh_amyg)
m2 <- glm.nb(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = hipp_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no


# acc_amyg-rh -----------------------------------------------------.

# cdi
m1 <- lm(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_panic
m1 <- lm(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_gad
m1 <- lm(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

#external
m1 <- lm(external_t2 ~ gPPI + external + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(external_t2 ~ gPPI + external + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
exp((856.7575-858.9657)/2) # 33% chance normal is better, go with negbin
summary(m2) # no

# pstd
m1 <- lm(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = acc_rh_amyg)
m2 <- glm.nb(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = acc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.0468, b =-0.45811

# pcc_amyg-rh -----------------------------------------------------.

# cdi
m1 <- lm(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_panic
m1 <- lm(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_gad
m1 <- lm(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

#external
m1 <- lm(external_t2 ~ gPPI + external + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(external_t2 ~ gPPI + external + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
exp((855.0848-857.3946)/2) # 31% chance normal is better, go with negbin
summary(m2) # 0.0793

# pstd
m1 <- lm(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = pcc_rh_amyg)
m2 <- glm.nb(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = pcc_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# vvs_amyg-rh -----------------------------------------------------.

# cdi
m1 <- lm(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(cdi_t2 ~ gPPI + cdi + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# scared_panic
m1 <- lm(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(scared_panic_t2 ~ gPPI + scared_panic + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # p = 0.0946

# scared_gad
m1 <- lm(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(scared_gad_t2 ~ gPPI + scared_gad + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

#external
m1 <- lm(external_t2 ~ gPPI + external + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(external_t2 ~ gPPI + external + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
exp((855.4379-857.5904)/2) # 34% chance normal is better, go with negbin
summary(m2) # 0.0793

# pstd
m1 <- lm(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = vvs_rh_amyg)
m2 <- glm.nb(ptsd_t2 ~ gPPI + ptsd + white + inc_needs, data = vvs_rh_amyg)
AIC(m1, m2) # m2 better fit
summary(m2) # no

# added variable plots of significant results ------------------------------------


############################################################ .
# gPPI-PSY SIGNIFICANT REGRESSION PLOTS ---------------------.
############################################################ .
# Added variable plots (aka partial regression plots, adjusted variable plots, individual coefficient plots), 
# are “results” plots. They are plots showing the estimated relationship between the response and an explanatory 
# variable after accounting for the other variables in the model. Use these to plot the findings, 
# holding covariates white and inc_needs constant (plotting for the median income, white participant).
# https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
# If you plot like this it doesn't account for the influence of race or inc_needs on the model:
# ggplot(brain_data, aes(x = brain_data, y = symptom)) + geom_point() + theme_classic() +
# labs(x = "Brain data", y = "Symptom") + geom_smooth(method = MASS::glm.nb, se = TRUE)

# RE-RUN CORRECT MODELS FIRST; NAMES OVERWRITE IN EACH ANALYSIS SECTION

# SIGNIFICNAT MODELS ONLY (NOT TRENDS)

# baseline ---------------------------------------------------.

# hipp-amyg-rh ----------------------------------------------.
# cdi # p = 0.00159, b = -0.44585
newdata <- data.frame(gPPI = hipp_rh_amyg$gPPI, white = factor(1), inc_needs = median(hipp_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  cdi <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, cdi)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = hipp_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral hippocampus connectivity", y = "Depression")
#ggsave("plots/symptom-regressions/gPPI/hipp_amyg-rh_cdi.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# scared_panic # p = 0.00097, b= -0.51407
newdata <- data.frame(gPPI = hipp_rh_amyg$gPPI, white = factor(1), inc_needs = median(hipp_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  scared_panic <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, scared_panic)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = hipp_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral hippocampus connectivity", y = "Panic")
#ggsave("plots/symptom-regressions/gPPI/hipp_amyg-rh_panic.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# scared_gad # p = 0.0103, b = -0.42962
newdata <- data.frame(gPPI = hipp_rh_amyg$gPPI, white = factor(1), inc_needs = median(hipp_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  scared_gad <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, scared_gad)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = hipp_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral hippocampus connectivity", y = "Generalized Anxiety")
#ggsave("plots/symptom-regressions/gPPI/hipp_amyg-rh_gad.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# ptsd # p = 0.017074 b= -0.62304
newdata <- data.frame(gPPI = hipp_rh_amyg$gPPI, white = factor(1), inc_needs = median(hipp_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  ptsd <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, ptsd)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = hipp_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral hippocampus connectivity", y = "PTSD")
#ggsave("plots/symptom-regressions/gPPI/hipp_amyg-rh_ptsd.pdf", width = 4.5, height = 4, units = "in", dpi = 600)


# acc-amyg-rh ----------------------------------------------.

# scared_panic # p = 0.021, b = 0.28699
newdata <- data.frame(gPPI = acc_rh_amyg$gPPI, white = factor(1), inc_needs = median(acc_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  scared_panic <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, scared_panic)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = acc_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-anterior cingulate cortex connectivity", y = "Panic")
#ggsave("plots/symptom-regressions/gPPI/acc_amyg-rh_panic.pdf", width = 4.5, height = 4, units = "in", dpi = 600)


# pcc-amyg-rh ----------------------------------------------.
# none

# vvs-amyg-rh ----------------------------------------------.

# cdi # p = 0.003, b = -0.29475
newdata <- data.frame(gPPI = vvs_rh_amyg$gPPI, white = factor(1), inc_needs = median(vvs_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  cdi <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, cdi)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = vvs_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior parahippocampal gyrus connectivity", y = "Depression")
#ggsave("plots/symptom-regressions/gPPI/vvs_amyg-rh_cdi.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# scared_panic # p = 0.000553, b = -0.38092
newdata <- data.frame(gPPI = vvs_rh_amyg$gPPI, white = factor(1), inc_needs = median(vvs_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  scared_panic <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, scared_panic)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = vvs_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior parahippocampal gyrus connectivity", y = "Panic")
#ggsave("plots/symptom-regressions/gPPI/vvs_amyg-rh_panic.pdf", width = 4.5, height = 4, units = "in", dpi = 600)


# scared_gad # p = 0.0037, b = -0.34411
newdata <- data.frame(gPPI = vvs_rh_amyg$gPPI, white = factor(1), inc_needs = median(vvs_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  scared_gad <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, scared_gad)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = vvs_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior parahippocampal gyrus connectivity", y = "Generalized Anxiety")
#ggsave("plots/symptom-regressions/gPPI/vvs_amyg-rh_gad.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# ptsd # p = 0.034798, b = -0.38893
#m2 <- glm.nb(ptsd ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
newdata <- data.frame(gPPI = vvs_rh_amyg$gPPI, white = factor(1), inc_needs = median(vvs_rh_amyg$inc_needs))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  ptsd <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, ptsd)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = vvs_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior parahippocampal gyrus connectivity", y = "PTSD")
#ggsave("plots/symptom-regressions/gPPI/vvs_amyg-rh_ptsd.pdf", width = 4.5, height = 4, units = "in", dpi = 600)



################################################################################.
################################################################################.
# FOR NORMAL MODELS  ############.
# m1... doesn't work the same way!

# hipp-amyg, external # p = 0.01187, b= -3.8625 ############################################################.

#m1 <- lm(external ~ gPPI + white + inc_needs, data = hipp_rh_amyg)
#summary(m1)
m1 <- glm(external ~ gPPI + white + inc_needs, data = hipp_rh_amyg) #just making as glm object
summary(m1)

newdata <- data.frame(gPPI = hipp_rh_amyg$gPPI, white = factor(1), inc_needs = median(hipp_rh_amyg$inc_needs))
newdata$p.hat <- predict(m1, newdata, type = "response")
newdata <- cbind(newdata, predict(m1, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  external <- fit
  LL <- fit - 1.96 * se.fit
  UL <- fit + 1.96 * se.fit
})

ggplot(newdata, aes(gPPI, external)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = hipp_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral hippocampus connectivity", y = "Externalizing")
#ggsave("plots/symptom-regressions/gPPI/hipp_amyg-rh_external.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# pcc-amyg, external # p = 0.009, b = -2.6575 
#m1 <- lm(external ~ gPPI + white + inc_needs, data = pcc_rh_amyg)
m1 <- glm(external ~ gPPI + white + inc_needs, data = pcc_rh_amyg)

newdata <- data.frame(gPPI = pcc_rh_amyg$gPPI, white = factor(1), inc_needs = median(pcc_rh_amyg$inc_needs))
newdata$p.hat <- predict(m1, newdata, type = "response")
newdata <- cbind(newdata, predict(m1, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  external <- fit
  LL <- fit - 1.96 * se.fit
  UL <- fit + 1.96 * se.fit
})

ggplot(newdata, aes(gPPI, external)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = pcc_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior cingulate cortex connectivity", y = "Externalizing")
#ggsave("plots/symptom-regressions/gPPI/pcc_amyg-rh_external.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# vvs-amyg, external # p = 0.02113, b = -2.5072
#m1 <- lm(external ~ gPPI + white + inc_needs, data = vvs_rh_amyg)
m1 <- glm(external ~ gPPI + white + inc_needs, data = vvs_rh_amyg)

newdata <- data.frame(gPPI = vvs_rh_amyg$gPPI, white = factor(1), inc_needs = median(vvs_rh_amyg$inc_needs))
newdata$p.hat <- predict(m1, newdata, type = "response")
newdata <- cbind(newdata, predict(m1, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  external <- fit
  LL <- fit - 1.96 * se.fit
  UL <- fit + 1.96 * se.fit
})

ggplot(newdata, aes(gPPI, external)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = vvs_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-bilateral posterior parahippocampal gyrus connectivity", y = "Externalizing")
#ggsave("plots/symptom-regressions/gPPI/vvs_amyg-rh_external.pdf", width = 4.5, height = 4, units = "in", dpi = 600)

# T2 ---------------------------------------------------------.
# also control for the baseline sxs

# acc_amyg-rh ---------------------------------------------------.
# ptsd # p = 0.0468, b =-0.45811
newdata <- data.frame(gPPI = acc_rh_amyg$gPPI, white = factor(1), inc_needs = median(acc_rh_amyg$inc_needs), ptsd = mean(acc_rh_amyg$ptsd, na.rm = T))
newdata$p.hat <- predict(m2, newdata, type = "response")
newdata <- cbind(newdata, predict(m2, newdata, type = "link", se.fit = TRUE))
newdata <- within(newdata, {
  ptsd_t2 <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(gPPI, ptsd_t2)) + 
  theme_classic() +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.25) +
  geom_line(colour = "#3366FF", size = 1) + 
  geom_point(data = acc_rh_amyg, colour = 'black') +
  labs(x = "Right amygdala-anterior cingulate cortex connectivity", y = "PTSD at 2-year follow-up")
#ggsave("plots/symptom-regressions/gPPI/acc_amyg-rh_ptsd_t2.pdf", width = 4.5, height = 4, units = "in", dpi = 600)



