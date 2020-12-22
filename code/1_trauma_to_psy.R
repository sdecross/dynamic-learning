####################################################################.
# SOURCE AND LOAD DATA
####################################################################.
rm(list=ls())

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData") # object "mt_fearcond_data"
mt_fearcond_data <- data
if (!require("pacman")) {install.packages("pacman"); require("pacman")}
p_load("MASS")


####################################################################.
# Visualizing outcome variable
####################################################################.
# all are right-skewed in whole sample except external,
# therefore expect neg bin distribs to be best fit for everything except externalizing
hist(mt_fearcond_data$cdi, breaks = 20, xlab = "CDI counts", main = "Histogram CDI") # right-skewed
hist(mt_fearcond_data$scared_panic, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC") # right-skewed
hist(mt_fearcond_data$scared_gad, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD") # right-skewed
hist(mt_fearcond_data$external, breaks = 20, xlab = "External counts", main = "Histogram external") ## THIS ONE IS NORMAL DISTRIB
hist(mt_fearcond_data$ptsd, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD sev") # right-skewed

hist(mt_fearcond_data$cdi_t2, breaks = 20, xlab = "CDI counts", main = "Histogram CDI_T2") # right-skewed
hist(mt_fearcond_data$scared_panic_t2, breaks = 20, xlab = "SCARED_PANIC counts", main = "Histogram SCARED_PANIC_T2") # right-skewed
hist(mt_fearcond_data$scared_gad_t2, breaks = 20, xlab = "SCARED_GAD counts", main = "Histogram SCARED_GAD_T2") # right-skewed
hist(mt_fearcond_data$external_t2, breaks = 20, xlab = "External counts", main = "Histogram external_t2") ## THIS ONE IS NORMAL DISTRIB
hist(mt_fearcond_data$ptsd_t2, breaks = 20, xlab = "PTSD sev counts", main = "Histogram PTSD_T2") # right-skewed


########################################################################.
# Does trauma (assigned group; binary 0/1) predict psychopathology? ####
########################################################################.

# CDI - depression
plot(mt_fearcond_data$trauma, mt_fearcond_data$cdi)
reg_cdi_nor <- lm(cdi ~ trauma + white + inc_needs, data = mt_fearcond_data)
reg_cdi_nb <- glm.nb(cdi ~ trauma + white + inc_needs, data = mt_fearcond_data)
AIC(reg_cdi_nor, reg_cdi_nb)
summary(reg_cdi_nb) # yes, p < 0.001

# SCARED_PANIC - panic 
plot(mt_fearcond_data$trauma, mt_fearcond_data$scared_panic)
reg_scared_panic_nor <- lm(scared_panic ~ trauma + white + inc_needs, data = mt_fearcond_data)
reg_scared_panic_nb <- glm.nb(scared_panic ~ trauma + white + inc_needs, data = mt_fearcond_data)
AIC(reg_scared_panic_nor, reg_scared_panic_nb)
summary(reg_scared_panic_nb) # yes, p = 0.003

# SCARED_GAD - generalized anxiety
plot(mt_fearcond_data$trauma, mt_fearcond_data$scared_gad)
reg_scared_gad_nor <- lm(scared_gad ~ trauma + white + inc_needs, data = mt_fearcond_data)
reg_scared_gad_nb <- glm.nb(scared_gad ~ trauma + white + inc_needs, data = mt_fearcond_data)
AIC(reg_scared_gad_nor, reg_scared_gad_nb)
summary(reg_scared_gad_nb) # no but trend, p = 0.083

# EXTERNAL - externalizing
plot(mt_fearcond_data$trauma, mt_fearcond_data$external)
reg_external_nor <- lm(external ~ trauma + white + inc_needs, data = mt_fearcond_data)
reg_external_nb <- glm.nb(external ~ trauma + white + inc_needs, data = mt_fearcond_data)
AIC(reg_external_nor, reg_external_nb)
exp((1020.562-1022.035)/2) #  47% chance that normal is best fitting; 53% chance negbin best fitting
summary(reg_external_nb) # go with negbin; yes, p < 0.001
summary(reg_external_nor) # other one has same results

# ptsd - PTSD
plot(mt_fearcond_data$trauma, mt_fearcond_data$ptsd)
reg_ptsd_nor <- lm(ptsd ~ trauma + white + inc_needs, data = mt_fearcond_data)
reg_ptsd_nb <- glm.nb(ptsd ~ trauma + white + inc_needs, data = mt_fearcond_data)
AIC(reg_ptsd_nor, reg_ptsd_nb)
summary(reg_ptsd_nb) # yes, p < 0.001


##########################################################################.
# Does trauma (assigned group; binary 0/1) predict T2 psychopathology ####
# at T2, controlling for symptoms at T1?
##########################################################################.

# CDI - depression
plot(mt_fearcond_data$trauma, mt_fearcond_data$cdi_t2)
reg_cdi_nor_t2 <- lm(cdi_t2 ~ trauma + cdi + white + inc_needs, data = mt_fearcond_data)
reg_cdi_nb_t2 <- glm.nb(cdi_t2 ~ trauma + cdi + white + inc_needs, data = mt_fearcond_data)
AIC(reg_cdi_nor_t2, reg_cdi_nb_t2)
summary(reg_cdi_nb_t2) # no

# SCARED_PANIC - panic 
plot(mt_fearcond_data$trauma, mt_fearcond_data$scared_panic_t2)
reg_scared_panic_nor_t2 <- lm(scared_panic_t2 ~ trauma + scared_panic + white + inc_needs, data = mt_fearcond_data)
reg_scared_panic_nb_t2 <- glm.nb(scared_panic_t2 ~ trauma + scared_panic + white + inc_needs, data = mt_fearcond_data)
AIC(reg_scared_panic_nor_t2, reg_scared_panic_nb_t2)
summary(reg_scared_panic_nb_t2) # no

# SCARED_GAD - generalized anxiety
plot(mt_fearcond_data$trauma, mt_fearcond_data$scared_gad_t2)
reg_scared_gad_nor_t2 <- lm(scared_gad_t2 ~ trauma + scared_gad + white + inc_needs, data = mt_fearcond_data)
reg_scared_gad_nb_t2 <- glm.nb(scared_gad_t2 ~ trauma + scared_gad + white + inc_needs, data = mt_fearcond_data)
AIC(reg_scared_gad_nor_t2, reg_scared_gad_nb_t2)
summary(reg_scared_gad_nb_t2) # no

# EXTERNAL - externalizing
plot(mt_fearcond_data$trauma, mt_fearcond_data$external_t2)
reg_external_nor_t2 <- lm(external_t2 ~ trauma + external + white + inc_needs, data = mt_fearcond_data)
reg_external_nb_t2 <- glm.nb(external_t2 ~ trauma + external + white + inc_needs, data = mt_fearcond_data)
AIC(reg_external_nor_t2, reg_external_nb_t2)
exp((852.5761-854.8181)/2) #  32% chance that normal is best fitting
summary(reg_external_nb_t2) # 0.038

# ptsd - PTSD
plot(mt_fearcond_data$trauma, mt_fearcond_data$ptsd_t2)
reg_ptsd_nor_t2 <- lm(ptsd_t2 ~ trauma + ptsd + white + inc_needs, data = mt_fearcond_data)
reg_ptsd_nb_t2 <- glm.nb(ptsd_t2 ~ trauma + ptsd + white + inc_needs, data = mt_fearcond_data)
AIC(reg_ptsd_nor_t2, reg_ptsd_nb_t2)
summary(reg_ptsd_nb_t2) # 0.042




####################################################################.
# SUMMARY
####################################################################.

#-------------------------------------------------------------------.
# DOES TRAUMA (GROUP) PREDICT PSYCHOPATHOLOGY?

# Negbin models best fit for all psychopathology (external was close)

# Trauma predicts cdi, scared_panic, external, ptsd
# Trend for scared_gad (p = 0.08)


#-------------------------------------------------------------------.
# DOES TRAUMA (GROUP) PREDICT PSYCHOPATHOLOGY
# AT FOLLOW-UP, CONTROLLING FOR BASELINE SYMPTOMS?

# Yes, trauma (controlling for baseline symptoms) predicts
# externalizing (p = 0.038) and PTSD (p = 0.042) at follow-up

# Negbin models best fit


