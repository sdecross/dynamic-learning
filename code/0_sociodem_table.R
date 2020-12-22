############################################################
# SOURCE & LOAD DATA
############################################################

rm(list=ls())

if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}

load("data/behavioral/FearLearning_n147_Demogs_Symptoms_T1_T2.RData")


############################################################
# SOCIODEMS FOR TABLE
############################################################

# n control (70) vs. trauma (77) #table(data$trauma)
n <- data %>% 
  group_by(trauma) %>%
  summarise(Count = n())
n

# age ------------------------------------------------------------------------------------------------
hist(data$age, breaks = 20, xlab = "Age", main = "Histogram")
df1 <- data.frame(trauma=c('control','trauma'),
                       mean=tapply(data$age, data$trauma, mean),
                       n=tapply(data$age, data$trauma, length),
                       sd=tapply(data$age, data$trauma, sd))
df1

t.test(age ~ trauma, data = data) #p-value = 0.5014

# sex ------------------------------------------------------------------------------------------------
with(data, table(trauma, sex))               

sex_test <- table(data$trauma, data$sex) # create table with the needed variables and do chi-square
chisq.test(sex_test) #p-value = 0.6768

#library(pastecs)
#stat.desc(data$age)

# race & white ------------------------------------------------------------------------------------------------
with(data, table(trauma, race)) 
race_test <- table(data$trauma, data$race) # create table with the needed variables and do chi-square
chisq.test(race_test) #p-value = 1.383e-08

# race = 1 is white; ref group
# 2 is black
# 3 is latino
# 4 is asian or pacific islander
# 5 is biracial or other

# for controls:
47/70 #white
4/70 #black
11/70 #asian
5/70 #latino
3/70 #other/biracial

# for trauma:
18/77 #white
31/77 #black
6/77 #asian
9/77 #latino
13/77 #other/biracial



with(data, table(trauma, white))
white_test <- table(data$trauma, data$white)  # create table with the needed variables and do chi-square
chisq.test(white_test) #p-value = 2.342e-07

# inc_needs ------------------------------------------------------------------------------------------------
df2 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$inc_needs, data$trauma, mean),
                  n=tapply(data$inc_needs, data$trauma, length),
                  sd=tapply(data$inc_needs, data$trauma, sd))
df2
t.test(inc_needs ~ trauma, data = data) #p-value < 2.2e-16


# cdi ------------------------------------------------------------------------------------------------
df3 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$cdi, data$trauma, mean),
                  n=tapply(data$cdi, data$trauma, length),
                  sd=tapply(data$cdi, data$trauma, sd))
df3
t.test(cdi ~ trauma, data = data) #p-value = 1.693e-07

# scared_panic ------------------------------------------------------------------------------------------------
df4 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$scared_panic, data$trauma, mean, na.rm=TRUE),
                  sd=tapply(data$scared_panic, data$trauma, sd, na.rm=TRUE))
df4

tmp <- data %>% 
  filter(!is.na(scared_panic))
tmp_df <- data.frame(trauma=c('control','trauma'),
                     length=tapply(tmp$scared_panic, tmp$trauma, length))
tmp_df # shows total number that have data
n$Count[1] - tmp_df$length[1] # number inc_needs missing for control (0)
n$Count[2] - tmp_df$length[2] # number inc_needs missing for trauma (2)

t.test(scared_panic ~ trauma, data = data, na.rm=TRUE) #p-value = 3.944e-05


# scared_gad ------------------------------------------------------------------------------------------------
df5 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$scared_gad, data$trauma, mean, na.rm=TRUE),
                  sd=tapply(data$scared_gad, data$trauma, sd, na.rm=TRUE))
df5

tmp <- data %>% 
  filter(!is.na(scared_gad))
tmp_df <- data.frame(trauma=c('control','trauma'),
                     length=tapply(tmp$scared_gad, tmp$trauma, length))
tmp_df # shows total number that have data
n$Count[1] - tmp_df$length[1] # number inc_needs missing for control (0)
n$Count[2] - tmp_df$length[2] # number inc_needs missing for trauma (2)

t.test(scared_gad ~ trauma, data = data, na.rm=TRUE) #p-value = 0.0008868


# external ------------------------------------------------------------------------------------------------
df6 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$external, data$trauma, mean),
                  n=tapply(data$external, data$trauma, length),
                  sd=tapply(data$external, data$trauma, sd))
df6
t.test(external ~ trauma, data = data)

# ptsd ------------------------------------------------------------------------------------------------
df7 <- data.frame(trauma=c('control','trauma'),
                  mean=tapply(data$ptsd, data$trauma, mean),
                  n=tapply(data$ptsd, data$trauma, length),
                  sd=tapply(data$ptsd, data$trauma, sd))
df7
t.test(ptsd ~ trauma, data = data)




