# Biostatistics 624 Project
# Created by: Taylor Martin 2018


# Contents:

# Part a.  Input the raw data
# 
#install.packages("Hmisc")
#install.packages("moments")
#install.packages("samplesize")
#install.packages("pwr")
#install.packages("tidyverse")
#install.packages("haven")

# Load libraries needed
library(Hmisc)       # for describe() function
library(moments)     # to calculate skewness and kurtosis
library(samplesize)  # for sample size calculations
library(pwr)         # for sample size calculations
library(tidyverse)   # general functions for working with data
library(haven)       # for opening .dta

# Check Working Directory -- .R script and any data must be in the same working directory
#                              Open .R file from folder using RStudio or R sets the working directory
getwd()


sink(file="my_R_log.txt")

#### Part a: Input the raw data
################################################################################

data.stata <- read_dta(file = "T.Martin.dta") 

# Look at structure of data
head(data.stata)
glimpse(data.stata)


# Counting missing values in each column
data.stata %>% summarize_all(~sum(is.na(.)))
# New column for overall change pos
data.stata <- transform(data.stata, POS_dif = POS_ENCOUR_6M - POS_ENCOUR)
# New column for overall change neg
data.stata <- transform(data.stata, NEG_dif = NEG_SABOTAG_6M - NEG_SABOTAG)
# New column for overall change BMI
data.stata <- transform(data.stata, BMI_dif = bmi_6M - bmi)
# New column for overall change BMI
data.stata$type <- ifelse(grepl("-1", data.stata$id, ignore.case = T), "Peer", 
                         ifelse(grepl("-2", data.stata$id, ignore.case = T), "Sidekick", "Other"))


## Convert to long format
longdata <- pivot_longer( data.stata, c(POS_ENCOUR, POS_ENCOUR_3M, POS_ENCOUR_6M, NEG_SABOTAG, NEG_SABOTAG_3M, NEG_SABOTAG_6M))
##New column for positive or negative
longdata$Visit <- ifelse(grepl("POS", longdata$name, ignore.case = T), "Positive", 
                           ifelse(grepl("NEG", longdata$name, ignore.case = T), "Negative", "Other"))
##New column for visit increment
longdata$Visit_inc <- ifelse(grepl("3M", longdata$name, ignore.case = T), "3 Month", 
                         ifelse(grepl("6M", longdata$name, ignore.case = T), "6 Month", "Baseline"))

library(ggplot2)

summary(longdata$Visit == "Positive")
summary(longdata$type == "Peer")
describe(longdata$Visit == "Positive")

### Boxplots showing comparability in score by BMI category
#baseline
data.stata %>% 
  ggplot(aes(factor (BMI_cat), POS_ENCOUR)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Positive Encouragement", title = "Boxplots: Baseline Positive Encouragement")

data.stata %>% 
  ggplot(aes(factor (BMI_cat), POS_ENCOUR)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Positive Encouragement", title = "Boxplots: Baseline Positive Encouragement")


ggsave("fig_0_boxplot_bmi_p.png", width = 6, height = 4, units = "in")
data.stata %>% 
  ggplot(aes(factor (BMI_cat), NEG_SABOTAG)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Negative Sabotage", title = "Boxplots: Baseline Negative Sabotage")

ggsave("fig_0_boxplot_bmi_n.png", width = 6, height = 4, units = "in")
#3month
data.stata %>% 
  ggplot(aes(factor (BMI_cat), POS_ENCOUR_3M)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Positive Encouragement", title = "Boxplots: 3 Month Positive Encouragement")

ggsave("fig_0_boxplot_bmi_p3.png", width = 6, height = 4, units = "in")
data.stata %>% 
  ggplot(aes(factor (BMI_cat), NEG_SABOTAG_3M)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Negative Sabotage", title = "Boxplots: 3 Month Negative Sabotage")

ggsave("fig_0_boxplot_bmi_n3.png", width = 6, height = 4, units = "in")
#6month
data.stata %>% 
  ggplot(aes(factor (BMI_cat), POS_ENCOUR_6M)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Positive Encouragement", title = "Boxplots: 6 Month Positive Encouragement")

ggsave("fig_0_boxplot_bmi_p6.png", width = 6, height = 4, units = "in")
data.stata %>% 
  ggplot(aes(factor (BMI_cat), NEG_SABOTAG_6M)) +
  geom_boxplot() + 
  labs(x = "BMI Category", y = "Negative Sabotage", title = "Boxplots: 6 Month Negative Sabotage")

ggsave("fig_0_boxplot_bmi_n6.png", width = 6, height = 4, units = "in")


## Part a1.  Scatterplot matrix, outcomes by visit

# Using the pairs() function
# Use wide data format to plot visit-specific responses
data.stata %>%
  select(POS_ENCOUR, NEG_SABOTAG, POS_ENCOUR_3M, NEG_SABOTAG_3M, POS_ENCOUR_6M, NEG_SABOTAG_6M) %>%
  pairs(upper.panel=NULL, main="Scatterplot Matrix")

# To save plot
png("fig_a1.png")
data.stata %>%
  select(POS_ENCOUR, NEG_SABOTAG, POS_ENCOUR_3M, NEG_SABOTAG_3M, POS_ENCOUR_6M, NEG_SABOTAG_6M) %>%
  pairs(upper.panel=NULL, main="Scatterplot Matrix")
dev.off()

###########3
# Part a2.  Individual response profiles -- separate time plots for each unit

# Separate graph for each person, separate graphs for the two treatment groups
# Placebo group
longdata %>%
  filter(name == "POS_ENCOUR" |name == "POS_ENCOUR_3M" |name == "POS_ENCOUR_6M",) %>%
  ggplot(aes(x = Visit_inc, y = value)) +
  geom_point(aes()) +
  geom_line(aes(group = Visit, )) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) +
  facet_wrap(~ id) +
  scale_x_discrete(limits = c("Baseline", "3 Month", "6 Month")) +
  labs(x = "Visit number", y = "Positive Encouragement Score", 
       title = "Positive Encouragment Score", subtitle = "Graphs by subject")

ggsave("fig_a2_1.png", width = 6, height = 4, units = "in")
#######
longdata %>%
  filter(name == "NEG_SABOTAG" |name == "NEG_SABOTAG_3M" |name == "NEG_SABOTAG_6M",) %>%
  ggplot(aes(x = Visit_inc, y = value)) +
  geom_point(aes()) +
  geom_line(aes(group = Visit, )) +
  theme(text = element_text(size=8),
        axis.text.x = element_text(angle=90, hjust=1)) +
  facet_wrap(~ id) +
  scale_x_discrete(limits = c("Baseline", "3 Month", "6 Month")) +
  labs(x = "Visit number", y = "Negative Sabotage Score", 
       title = "Negative Sabotage Score", subtitle = "Graphs by subject")


ggsave("fig_a2_2.png", width = 6, height = 4, units = "in")
######
#####3


#complete
# Part a3.  Mean response profiles by treatment

# Create table with mean responses
means_table3 <- longdata %>%
  group_by(Visit_inc, Visit) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), 
            sd_score = sd(value, na.rm = TRUE),
            n = sum(!is.na(value)))


# Print this table
means_table3 



# Plot mean profiles for drug tx and placebo on same graph, connecting x ascending points
means_table3 %>%
  ggplot(aes(x = Visit_inc, y = mean_value)) + 
  geom_point(aes(color = Visit)) + 
  geom_line(aes(group = Visit, color = Visit)) + 
  coord_cartesian(ylim = c(12,32)) +   # give a wider y-axis range than default
  scale_x_discrete(limits = c("Baseline", "3 Month", "6 Month")) +
  labs(x = "Visit increment", y = "Score", 
       title = "Mean score of Positive Encouragement and Negative Sabatoge", color = "Positive/Negative")

ggsave("fig_a3.png", width = 6, height = 4, units = "in")


# Part a5.  Boxplots by visit and treatment

longdata %>%
  ggplot(aes(x = as_factor(Visit_inc), y = value, color = Visit)) + 
  geom_boxplot() +
  scale_x_discrete(limits = c("Baseline", "3 Month", "6 Month")) +
  labs(x = "Visit number", y = "Score", 
       title = "Scores by Pos/Neg and Visit",
       color = "Treatment status")

ggsave("fig_a5.png", width = 6, height = 4, units = "in")


#####
####
longdata_1 <- longdata %>%
  filter(Visit == "Positive")

model_1 <- lm(POS_dif ~ gender + race + age, data = longdata_1)
summary(model_1)

model_2 <- lm(POS_dif ~ gender + race + age + educ, data = longdata_1)
summary(model_2)

model_3 <- lm(POS_dif ~ bmi + bmi_3M + bmi_6M, data = longdata_1)
summary(model_3)

model_4 <- lm(POS_dif ~ bmi + BMI_dif, data = longdata_1)
summary(model_4)

model_5 <- lm(POS_dif ~ BMI_cat + BMI_cat_3M + BMI_cat_6M, data = longdata_1)
summary(model_5)

model_6 <- lm(POS_dif ~ BMI_cat + BMI_dif, data = longdata_1)
summary(model_6)



longdata_2 <- longdata %>%
  filter(Visit == "Negative")

#Negative
model_1.n <- lm(NEG_dif ~ gender + race + age, data = longdata_2)
summary(model_1.n)

model_2.n <- lm(NEG_dif ~ gender + race + age + employ + educ, data = longdata_2)
summary(model_2.n)

model_3.n <- lm(NEG_dif ~ bmi + bmi_3M + bmi_6M, data = longdata_2)
summary(model_3.n)

model_4.n <- lm(NEG_dif ~ bmi + BMI_dif, data = longdata_2)
summary(model_4.n)

model_5 <- lm(NEG_dif ~ BMI_cat + BMI_cat_3M + BMI_cat_6M, data = longdata_2)
summary(model_5)

model_6 <- lm(NEG_dif ~ BMI_cat + BMI_dif, data = longdata_2)
summary(model_6)

#BMI Difference
model_1.b <- lm(BMI_dif ~ gender + race + age, data = longdata)
summary(model_1.b)

model_2.b <- lm(BMI_dif ~ gender + race + age + employ + educ, data = longdata)
summary(model_2.b)

model_3.b <- lm(BMI_dif ~ bmi + bmi_3M + bmi_6M, data = longdata)
summary(model_3.b)

model_4.b <- lm(BMI_dif ~ bmi + BMI_dif, data = longdata)
summary(model_4.b)



sink()
