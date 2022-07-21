# Biostatistics 624 e-Quiz 4 R Script
# Created by: Leah Jager 2020

# Raw data:        chemo.txt   (text file)

# CSV dataset:   eq4.csv   (created by this program from raw data)


###  Contents:

# Part 0.  Data management for LDA:
#            Convert "wide" format to "long" format
#            Convert missing code from -9 to .
#            Convert baseline rates from 8-week counts to two-week averages
#            Descriptive stats on baseline and follow-up measures by treatment group

# Part a. Exploratory Displays for LDA

# Part a1.  Scatterplot matrix, outcomes by visit
# Part a2.  Individual response profiles -- separate time plots for each unit
# Part a3.  Mean response profiles by treatment
# Part a4.  Boxplots showing comparability in age and baseline seizure counts by treatment group
# Part a5.  Boxplots by visit and treatment
# Part a6.  Boxplots by visit and treatment -- log(count+1) transform
# Part a7.  Boxplots by visit and treatment -- square root(count) transform
# Part a8.  Boxplots by visit and treatment -- log ( (count+1) / BL_count )
# Part a9.  Boxplots for summary measures by treatment


# Part b.   Use MLR on summary measure: average seizures/visit during follow-up

# Part b1.  Generate means within clusters, change from baseline,
#             ratio to baseline, log transforms
# Part b2.  Fit variety of MLRs for summary measures and nonparametric tests
# Part b3   Shapiro-Wilk tests and points with high influence on MLRs for
#             ave. seizures and log ratio
# Part b4   Sensitivity analysis:  log ratio MLR


# Part c.  LDA on repeated measures on log ratio Y = log ( (seiz+1) / (bl_seiz/2) )
#            with adjustments for within-cluster correlations using either:
#              (1) MLR with GEE and robust variance estimation using -- regress
#              (2) Mixed effects / Multilevel model -- xtmixed

# Part c1.  Fit MLR using GEE and robust variance estimation
# Part c2.  Shapiro-Wilk test
# Part c3.  List points with highest influence (dfits)
# Part c5.  Sensitivity analysis:  Re-fit LDA/GEE model excluding subj #49
# Part c6.  Mixed effects / multilevel for repeated measures - xtmixed


# Part d.  Fit log-linear Poisson model using GEE and robust variance estimation,
#            adjust for bl_seiz and age -- excluding subj #49

# Part d1.  Fit log-linear model with poisson command with cluster() robust options
# Part d2.  Get Relative Risks (RR)
# Part d3.  Plot residuals vs visit
# Part d4.  Goodness-of-fit for Poisson model -- check for over/under dispersion
# Part d5.  Does trt effect vary with time or BL levels -- interactions


# Part e.  Fit log-linear model with xtgee, using scaled Pearson's X2 to adjust for
#            Poisson over-dispersion
#           -- alternative to GEE when Poisson model does not fit -- excluding subj #49

# Part e1.  Check for over-dispersed Poisson distribution - visually compare means
#             and variances by treatment group within visits for equality

# Part e2.  Fit GEE Poisson model, excluding Pt # 49
#
# Part e3.  Fit over-dispersed Poisson model, correcting SEs for the
#             the over-dispersion factor (scale from Pearson's X2 per McCullagh & Nelder)

# Part e4.  Show estimate of within group correlation -- intraclass correlation

# Part e5.  Sensitivity analyses -- effects of confounders


# Part f.  Fit log-linear model with BL seizures accounted for as a response
#             Illustrates another way to develop a model that adjusts for baseline measures
#             See analysis of epileptic data in DHLZ text
# Part f1.  Make dataset with records for BL seizures
# Part f2.  Generate offsets.  Models must include an offset, since observation period
#              for BL seizures is 8wks vs. 2wks for follow-up points
# Part f2.  Fit a variety of log-linear models, with/without corrections for over dispersion


# Part g.  Fit random effects models (negative binomial) for count data
#             Illustrates another way to correct for over dispersion /lack of fit
#             when the Poisson model does not fit count data
#
# Part g1. Fit a variety of random effects (negative binomial)
#
#
# Part h.  meqrepoisson -- Multilevel mixed-effects Poisson regression
#             Analyze clustered count data as a multilevel/mixed effects
#             Poisson model



# Load libraries needed
# May need to install these with install.packages() if you don't have them already
library(car)         # for added variable plots, VIF
library(geepack)     # for fitting GEE models
library(lme4)        # for fitting mixed effects models
library(lmerTest)    # for getting p-values for mixed effects models
library(Hmisc)       # describe function, correlation test matrix
library(tidyverse)   # general functions for working with data

# Load the tidyverse last because both the dplyr package from the tidyverse and the Hmisc package have a filter() 
# and summarize() functions and the one loaded last is the one that will be used.  Also both dplyr and MASS have
# select() functions.
# Alternatively, can use dplyr::filter() and dplyr::summarize() to specify the package for the function
# This is what the code below does


#### Part 0: Data management for LDA -- convert "wide" format to "long" format
################################################################################

# Input data from file
# Use read_tsv since data in file are separated by tab characters
# Use meaningful column names
# Specify -9 values be coded as NA
chemo_data_wide <- read_tsv("chemo.dat", 
                            col_names = c("subj", "seiz1", "seiz2", "seiz3", "seiz4", "trt", "bl_seiz8", "age"),
                            na = c("-9"))


# Recode trt as a factor variable (0 = Plbo, 1 = Active)
# This in only non-continuous variable needed to recode
chemo_data_wide <- chemo_data_wide %>%
  mutate(trt = factor(trt, levels=c(0,1), labels=c("Plbo", "Active")))

# Convert baseline rates from 8-week counts to two-week averages
chemo_data_wide <- chemo_data_wide %>%
  mutate(bl_seiz = bl_seiz8/4)


# List first 8 records -- first two patients
head(chemo_data_wide, n = 8)

# Save "Wide" format
write_csv(chemo_data_wide, "eq4wide.csv")

# Convert from "Wide" format into "Long" format
chemo_data_long <- chemo_data_wide %>%
  pivot_longer(cols = starts_with("seiz") ,              # the columns we want to expand
               names_to = "visit" ,                      # name of *new* column to store the visit numbers
               names_prefix = "seiz",                    # what prefix to split off the visit numbers (if not included would get seiz1, seiz2, etc instead of 1, 2, etc)
               ##names_ptypes = list(visit = integer()),   # be sure to convert visit number to integer not character
               values_to = "seiz" )                      # name of *new* column to store values for visits

# Calculate summary measures for each patient

# Average counts during followup
chemo_data_long <- chemo_data_long %>%
  group_by(subj) %>%
  mutate(aveseiz = mean(seiz)) %>%
  ungroup()

# Change from baseline
# Log transform seiz1=seiz+1
# Log relative change from baseline
chemo_data_long <- chemo_data_long %>%
  mutate(diffseiz = aveseiz - bl_seiz,
         logave = log(aveseiz + 1),
         logratio = log((aveseiz + 1)/bl_seiz) )

# List records for first 2 subjects
head(chemo_data_long, n=8)

# Save "Long" format    
write_csv(chemo_data_long, "eq4.csv")

# Rename this long data set to chemo_data to save typing time!
chemo_data <- chemo_data_long

# Look at entire dataset
names(chemo_data)
head(chemo_data)
glimpse(chemo_data)
chemo_data
chemo_data %>% as.data.frame()    # to see whole data set
summary(chemo_data)
describe(chemo_data)

# Generate univarate stats, treated vs placebo

# bl_seiz by trt
chemo_data %>%
  group_by(trt) %>%
  dplyr::summarize(n=n(), 
                   mean = mean(bl_seiz), 
                   sd = sd(bl_seiz), 
                   se = sd/n, 
                   p25 = quantile(bl_seiz, .25), 
                   med = median(bl_seiz),
                   p75 = quantile(bl_seiz, .75),
                   iqr = IQR(bl_seiz),
                   min = min(bl_seiz), 
                   max = max(bl_seiz))

# age by trt
chemo_data %>%
  group_by(trt) %>%
  dplyr::summarize(n=n(), 
                   mean = mean(age), 
                   sd = sd(age), 
                   se = sd/n, 
                   p25 = quantile(age, .25), 
                   med = median(age),
                   p75 = quantile(age, .75),
                   iqr = IQR(age),
                   min = min(age), 
                   max = max(age))

# seiz by trt
chemo_data %>%
  group_by(trt) %>%
  dplyr::summarize(n=n(), 
                   mean = mean(seiz), 
                   sd = sd(seiz), 
                   se = sd/n, 
                   p25 = quantile(seiz, .25), 
                   med = median(seiz),
                   p75 = quantile(seiz, .75),
                   iqr = IQR(seiz),
                   min = min(seiz), 
                   max = max(seiz))


# Boxplots showing comparability in age and baseline seizure counts by treatment group

chemo_data %>% 
  ggplot(aes(x = trt, y = age)) +
  geom_boxplot() + 
  labs(x = "Treatment group", y = "Age (in years)", title = "Boxplots: Baseline Age")

ggsave("fig_0a.png", width = 6, height = 4, units = "in")

chemo_data %>% 
  ggplot(aes(x = trt, y = bl_seiz)) +
  geom_boxplot() + 
  labs(x = "Treatment group", y = "Number of seizures", title = "Boxplots: Baseline seizure rate")

ggsave("fig_0b.png", width = 6, height = 4, units = "in")

# Generate univarate stats for seizure counts, treated vs placebo by visit

# seiz by trt by visit
chemo_data %>%
  group_by(visit, trt) %>%
  dplyr::summarize(n=n(), 
                   mean = mean(seiz), 
                   sd = sd(seiz), 
                   se = sd/n, 
                   p25 = quantile(seiz, .25), 
                   med = median(seiz),
                   p75 = quantile(seiz, .75),
                   iqr = IQR(seiz),
                   min = min(seiz), 
                   max = max(seiz))


# List 3 subjects with highest follow-up seizure counts in each group
chemo_data %>%
  group_by(trt) %>%
  top_n(3, seiz) %>%
  arrange(trt, seiz)


#### Part a: Exploratory Displays for LDA
######################################################################################  

# Part a1.  Scatterplot matrix, outcomes by visit

# Using the pairs() function
# Use wide data format to plot visit-specific responses
chemo_data_wide %>%
  select(trt, bl_seiz, age, seiz1, seiz2, seiz3, seiz4) %>%
  pairs(upper.panel=NULL, main="Scatterplot Matrix")

# To save plot
png("fig_a1.png")
chemo_data_wide %>%
  select(trt, bl_seiz, age, seiz1, seiz2, seiz3, seiz4) %>%
  pairs(upper.panel=NULL, main="Scatterplot Matrix")
dev.off()

# Part a2.  Individual response profiles -- separate time plots for each unit

# Separate graph for each person, separate graphs for the two treatment groups
# Placebo group
chemo_data %>%
  filter(trt == "Plbo") %>%
  ggplot(aes(x = visit, y = seiz)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ subj) +
  labs(x = "Visit number", y = "Number of seizures", 
       title = "Response profiles for Placebo treatment group", subtitle = "Graphs by subject")

ggsave("fig_a2_1.png", width = 6, height = 4, units = "in")

# Active group
chemo_data %>%
  filter(trt == "Active") %>%
  ggplot(aes(x = visit, y = seiz)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ subj) +
  labs(x = "Visit number", y = "Number of seizures", 
       title = "Response profiles for Active treatment group", subtitle = "Graphs by subject")

ggsave("fig_a2_2.png", width = 6, height = 4, units = "in")

# Part a3.  Mean response profiles by treatment

# Create table with mean responses
means_table <- chemo_data %>%
  group_by(visit, trt) %>%
  summarize(mean_seiz = mean(seiz), 
            sd_seiz = sd(seiz),
            n = sum(!is.na(seiz)),
            mean = mean(seiz), 
            sd = sd(seiz), 
            se = sd/n, 
            p25 = quantile(seiz, .25), 
            med = median(seiz),
            p75 = quantile(seiz, .75),
            iqr = IQR(seiz),
            min = min(seiz), 
            max = max(seiz))


# Print this table
means_table


# Plot mean profiles for drug tx and placebo on same graph, connecting x ascending points
means_table %>%
  ggplot(aes(x = visit, y = mean_seiz)) + 
  geom_point(aes(color = trt)) + 
  geom_line(aes(group = trt, color = trt)) + 
  coord_cartesian(ylim = c(0,12)) +   # give a wider y-axis range than default
  labs(x = "Visit number", y = "Number of seizures", 
       title = "Mean response profiles by treatment status", color = "Treatment status")

ggsave("fig_a3.png", width = 6, height = 4, units = "in")


# Part a5.  Boxplots by visit and treatment

chemo_data %>%
  ggplot(aes(x = as_factor(visit), y = seiz, color = trt)) + 
  geom_boxplot() +
  labs(x = "Visit number", y = "Number of seizures", 
       title = "Boxplots: Seizure Rates by Treatment Group and Visit",
       color = "Treatment status")

ggsave("fig_a5.png", width = 6, height = 4, units = "in")

# Part a6.  Boxplots by visit and treatment -- log(count+1) transform

chemo_data <- chemo_data %>%
  mutate(lncnt = log(seiz + 1))

chemo_data %>%
  ggplot(aes(x = as_factor(visit), y = lncnt, color = trt)) + 
  geom_boxplot() +
  labs(x = "Visit number", y = "Log(Seizures + 1)", 
       title = "Boxplots: Log Seizure Rates by Treatment Group and Visit",
       color = "Treatment status")

ggsave("fig_a6.png", width = 6, height = 4, units = "in")

# Part a7.  Boxplots by visit and treatment -- square root(count) transform

chemo_data <- chemo_data %>%
  mutate(sqrtcnt = sqrt(seiz))

chemo_data %>%
  ggplot(aes(x = as_factor(visit), y = sqrtcnt, color = trt)) + 
  geom_boxplot() +
  labs(x = "Visit number", y = "Sqrt(Seizures)", 
       title = "Boxplots: Sqrt Seizure Rates by Treatment Group and Visit",
       color = "Treatment status")

ggsave("fig_a7.png", width = 6, height = 4, units = "in")

# Part a8.  Boxplots by visit and treatment -- log ( (count+1) / BL_count )

chemo_data <- chemo_data %>%
  mutate(logratios = log((seiz + 1)/bl_seiz))

chemo_data %>%
  ggplot(aes(x = as_factor(visit), y = logratios, color = trt)) + 
  geom_boxplot() +
  labs(x = "Visit number", y = "Log ratio", 
       title = "Boxplots: Log Ratios 2-Week Seizure Rates, Follow-up vs. Baseline",
       color = "Treatment status")

ggsave("fig_a8.png", width = 6, height = 4, units = "in")

# Part a9.  Boxplots for summary measures by treatment

chemo_data %>%
  ggplot(aes(x = trt, y = aveseiz)) + 
  geom_boxplot() +
  labs(x = "Treatment group", y = "Number of seizures", 
       title = "Boxplots: Mean Seizure Rate")

ggsave("fig_a9a.png", width = 6, height = 4, units = "in")

chemo_data %>%
  ggplot(aes(x = trt, y = logave)) + 
  geom_boxplot() +
  labs(x = "Treatment group", y = "Log(Number of seizures)", 
       title = "Boxplots: Log Mean Seizure Rate")

ggsave("fig_a9b.png", width = 6, height = 4, units = "in")

chemo_data %>%
  ggplot(aes(x = trt, y = diffseiz)) + 
  geom_boxplot() +
  labs(x = "Treatment group", y = "Change from baseline rate", 
       title = "Boxplots: Change from Baseline Seizure Rate")

ggsave("fig_a9c.png", width = 6, height = 4, units = "in")

chemo_data %>%
  ggplot(aes(x = trt, y = logratio)) + 
  geom_boxplot() +
  labs(x = "Treatment group", y = "Log(Ratio to Baseline)", 
       title = "Boxplots: Log ratio of seizure rate to Baseline rate")

ggsave("fig_a9d.png", width = 6, height = 4, units = "in")



#### Part b: Use MLR on summary response measures:
#
#              Average seizure rate
#              log average seizure rate
#              Change from BL in average seizure rate
#              log ratio follow-up average vs. baseline
######################################################################################  

# Part b2.  Fit variety of MLRs for summary measures and nonparametric tests

# For summary measures -- Must use 1st record only - 1 measure per patient
# To make it easier, let's create a smaller dataset for only visit 1
chemo_data_1 <- chemo_data %>%
  filter(visit == 1)

model_1 <- lm(aveseiz ~ trt + bl_seiz + age, data = chemo_data_1)
summary(model_1)

model_2 <- lm(diffseiz ~ trt + bl_seiz + age, data = chemo_data_1)
summary(model_2)

model_3 <- lm(logave ~ trt + bl_seiz + age, data = chemo_data_1)
summary(model_3)

model_4 <- lm(logratio ~ trt + bl_seiz + age, data = chemo_data_1)
summary(model_4)

# Nonparametric tests:  Wilcoxon ranksum and median tests

wilcox.test(aveseiz ~ trt, data = chemo_data_1)

wilcox.test(logratio ~ trt, data = chemo_data_1)

# Non parametric median test
# Using Fisher's exact test
fisher.test(chemo_data_1$logratio > median(chemo_data_1$logratio), chemo_data_1$trt)


# Part b3   Shapiro-Wilk tests and points with high influence on MLRs for
#             ave. seizures and log ratio

# AVERAGE SEIZURES
model_5 <- lm(aveseiz ~ trt + age, data = chemo_data_1)
summary(model_5)

chemo_data_1 <- chemo_data_1 %>%
  mutate(rstud = rstudent(model_5),
         dfits = dffits(model_5))

shapiro.test(chemo_data_1$rstud)

chemo_data_1 %>% 
  dplyr::select(subj, trt, visit, rstud, dfits) %>%
  filter(dfits <= quantile(dfits, .25) - 1.5*IQR(dfits) | 
           dfits >= quantile(dfits, .75) + 1.5*IQR(dfits)) %>%
  arrange(trt)

# remove rstud and dfits from dataset
chemo_data_1 <- chemo_data_1 %>%
  select(-rstud, -dfits)


# LOG RATIO OF Follow-up to baseline
model_6 <- lm(logratio ~ trt + age, data = chemo_data_1)
summary(model_6)

chemo_data_1 <- chemo_data_1 %>%
  mutate(rstud = rstudent(model_6),
         dfits = dffits(model_6))

shapiro.test(chemo_data_1$rstud)

chemo_data_1 %>% 
  dplyr::select(subj, trt, visit, rstud, dfits) %>%
  filter(dfits <= quantile(dfits, .25) - 1.5*IQR(dfits) | 
           dfits >= quantile(dfits, .75) + 1.5*IQR(dfits)) %>%
  arrange(trt)


# Part b4   Sensitivity analyses:  log ratio MLR
model_7 <- lm(logratio ~ trt + age, data = chemo_data_1, subset = (subj != 58) )
summary(model_7)



#### Part c: LDA on repeated measures on log ratio Y = log ( (seiz+1) / (bl_seiz/2) )
#            with adjustments for within-cluster correlations using either:
#              (1) MLR with GEE and robust variance estimation using -- geeglm() from geepack package
#              (2) Mixed effects / Multilevel model -- lmer() from lme4 package and lmerTest package to get p-values
######################################################################################  

# Part c1.  Fit MLR using GEE and robust variance estimation
# Robust SE estimates are slightly different than ones from Stata
model_gee_1 <- geeglm(logratios ~ trt + bl_seiz + age + visit, id = subj, data = chemo_data)
summary(model_gee_1)

# Get diagnostics
# Re-fit model without cluster adjustments --  diagnostics are approximate
model_gee_2 <- lm(logratios ~ trt + bl_seiz + age + visit, data = chemo_data)

chemo_data <- chemo_data %>%
  mutate(rstud = rstudent(model_gee_2),
         dfits = dffits(model_gee_2))

# Part c2.  Shapiro-Wilk test
shapiro.test(chemo_data$rstud)

# Part c3.  List points with highest influence (dfits)
# List outliers for influence by visit
# slightly different than Stata due to different way of calculating quartiles
chemo_data %>%
  dplyr::select(subj, trt, visit, rstud, dfits) %>%
  filter(dfits <= quantile(dfits, .25) - 1.5*IQR(dfits) | 
           dfits >= quantile(dfits, .75) + 1.5*IQR(dfits)) %>%
  arrange(trt)

# Part c4.  Variance inflation factors
# use vif() function from car package
vif(model_gee_2)


# Part c5.  Sensitivity analysis:  Re-fit LDA/GEE model excluding subj #49
model_gee_3 <- geeglm(logratios ~ trt + bl_seiz + age + visit, id = subj, data = chemo_data, subset = (subj != 49))
summary(model_gee_3)

# Part c6.  Mixed effects / multilevel model for repeated measures -- lmer()

model_lme_1 <- lmer(logratios ~ trt + bl_seiz + age + visit + (1 | subj), data = chemo_data, subset = (subj != 49))
summary(model_lme_1)



#### Part d: Fit log-linear Poisson model using GEE and robust variance estimation,
#            adjust for bl_seiz and age -- comapare to Part c.
######################################################################################  

# Part d1.  Fit log-linear model with poisson command using GEE 
model_pois_1 <- geeglm(seiz ~ trt + bl_seiz + age + visit, id = subj, data = chemo_data, 
                       family = poisson(link = "log"))
summary(model_pois_1)


# Part d2.  Get Relative Risks (RR) and associated confidence intervals by exponentiating estimates and CIs
exp(coef(model_pois_1))
exp(confint.default(model_pois_1))

# Part d3.  Plot residuals vs visit

# Residuals
chemo_data <- chemo_data %>%
  mutate(r_poi = seiz - fitted(model_pois_1))

# Graph residuals vs visits
chemo_data %>%
  ggplot(aes(x = as_factor(visit), y = r_poi)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, color = "red") + 
  labs(x = "Visit number", y = "Residual", title = "Boxplots: Residuals from Poisson Regression by Visit")

ggsave("fig_d4.png", width = 6, height = 4, units = "in")

# Part d4.  Goodness-of-fit for Poisson model -- check for over/under dispersion

# Pearson goodness-of-fit test
pearson.resid <- resid(model_pois_1, type="pearson")    # pearson residuals
x.squared <- sum(pearson.resid^2)                       # test statistic
df <- model_pois_1$df.residual                          # degrees of freedom
pchisq(x.squared, df, lower.tail=FALSE)                 # p-value


# Part d5.  Does trt effect vary with time or BL levels -- interactions

model_pois_2 <- geeglm(seiz ~ trt*visit + bl_seiz + age, id = subj, data = chemo_data,
                       family = poisson(link = "log"))
summary(model_pois_2)

# To test interaction term look at the p-value/CI for the interaction term
# Could also do a LRT comparing to the model without the interaction term
model_pois_0 <- geeglm(seiz ~ trt+visit + bl_seiz + age, id = subj, data = chemo_data,
                       family = poisson(link = "log"))
anova(model_pois_0, model_pois_2, test = "LRT")


model_pois_3 <- geeglm(seiz ~ trt*bl_seiz + age, id = subj, data = chemo_data, 
                       family = poisson(link = "log"))
summary(model_pois_3)

# To test interaction term look at the p-value/CI for the interaction term
# Could also do a LRT comparing to the model without the interaction term
model_pois_0 <- geeglm(seiz ~ trt+bl_seiz + age, id = subj, data = chemo_data,
                       family = poisson(link = "log"))
anova(model_pois_0, model_pois_3, test = "LRT")


#### Part e: Fit log-linear model using scaled Pearson's X2 to adjust for Poisson over-dispersion
#               -- alternative to GEE when Poisson model does not fit -- excluding subj #49
######################################################################################  

# Part e1.  Check for over-dispersed Poisson distribution - visually compare means
#             and variances by treatment group within visits for equality

chemo_data %>%
  group_by(trt, visit) %>%
  summarize(n=n(), mean_seiz = mean(seiz), var_seiz = var(seiz))


# Part e2.  Fit GEE Poisson model, excluding Pt # 49
model_pois_4 <- geeglm(seiz ~ trt + bl_seiz + age + visit, id = subj, data = chemo_data, corstr = "exchangeable", 
                       family = poisson(link = "log"), subset = (subj != 49))
summary(model_pois_4)
summary(model_pois_4)$dispersion   # We see the estimate of the dispersion is 5.03; this is the value we want to correct for in next model


# Part e3.  Fit over-dispersed Poisson model, correcting SEs for the
#             the over-dispersion factor (scale from Pearson's X2 per McCullagh & Nelder)
#           This is done automatically with the geeglm() function through the robust SE estimation
#           So the results above will agree with the corrected version from Stata
summary(model_pois_4)


# Part e4.  Show estimate of within group correlation -- intraclass correlation
summary(model_pois_4)$corr   # this also comes out at the end of summary(model_pois_4)


# Part e5.  Sensitivity analyses -- effects of confounders
model_pois_5 <- geeglm(seiz ~ trt, id = subj, data = chemo_data, corstr = "exchangeable", 
                       family = poisson(link = "log"), subset = (subj != 49))
summary(model_pois_5)


#### Part f: Fit log-linear model with BL seizures accounted for as a response
#             Illustrates another way to develop a model that adjusts for baseline measures
#             See analysis of epileptic data in DHLZ text
######################################################################################  

# Part f1.  Make dataset with records for BL seizures
# Go back to wide data
# Add a seiz0 to represent number of seizures at visit 0 (baseline)
chemo_data_wide <- chemo_data_wide %>%
  mutate(seiz0 = bl_seiz8)

# Part f2.  Generate offsets.  Models must include an offset, since observation period
#    for BL seizures is 8wks vs. 2wks for follow-up points

chemo_data_wide <- chemo_data_wide %>%
  mutate(offset0 = 8,
         offset1 = 2,
         offset2 = 2,
         offset3 = 2,
         offset4 = 2)

# Return to long format
# Now we are making 2 variables (seiz and offset) long from wide
# This is easiest to do if we first put a separator (like a "." or "_") between the stub and the number
chemo_data_wide <- chemo_data_wide %>%
  rename(seiz_0 = seiz0, seiz_1 = seiz1, seiz_2 = seiz2, seiz_3 = seiz3, seiz_4 = seiz4,
         offset_0 = offset0, offset_1 = offset1, offset_2 = offset2, offset_3 = offset3, offset_4 = offset4)

chemo_data_long_2 <- chemo_data_wide %>%
  pivot_longer(cols = -c(subj, trt, bl_seiz8, age, bl_seiz),  # the columns we want to expand, everything minus the ones listed
               names_to = c(".value", "visit") ,               # name of *new* column to store values for visits
               names_sep = "_",                                # the character to separate on
               names_ptypes = list(visit = integer()))         # be sure to convert visit number to integer not character

# Add indicator for follow-up visit
chemo_data_long_2 <- chemo_data_long_2 %>%
  mutate(fup_ind = if_else(visit > 0, 1, 0))



# Part f2. Fit a variety of log-linear models, with/without corrections
#            for over dispersion

# Automatically corrected for overdispersion (compare to the scale(x2) option in Stata)
# Include the offset of log(N) in the formula itself within the offset() function
model_ll_1 <- geeglm(seiz ~ trt*fup_ind + offset(log(offset)), id = subj, data = chemo_data_long_2, corstr = "exchangeable", 
                     family = poisson(link = "log"))
summary(model_ll_1)

model_ll_2 <- geeglm(seiz ~ trt*fup_ind + offset(log(offset)), id = subj, data = chemo_data_long_2, corstr = "exchangeable", 
                     family = poisson(link = "log"), subset = (subj != 49))
summary(model_ll_2)


#### Part g: Fit random effects model (negative binomial) for count data
#             Illustrates another way to correct for over dispersion /lack of fit
#             when the Poisson model does not fit count data
# NOTE: These results do not match the results from Stata
######################################################################################  

# Part g1. Fit a variety of random effects (negative binomial)

model_nb_1 <- glmer.nb(seiz ~ trt*fup_ind + (1 | subj) + offset(log(offset)), data = chemo_data_long_2 )
summary(model_nb_1)

# Exponentiate to get IRR
# Note: The coef() and confint.default() functions don't work for glmer.nb
results <- as_tibble(summary(model_nb_1)$coefficients, rownames = "Variable")
results

results <- results %>%
  mutate(IRR = exp(Estimate),
         CI_IRR_lower = exp(Estimate - 1.96*`Std. Error`),
         CI_IRR_upper = exp(Estimate + 1.96*`Std. Error`))
results

model_nb_2 <- glmer.nb(seiz ~ trt*fup_ind + (1 | subj) + offset(log(offset)), data = chemo_data_long_2, subset = (subj != 49) )
summary(model_nb_2)

# Exponentiate to get IRR
# Note: The coef() and confint.default() functions don't work for glmer.nb
results <- as_tibble(summary(model_nb_2)$coefficients, rownames = "Variable")
results

results <- results %>%
  mutate(IRR = exp(Estimate),
         CI_IRR_lower = exp(Estimate - 1.96*`Std. Error`),
         CI_IRR_upper = exp(Estimate + 1.96*`Std. Error`))
results


#### Part h: glmer -- Multilevel mixed-effects Poisson regression
#             Analyze clustered count data as a multilevel/mixed effects
#             Poisson model
###################################################################################### 

model_pois_6 <- glmer(seiz ~ trt*fup_ind + (1 | subj) + offset(log(offset)), data = chemo_data_long_2, 
                      family = poisson(link = "log") )
summary(model_pois_6)

# Exponentiate to get IRR
# Note: The coef() and confint.default() functions don't work for glmer.nb
results <- as_tibble(summary(model_pois_6)$coefficients, rownames = "Variable")
results <- results %>%
  mutate(IRR = exp(Estimate),
         CI_IRR_lower = exp(Estimate - 1.96*`Std. Error`),
         CI_IRR_upper = exp(Estimate + 1.96*`Std. Error`))
results


model_pois_7 <- glmer(seiz ~ trt*fup_ind + (1 | subj) + offset(log(offset)), data = chemo_data_long_2, subset = (subj != 49),
                      family = poisson(link = "log") )
summary(model_pois_7)

# Exponentiate to get IRR
# Note: The coef() and confint.default() functions don't work for glmer.nb
results <- as_tibble(summary(model_pois_7)$coefficients, rownames = "Variable")
results <- results %>%
  mutate(IRR = exp(Estimate),
         CI_IRR_lower = exp(Estimate - 1.96*`Std. Error`),
         CI_IRR_upper = exp(Estimate + 1.96*`Std. Error`))
results


