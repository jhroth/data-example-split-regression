# Jeremy Roth
rm(list=ls())
library(DevTreatRules)
library(glmnet)
library(dplyr)

## LOAD WHI-OS data
# load .dat files
data.dir <- "~/Downloads/Documents/WHI_OS_2014b/WHI_OS_2014a/data/whi/ascii/"
read.dat.file <- function(dir, file.name) {
    read.delim(paste0(dir, file.name, ".dat"), sep="\t", header=TRUE)
}
demographics <- read.dat.file(dir=data.dir, "dem_os_pub")
eligibility <- read.dat.file(dir=data.dir, "f2_os_pub")
#personal.information <- read.dat.file(dir=data.dir, "f20_os_pub")
medical.history <- read.dat.file(dir=data.dir, "f30_os_pub")
reproductive.history <- read.dat.file(dir=data.dir, "f31_os_pub")
family.history <- read.dat.file(dir=data.dir, "f32_os_pub")
personal.habits <- read.dat.file(dir=data.dir, "f34_os_pub")
thoughts.and.feelings <- read.dat.file(dir=data.dir, "f37_os_pub")
baseline.questionnaire <- read.dat.file(dir=data.dir, "f42_os_pub")
hormone.use <- read.dat.file(dir=data.dir, "f43_os_pub")
current.medications <- read.dat.file(dir=data.dir, "f44_os_pub")
current.supplements <- read.dat.file(dir=data.dir, "f45_os_pub")
current.supplements <- current.supplements[current.supplements$F45VTYP==1, ] ## restricting to baseline (screening) measurements
outcome.adj <- read.dat.file(dir=data.dir, "outc_adj_os_pub")

apply(demographics, 2, function(x) sum(is.na(x)))
apply(eligibility, 2, function(x) sum(is.na(x)))
apply(medical.history, 2, function(x) sum(is.na(x)))
apply(reproductive.history, 2, function(x) sum(is.na(x)))
apply(family.history, 2, function(x) sum(is.na(x)))
apply(personal.habits, 2, function(x) sum(is.na(x)))
apply(thoughts.and.feelings, 2, function(x) sum(is.na(x)))
apply(hormone.use, 2, function(x) sum(is.na(x)))
apply(current.medications, 2, function(x) sum(is.na(x)))
apply(current.supplements, 2, function(x) sum(is.na(x)))

# create necessary variables
summary(reproductive.history$MENO)
reproductive.history$menopause.before.age.40 <- as.numeric(reproductive.history$MENO < 40)
table(reproductive.history$OOPH, reproductive.history$OOPHA, useNA="always")
reproductive.history$ovaries.removed.before.age.40 <- rep(NA, length(reproductive.history$OOPH))
reproductive.history$ovaries.removed.before.age.40[(reproductive.history$OOPH %in% 0) | (reproductive.history$OOPHA %in% 4:8)] <- 0
reproductive.history$ovaries.removed.before.age.40[reproductive.history$OOPHA %in% 1:3] <- 1
hormone.use$current_HRT <- ifelse(hormone.use$RECENCY == 4, 1, 0)

# formatting CHD variables
outcome.adj$no_CHD <- 1 - outcome.adj[, "CHD"]
outcome.adj$CHD_after_5_yr <- outcome.adj$CHD & outcome.adj$CHDDY < 365*5
outcome.adj$no_CHD_after_5_yr <- 1 - outcome.adj$CHD_after_5_yr
outcome.adj$CHD_after_10_yr <- outcome.adj$CHD & outcome.adj$CHDDY < 365*10
outcome.adj$no_CHD_after_10_yr <- 1 - outcome.adj$CHD_after_10_yr
outcome.adj$CHDDY_replaced <- ifelse(is.na(outcome.adj[, "CHDDY"]), -9999, outcome.adj[, "CHDDY"])

# formatting BREAST variables
outcome.adj$no_BREAST <- 1 - outcome.adj[, "BREAST"]
outcome.adj$BREAST_after_5_yr <- outcome.adj$BREAST & outcome.adj$BREASTDY < 365*5
outcome.adj$no_BREAST_after_5_yr <- 1 - outcome.adj$BREAST_after_5_yr
outcome.adj$BREAST_after_10_yr <- outcome.adj$BREAST & outcome.adj$BREASTDY < 365*10
outcome.adj$no_BREAST_after_10_yr <- 1 - outcome.adj$BREAST_after_10_yr
outcome.adj$BREASTDY_replaced <- ifelse(is.na(outcome.adj[, "BREASTDY"]), -9999, outcome.adj[, "BREASTDY"])

# get names of relevant variables
names.eligibility <- c("CHF_F2", "STROKE", "HYST", "DIAB")
names.factors.eligibility <- c("HEARSTDY")
df.factors.eligibility <- data.frame(apply(as.matrix(eligibility[, names.factors.eligibility]), 2, as.factor), stringsAsFactors=TRUE)
names(df.factors.eligibility) <- names.factors.eligibility
names.demographics <- c("AGE")
names.demographics.missingness <- c("REGION", "BMDFLAG", "LANG") # c("AGE")
names.factors.demographics <- c("ETHNIC", "EDUC", "INCOME")
df.factors.demographics <- data.frame(apply(as.matrix(demographics[, names.factors.demographics]), 2, as.factor), stringsAsFactors=TRUE)
names(df.factors.demographics) <- names.factors.demographics
#names.personal.information <- c("CERVDYS")
names.medical.history <- c("HICHOLRP", "CVD", "CARDCATH", "CABG", "PTCA", "CAROTID",
                                     "AORTICAN", "ANGINA", "PAD", "CANC_F30", "BRCA_F30",
                                     "CERVCA", "OSTEOPOR")
names.factors.medical.history <- c("HTNTRT")
df.factors.medical.history <- data.frame(apply(as.matrix(medical.history[, names.factors.medical.history]), 2, as.factor), stringsAsFactors=TRUE)
names(df.factors.medical.history) <- names.factors.medical.history
names.thoughts.and.feelings <- c("LIFEQUAL")
names.factors.thoughts.and.feelings <- c("GENHEL", "VIGACT", "MODACT", "INTSOC", "ENERGY", "HOTFLASH")
df.factors.thoughts.and.feelings <- data.frame(apply(as.matrix(thoughts.and.feelings[, names.factors.thoughts.and.feelings]), 2, as.factor), stringsAsFactors=TRUE)
names(df.factors.thoughts.and.feelings) <- names.factors.thoughts.and.feelings
names.hormone.use <- c("TOTH", "current_HRT")
names.reproductive.history <- c("MENPSYMP", "PREG", "menopause.before.age.40", "ovaries.removed.before.age.40")
names.factors.reproductive.history <- c("OOPH")
df.factors.reproductive.history <- data.frame(apply(as.matrix(reproductive.history[, names.factors.reproductive.history]), 2, as.factor), stringsAsFactors=TRUE)
names(df.factors.reproductive.history) <- names.factors.reproductive.history
names.current.supplements.missingness <- c("F45MULTI", "F45MVMIN", "F45IRON", "F45STRES")
names.outcome.adj <- c("no_CHD", "CHDDY_replaced", "no_CHD_after_5_yr", "no_CHD_after_10_yr",
                                  "no_BREAST", "BREASTDY_replaced", "no_BREAST_after_5_yr", "no_BREAST_after_10_yr")

# merge by ID into single data frame
df.combined <- left_join(demographics[, c("ID", names.demographics)], data.frame("ID"=demographics[, "ID"], df.factors.demographics)) %>%
                     left_join(., demographics[, c("ID", names.demographics.missingness)], by="ID") %>%
#                     left_join(., personal.information[, c("ID", names.personal.information)], by="ID") %>%
                     left_join(., eligibility[, c("ID", names.eligibility)], by="ID") %>%
                     left_join(., data.frame("ID"=eligibility[, "ID"], df.factors.eligibility), by="ID") %>%
                     left_join(., medical.history[, c("ID", names.medical.history)], by="ID") %>%
                     left_join(., data.frame("ID"=medical.history[, "ID"], df.factors.medical.history), by="ID") %>%
                     left_join(., thoughts.and.feelings[, c("ID", names.thoughts.and.feelings)], by="ID") %>%
                     left_join(., data.frame("ID"=thoughts.and.feelings[, "ID"], df.factors.thoughts.and.feelings), by="ID") %>%
                     left_join(., hormone.use[, c("ID", names.hormone.use)], by="ID") %>% 
                     left_join(., reproductive.history[, c("ID", names.reproductive.history)], by="ID") %>%
                     left_join(., current.supplements[, c("ID", names.current.supplements.missingness)], by="ID") %>% 
                     left_join(., data.frame("ID"=reproductive.history[, "ID"], df.factors.reproductive.history), by="ID") %>%
                     left_join(., outcome.adj[, c("ID", names.outcome.adj)], by="ID")
str(df.combined)

# PREDICTING MISSINGNESS
apply(df.combined, 2, function(x) sum(is.na(x)))
table(complete.cases(df.combined)) ## 16468 of 93676 observations (17.6%) have at least one missing value
df.combined$complete.case.binary <- as.numeric(complete.cases(df.combined))
names.complete.case.predictors <- c("REGION", "BMDFLAG", "LANG", "F45MULTI", "F45MVMIN", "F45IRON", "F45STRES")
logistic.model.complete.case <- glm(NamesToGlmFormula(name.response="complete.case.binary", names.features=names.complete.case.predictors, include.intercept=TRUE),
                              family="binomial", data=df.combined)
coef(logistic.model.complete.case)
predictions.logistic.complete.case <- predict(logistic.model.complete.case, newdata=df.combined, type="response")
summary(predictions.logistic.complete.case) # range is (0.46, 0.90)
df.combined$IPW.CC <- 1 / predictions.logistic.complete.case
summary(df.combined$IPW.CC) # so weights have range of (1.11, 2.15) -- no extreme values
          
# DEFINE ROLES OF VARIABLES
## Y - no CHD after 10 years / no breast cancer after 10 years
## T - current hormone therapy use
## X -  variables from demographics, hormone use, and medication history. 
## L -  Same stuff as X, plus study-specific eligibility and demographics
name.outcome <- "no_CHD_after_10_yr"
name.treatment <-  "current_HRT"
names.influencing.treatment <- c(names.eligibility, names.factors.eligibility,
                                              names.demographics, names.factors.demographics,
                                              names.medical.history, names.factors.medical.history,
                                              names.thoughts.and.feelings, names.factors.thoughts.and.feelings,
                                              names.reproductive.history, names.factors.reproductive.history)
names.influencing.treatment.only <- c("HEARSTDY", "ETHNIC", "EDUC", "INCOME")
names.influencing.rule <- names.influencing.treatment[!(names.influencing.treatment %in% names.influencing.treatment.only)]

# SPLIT THE DATA
df.combined.complete.cases <- df.combined[complete.cases(df.combined[, c(name.outcome, name.treatment, names.influencing.rule, names.influencing.treatment.only)]), ]
set.seed(123)
df.combined.with.split <- split.data(data=df.combined.complete.cases,
                                                n.sets=3, 
                                                split.proportions=c(0.50, 0.25, 0.25))
df.combined.complete.cases <- df.combined[complete.cases(df.combined), ]
training.data <-  df.combined.with.split %>% filter(partition == "training")
validation.data <-  df.combined.with.split %>% filter(partition == "validation")
test.data <-  df.combined.with.split %>% filter(partition == "test")
dim(training.data)
dim(validation.data)
dim(test.data)

## SAVE FORMATTED DATA FOR ANALYSIS
save(df.combined, df.combined.with.split, training.data, validation.data, test.data, 
       name.outcome, name.treatment, names.influencing.treatment, names.influencing.treatment.only, names.influencing.rule, names.complete.case.predictors,
       file="Results/datasets_for_data_example_with_missingness_weights.RData")
