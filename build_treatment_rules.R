# Jeremy Roth
rm(list=ls())
library(DevTreatRules)
library(glmnet)

# load data
load("Data/datasets_for_data_example_with_missingness_weights.RData")
# set parameters
name.outcome <- "no_BREAST_after_10_yr" #"no_CHD_after_10_yr"
bootstrap.CI.replications <- 100

# develop treatment rules on development dataset and compare performance on independent validation dataset
set.seed(123)
model.selection <- CompareRulesOnValidation(development.data=development.data,
                                            validation.data=validation.data,
                                            study.design.development="observational",
                                            vec.approaches="split.regression",
                                            vec.rule.methods=c("ridge", "lasso", "glm.regression"),
                                            vec.propensity.methods=c("ridge", "lasso", "logistic.regression"),
                                            name.outcome.development=name.outcome,
                                            type.outcome.development="binary",
                                            name.treatment.development=name.treatment,
                                            names.influencing.treatment.development=names.influencing.treatment,
                                            names.influencing.rule.development=names.influencing.rule,
                                            desirable.outcome.development=TRUE,
                                            additional.weights.development=development.data$IPW.CC,
                                            additional.weights.validation=validation.data$IPW.CC)
save(model.selection,
     file=paste0("Results/", "split.regression.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))

