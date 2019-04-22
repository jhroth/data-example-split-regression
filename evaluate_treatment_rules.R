# Jeremy Roth
rm(list=ls())
library(glmnet)
library(DevTreatRules)

# load data
load("Results/datasets_for_data_example_with_missingness_weights.RData")
my.bootstrap.CI.replications <- 10

## OUTCOME: No breast cancer after 10 years
date.for.evaluation.set.split.regression.BREAST <- "2019-04-22"
load(paste0("Results/", "split.regression.no_BREAST_after_10_yr_current_HRT_", date.for.evaluation.set.split.regression.BREAST, ".RData"))
list.build.rules.split.regression.BREAST <- model.selection$list.rules$split.regression
rm(model.selection)

### evalute locked-in model on evaluation set
set.seed(123)
split.regression.evaluation.BREAST.logistic <- EvaluateRule(data=evaluation.data,
                                                                BuildRule.object=list.build.rules.split.regression.BREAST[["propensity_logistic.regression_rule_glm.regression"]],
                                                                study.design="observational",
                                                                additional.weights=evaluation.data$IPW.CC,
                                                                name.outcome="no_BREAST_after_10_yr",
                                                                type.outcome="binary",
                                                                desirable.outcome=TRUE,
                                                                separate.propensity.estimation=TRUE,
                                                                clinical.threshold=0,
                                                                name.treatment=name.treatment,
                                                                names.influencing.treatment=names.influencing.treatment,
                                                                names.influencing.rule=names.influencing.rule,
                                                                propensity.method="logistic.regression",
                                                                bootstrap.CI=TRUE,
                                                                bootstrap.CI.replications=my.bootstrap.CI.replications)
save(split.regression.evaluation.BREAST.logistic, file=paste0("Results/split_regression_on_evaluation_BREAST_logistic_", Sys.Date(), ".RData"))

# develop treatment rules on development dataset and compare performance on independent validation dataset
set.seed(123)
split.regression.evaluation.BREAST.lasso <- EvaluateRule(data=evaluation.data,
                                                                build.rule.object=list.build.rules.split.regression.BREAST[["propensity_ridge_rule_lasso"]],
                                                                study.design="observational",
                                                                additional.weights=evaluation.data$IPW.CC,
                                                                name.outcome="no_BREAST_after_10_yr",
                                                                type.outcome="binary",
                                                                desirable.outcome=TRUE,
                                                                separate.propensity.estimation=TRUE,
                                                                clinical.threshold=0,
                                                                name.treatment=name.treatment,
                                                                names.influencing.treatment=names.influencing.treatment,
                                                                names.influencing.rule=names.influencing.rule,
                                                                propensity.method="logistic.regression",
                                                                bootstrap.CI=TRUE,
                                                                bootstrap.CI.replications=my.bootstrap.CI.replications)
save(split.regression.evaluation.BREAST.lasso, file=paste0("Results/split_regression_on_evaluation_BREAST_lasso_", Sys.Date(), ".RData"))

