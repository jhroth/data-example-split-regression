# Jeremy Roth
rm(list=ls())
library(glmnet)
library(DevTreatRules)
library(DynTxRegime)

# load data
load("Results/datasets_for_data_example_with_missingness_weights.RData")
my.bootstrap.CI.replications <- 1000

## OUTCOME: No breast cancer after 10 years
date.for.test.set.split.regression.BREAST <- "2019-02-19"
load(paste0("Results/", "list_split.regression.no_BREAST_after_10_yr_current_HRT_", date.for.test.set.split.regression.BREAST, ".RData"))
list.build.rules.split.regression.BREAST <- list.build.rules.split.regression
rm(list.build.rules.split.regression, list.evaluate.rules.split.regression)

### evalute locked-in model on test set
set.seed(123)
split.regression.test.BREAST.logistic <- EvaluateRule(data=test.data,
                                                                build.rule.object=list.build.rules.split.regression.BREAST[["propensity_logistic.regression_rule_glm.regression"]],
                                                                study.design="observational",
                                                                additional.weights=test.data$IPW.CC,
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
save(split.regression.test.BREAST.logistic, file=paste0("Results/split_regression_on_test_BREAST_logistic_", Sys.Date(), ".RData"))

set.seed(123)
split.regression.test.BREAST.lasso <- EvaluateRule(data=test.data,
                                                                build.rule.object=list.build.rules.split.regression.BREAST[["propensity_ridge_rule_lasso"]],
                                                                study.design="observational",
                                                                additional.weights=test.data$IPW.CC,
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
save(split.regression.test.BREAST.lasso, file=paste0("Results/split_regression_on_test_BREAST_lasso_", Sys.Date(), ".RData"))

## TREATING NO ONE
B.treat.noone <- rep(0, nrow(test.data))
set.seed(123)
evaluate.treat.noone.test.BREAST <- EvaluateRule(data=test.data,
                                                  build.rule.object=NULL,
                                                  B=B.treat.noone,
                                                  study.design="observational",
                                                  additional.weights=test.data$IPW.CC,
                                                  name.outcome="no_BREAST_after_10_yr",
                                                  type.outcome="binary",
                                                  desirable.outcome=TRUE,
                                                  separate.propensity.estimation=FALSE,
                                                  clinical.threshold=0,
                                                  name.treatment=name.treatment,
                                                  names.influencing.treatment=names.influencing.treatment,
                                                  names.influencing.rule=names.influencing.rule,
                                                  propensity.method="logistic.regression",
                                                  bootstrap.CI=TRUE,
                                                  bootstrap.CI.replications=bootstrap.CI.replications)
save(evaluate.treat.noone.test.BREAST, file=paste0("Results/treat_noone_on_test_BREAST_", Sys.Date(), ".RData"))
