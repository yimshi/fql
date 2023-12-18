# fql
# install the package
library(devtools) 

install_github('yimshi/fql')

# run an example
library(MASS)

library(fql)

library(mgcv)

data("quine")
# Fit a negative binomial regression model
nb_model <- glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine)
# Get a summary of the model
summary(nb_model)

fql_model <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
# 'estimation', which is a dataframe includes coefficients, standard error and p-value
fql_model[[1]]
# 'covariance', which is the matrix of variance-covariance of the regression coefficients
fql_model[[2]]
# 'meanerror', which is the mean of fitted error
fql_model[[3]]
# 'step', which shows the method used to estimate the variance. '0' means that the varicance is estimated by p-spline. '1' for kernel function. '2' uses another kernel method. '3' uses vi=1 for all subjects
fql_model[[4]]
# estimation of the mean
fql_model[[5]]
# the estimation of the variance.
fql_model[[6]]
