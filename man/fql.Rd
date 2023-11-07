\name{fql}
\alias{fql}
\title{flexible quasi-likelihood generalized linear model}
\usage{
fql(formula, data, na.action,preci = 0.00001,initial.beta = )
}
\description{
fql is used to fit quasi-likelihood nonparametric generalized linear model. It assummes the variance of response variable y is unknown and use p-spline to estimate the unknown variance function. As a nonparametric model fql is quite flexible on fitting data with complex distribution.
}
\arguments{
  \item{formula}{	an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{ an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{na.action}{	a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is 'na.fail' if that is unset.  Another possible value is NULL, no action.}
  \item{preci}{	a value or cutoff for iterations to stop.}
  \item{initial.beta}{a optional character. By default is "Negative Bionomial" which will perform a negative binomial regression to determine the intial values of coefficients beta. If changed to other characters, the initial values will be set as 0 for all the coefficients beta}
}


\value{
A list of result  with 6 elements.

The first element is the 'estimation', which is a dataframe includes coefficients, standard error and p-value;

The second element is the 'covariance', which is the matrix of variance-covariance of the regression coefficients;

The third element is the 'meanerror', which is the mean of fitted error;

The fourth element is the 'step', which shows the method used to estimate the variance. '0' means that the varicance is estimated by p-spline. '1' for kernel function. '2' uses another kernel method. '3' uses vi=1 for all subjects;

The fifth element is the estimation of the mean;

The sixth element is the estimation of the variance. The details of the variance estimation methods can be seen above in 'step' part.
}

\examples{
library(MASS)
library(fql)
library(mgcv)
data("quine")
# Fit a negative binomial regression model for comparison
nb_model <- glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine)
# Get a summary of the model
summary(nb_model)
fql_model <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#estimation, a dataframe includes coefficients, standard error and p-value
fql_model[[1]]
#covariance, matrix of variance-covariance of the regression coefficients
fql_model[[2]]
#meanerror
fql_model[[3]]
#step', shows the method used to estimate the variance. '0' means that the varicance is estimated by p-spline. '1' for kernel function. '2' uses another kernel method. '3' uses vi=1 for all subjects.
fql_model[[4]]
#estimation of the mean
fql_model[[5]]
#estimation of the variance. The details of the variance estimation methods can be seen above in 'step' part.
fql_model[[6]]
}
