\name{fql}
\alias{fql}
\title{flexible quasi-likelihood generalized linear model}
\usage{
fql(formula, data, na.action,preci=0.00001)
}
\description{
fql is used to fit quasi-likelihood nonparametric generalized linear model. It assummes the variance of response variable y is unknown and use p-spline to estimate the unknown variance function. As a nonparametric model fql is quite flexible on fitting data with complex distribution.
}
\arguments{
  \item{formula}{	an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{ an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{na.action}{	a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is 'na.fail' if that is unset.  Another possible value is NULL, no action.}
  \item{preci}{	a value or cutoff for iterations to stop.}
}


\value{
#coefficients=bn,covariance=DV1,meanerror=meanerror,step=step
A list of result.
Result[[1]] is the 'estimation', which is a dataframe includes coefficients, standard error and p-value;
Result[[2]] is the 'covariance', which is the matrix of variance-covariance of the regression coefficients;
Result[[3]] is the meanerror;
Result[[4]] is the 'step', which shows which method are used to estimate the variance. '0' means that the varicance is estimated by p-spline. '1' for kernel function. '2' uses another kernel method. '3' uses vi=1 for all subjects.
Result[[5]] is the estimation of the mean.
Result[[6]] is the estimation of the variance. The details of the variance estimation methods can be seen above in 'step' part.
}
