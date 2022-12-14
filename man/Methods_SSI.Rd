% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_utils_wC.R
\name{Methods_SSI}
\alias{summary.SSI}
\alias{plot.SSI}
\alias{coef.SSI}
\alias{fitted.SSI}
\title{SSI methods}
\usage{
\method{coef}{SSI}(object, ..., df=NULL, tst=NULL)

\method{fitted}{SSI}(object, ...)

\method{summary}{SSI}(object, ...)
  
\method{plot}{SSI}(..., title=NULL, py=c("accuracy","MSE"))
}
\arguments{
\item{object}{An object of the class 'SSI'. One or more objects must be passed as \code{...} in the method \code{plot}}

\item{df}{(numeric) Average (across testing individuals) number of non-zero regression coefficients}

\item{tst}{(integer vector) Which elements from vector \code{y} (stored in \code{object$y}) are in testing set. Default \code{tst=NULL} will consider all individuals in \code{object$tst}}

\item{py}{(character) Either 'accuracy' or 'MSE' to plot the correlation between observed and predicted values or the mean squared error, respectively, in the y-axis}

\item{title}{(character or expression) Title of the plot}

\item{...}{Arguments to be passed: 
\itemize{
 \item \code{object}: One or more objects of the class 'SSI' (for method \code{plot})
 \item Not needed for methods \code{summary} and \code{fitted}
}
}
}
\value{
Method \code{fitted} returns a matrix with the predicted values for each individual in the testing set (in rows) for each value of lambda (in columns).

Method \code{coef} (list of matrices) returns the regression coefficients for each testing set individual (elements of the list). Each matrix contains the coefficients for each value of lambda (in rows) associated to each training set individual (in columns). If \code{tst} is specified, the elements of the list will correspond only to the testing individuals given in \code{tst}. If \code{df} is specified, only the coefficients for the lambda associated to \code{df} are returned as a 'matrix' with testing individuals in rows. 

Method \code{summary} returns a list object containing:
\itemize{
  \item \code{lambda}: (vector) sequence of (across testing individuals) values of lambda used in the coefficients' estimation.
  \item \code{df}: (vector) degrees of freedom (across testing individuals) at each solution associated to each value of lambda.
  \item \code{accuracy}: (vector) correlation between observed and predicted values associated to each value of lambda.
  \item \code{MSE}: (vector) mean squared error associated to each value of lambda.
  \item \code{optCOR}: (vector) summary of the SSI with maximum accuracy.
  \item \code{optMSE}: (vector) summary of the SSI with minimum MSE.
}
 
Method \code{plot} creates a plot of either accuracy or MSE versus the (average across testing individuals) number of predictors (with non-zero regression coefficient) and versus lambda. 
}
\description{
Useful methods for retrieving, summarizing and visualizing important results from an object of the class 'SSI'
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  X = scale(X[1:200,])/sqrt(ncol(X))    # Subset and scale markers
  G = tcrossprod(X)                     # Genomic relationship matrix
  y = as.vector(scale(Y[1:200,"YLD"]))  # Subset ans scale response variable
  
  fm1 = SSI(y,K=G,tst=1:50,trn=51:length(y))
  
  yHat = fitted(fm1)              # Predicted values for each SSI
  out = summary(fm1)              # Useful function to get results
  corTST = out$accuracy           # Testing set accuracy (correlation cor(y,yHat))
  out$optCOR                      # SSI with maximum accuracy
  out$optMSE                      # SSI with minimum MSE
  B = coef(fm1)                   # Regression coefficients
  B = coef(fm1,df=out$optCOR$df)  # Regression coefficients associated with one 'df'
  plot(fm1,title=expression('corr('*y[obs]*','*y[pred]*') vs sparsity'))   
  plot(fm1,py="MSE",title='Mean Square Error vs sparsity')     
  
}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}
