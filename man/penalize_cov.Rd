
\name{Ridge penalization}
\alias{penalize_cov}
\title{Penalization of a covariance matrix}
\usage{
penalize_cov(A, a = 1, lambda, inplace = FALSE, 
             rows = NULL, cols = NULL, drop = TRUE, 
             unpack = TRUE, verbose = TRUE)

}
\arguments{
\item{A}{(numeric matrix) Variance-covariance matrix}

\item{inplace}{\code{TRUE} or \code{FALSE} to whether operate directly on the input matrix.
When \code{TRUE} no result is produced but the input \code{A} is modified. Default \code{inplace=FALSE}}

\item{a}{(numeric) A number to multiply the whole matrix \code{A} before penalization. Default \code{a=1}}

\item{lambda}{(numeric) Value to be added to the diagonal of \code{aA}}

\item{rows}{(integer vector) Index which rows are to be returned. Default \code{rows=NULL} will return all the rows}

\item{cols}{(integer vector) Index which columns are to be returned. Default \code{cols=NULL} return all the columns}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when data is a matrix with either 1 row or 1 column}

\item{unpack}{\code{TRUE} or \code{FALSE} to whether unpack the output to a full matrix when \code{A} is a triangular matrix and \code{a} is an scalar}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing matrix information}
}
\value{
Returns the penalized matrix
}
\description{
Ridge penalization of a covariance matrix
}
\details{
In Ridge Regression, the problem of a near-singular variance-covariance matrix 
(\ifelse{html}{\out{<b>A</b>}}{\eqn{\textbf{A}}{A}}) of dimension \code{n}
is alleviated by adding positive elements to the diagonal entries, thereby decreasing its condition number. If the matrix is multiplied first by a factor \code{a}, this is

\ifelse{html}{\out{<p style='text-align:center'>a<b>A</b> + <b>D</b>(&lambda;)</p>}}{\deqn{a\textbf{A}+\textbf{D}(\lambda)}{aA + D(lambda)}}

where \ifelse{html}{\out{<b>D</b>(&lambda;)}}{\eqn{\textbf{D}(\lambda)}{D(lambda)}}
is a diagonal matrix of dimension \code{n} containing the ridge parameter
\ifelse{html}{\out{&lambda;≥0}}{\eqn{\lambda\geq 0}{lambda >= 0}} working as a diagonal-shifting parameter. This ridge parameter can be generalized to contain diagonal specific values
\ifelse{html}{\out{<b>D</b>(&lambda;<sub>1</sub>,...,&lambda;<sub>n</sub>)}}{\eqn{\textbf{D}(\lambda_{1},...,\lambda_{n})}{D(lambda_1,...,lambda_n)}}.

For the multitrait case, \code{a} and \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} might take a matrix form and in this case, the penalization will adopt a kronecker (see \code{help(kronecker)}) form:

\ifelse{html}{\out{<p style='text-align:center'>kronecker(<b>a</b>,<b>A</b>) + kronecker(<b>&lambda;</b>,<b>D</b>)</p>}}{\deqn{kronecker(\textbf{a},\textbf{A}) + kronecker(\boldsymbol{\lambda},\textbf{D})}{kronecker(a,A) + kronecker(lambda,D)}}

}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  A = cov2cor2(tcrossprod(M[1:100,]))   # A 100x100 covariance matrix 
  A[1:5,1:5]
  
  # Penalization when 'a' and 'lambda' are scalars
  # Multiply the matrix A by 2 then sum to the diagonal 5
  K = penalize_cov(A, a=2, lambda=5)
  K[1:5,1:5]
  
  # Using specific shifting values
  lambda = rnorm(nrow(A))^2
  K = penalize_cov(A, a=2, lambda=lambda)
  K[1:5,1:5]
  
  # Kronecker product
  # Penalization when 'a' and 'lambda' are matrices
  a = cor(Y[,4:7])         # 4x4 matrix
  lambda = diag(rnorm(4))  # 4x4 matrix
  K = penalize_cov(A, a=a, lambda=lambda)
  # (it must equal (but faster) to:)
  K0 = kronecker(a,A) + kronecker(lambda,diag(nrow(A)))
  max(K-K0)
  \donttest{
  # Subsetting rows/columns from the kronecker
  n = 100
  a = crossprod(matrix(rnorm(20*n), ncol=n))      # 100x100 matrix
  lambda = crossprod(matrix(rnorm(20*n), ncol=n)) # 100x100 matrix
  rows = c(1,3,5,7)
  cols = c(10,20,30,50)
  K = penalize_cov(A, a=a, lambda=lambda, rows=rows, cols=cols)
  # (it must equal (but much faster) to:)
  K0 = (kronecker(a,A) + kronecker(lambda,diag(nrow(A))))[rows,cols]
  max(abs(K-K0))
  }
  
}
