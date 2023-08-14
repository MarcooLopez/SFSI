
\name{Kronecker product}
\alias{Kronecker}
\title{Kronecker product}
\usage{
Kronecker(A, B, rows = NULL, cols = NULL, drop = TRUE, 
          unpack = TRUE, verbose = TRUE) 

}
\arguments{
\item{A}{(numeric matrix) Left numeric matrix}

\item{B}{(numeric matrix) Right numeric matrix}

\item{rows}{(integer vector) Index which rows are to be returned. They must range from 1 to \code{nrow(A)*nrow(B)}. Default \code{rows=NULL} will return all the kroneckered rows}

\item{cols}{(integer vector) Index which columns are to be returned. They must range from 1 to \code{ncol(A)*ncol(B)}. Default \code{cols=NULL} return all the kroneckered columns}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when data is a matrix with either 1 row or 1 column}

\item{unpack}{\code{TRUE} or \code{FALSE} to whether unpack the kronecker to a full matrix when both \code{A} and \code{B} are triangular matrices}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing matrix information}
}
\value{
Returns a matrix containing the full kronecker product. When \code{rows} or \code{cols} are not NULL, the kronecker is formed only from the corresponding elements.
}
\description{
Computes the direct kronecker product between two matrices
}
\details{
For any two matrices
\ifelse{html}{\out{<b>A</b>=\{a<sub>ij</sub>\}}}{\eqn{\textbf{A}=\{a_{ij}\}}{A=\{a_ij\}}} of dimensions nrowA x ncolA and
\ifelse{html}{\out{<b>B</b>=\{b<sub>ij</sub>\}}}{\eqn{\textbf{B}=\{b_{ij}\}}{B=\{b_ij\}}} of dimensions nrowB x ncolB, the direct kronecker product between them is a matrix defined as the block matrix

\ifelse{html}{\out{<p style='text-align:center'>\{a<sub>ij</sub><b>B</b>\}</p>}}{\deqn{\{a_{ij}\textbf{B}\}}{\{a_ij B\}}}

which is of dimensions (nrowA*nrowB) x (ncolA*ncolB)
}
\examples{
  require(SFSI)
  
  # Kronecker product of 2 vectors
  A = rnorm(3)
  B = rnorm(4) 
  K = Kronecker(A, B)
  # it must equal to:
  K0 = kronecker(A, B)
  max(abs(K-K0))
  
  # Kronecker product of 2 matrices
  A = matrix(rnorm(12), ncol=3)
  B = matrix(rnorm(4), ncol=2)
  K = Kronecker(A, B)
  # (it must equal (but faster) to:)
  K0 = kronecker(A, B)
  max(abs(K-K0))
  \donttest{
  # Subsetting rows/columns from the kronecker
  A = matrix(rnorm(150*150), ncol=150)
  B = matrix(rnorm(200*120), ncol=120)
  rows = c(1,3,5,7)
  cols = c(10,20,30,50)
  K = Kronecker(A, B, rows=rows, cols=cols)
  # (it must equal (but much faster) to:)
  K0 = kronecker(A, B)[rows,cols]
  max(abs(K-K0))
  }

}
