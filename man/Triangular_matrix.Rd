
\name{2. Triangular matrix}
\alias{lowerTri}
\alias{upperTri}
\alias{unpackTri}
\title{Upper/lower triangular portion from a squared matrix}
\usage{
upperTri(A, diag = TRUE, byrow = FALSE)

lowerTri(A, diag = TRUE, byrow = FALSE)

unpackTri(A, n = NULL, uplo = NULL, byrow = NULL)

}
\arguments{
\item{A}{(numeric matrix) A squared matrix of dimension \code{n}}

\item{diag}{\code{TRUE} or \code{FALSE} to whether include the main diagonal when extracting the triangular matrix.
When \code{TRUE} the triangular matrix will include \code{n(n+1)/2}, otherwise, the \code{n(n-1)/2} below/above the diagonal are extracted}

\item{byrow}{\code{TRUE} or \code{FALSE} to whether the triangular portion is obtained row-wise or column-wise}

\item{n}{(integer) Dimension of the squared matrix}

\item{uplo}{(character) Either 'upper' or 'lower'}
}
\value{
Function 'upperTri' and 'lowerTri' return a vector with the entries above/below (and on if \code{diag=TRUE}) the main diagonal. The output will contain attributes 'n', 'diag', 'uplo'. If the input is a symmetric matrix, \code{upperTri(A, byrow=TRUE)} and \code{lowerTri(A, byrow=FALSE)} will produce the same result. The same is true for \code{upperTri(A, byrow=FALSE)} and \code{lowerTri(A, byrow=TRUE)}

Function 'unpackTri' converts a packed upper/lower triangular matrix into a full symmetric matrix of dimension \code{n}.
}
\description{
Upper/lower triangular portion from a squared matrix
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  A = cov2cor2(tcrossprod(M))     # A symmetric matrix 
  
  # Extracting the upper triangular portion
  B1 = upperTri(A)
  str(B1)
  # it must equal (but faster) to:
  B0 = A[upper.tri(A, diag=TRUE)]
  max(abs(B1-B0))
  
  # Extracting row-wise the lower triangular portion
  B2 = lowerTri(A, byrow=TRUE)
  # it must equal to its upper counterpart
  max(abs(B1-B2))
  
  # Recovering the full matrix
  B1 = upperTri(A, diag=TRUE)
  A0 = unpackTri(B1)
  A0[1:10,1:5] 
  
  # If the diagonal is not included, it is assumed to be 1
  B1 = upperTri(A, diag=FALSE)
  A0 = unpackTri(B1)
  A0[1:10,1:5] 
  
}
