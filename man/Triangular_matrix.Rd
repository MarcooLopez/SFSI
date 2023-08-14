
\name{2. Triangular matrix}
\alias{lowerTri}
\alias{upperTri}
\alias{subsetTri}
\alias{unpackTri}
\title{Upper/lower triangular portion from a squared matrix}
\usage{
upperTri(A, diag = TRUE, byrow = FALSE)

lowerTri(A, diag = TRUE, byrow = FALSE)

unpackTri(A, rows = NULL, cols = NULL, drop = TRUE, verbose = TRUE)

subsetTri(A, rowcols, verbose = TRUE)

}
\arguments{
\item{A}{(numeric matrix) A squared matrix of dimension \code{n} or its upper/lower triangular portion}

\item{diag}{\code{TRUE} or \code{FALSE} to whether include the main diagonal when extracting the triangular matrix.
When \code{TRUE} the triangular matrix will include \code{n(n+1)/2}, otherwise, the \code{n(n-1)/2} below/above the diagonal are extracted}

\item{byrow}{\code{TRUE} or \code{FALSE} to whether the triangular portion is obtained row-wise or column-wise}

\item{rows}{(integer vector) Index which rows are to be returned in the full matrix. Default \code{rows=NULL} will return all the \code{n} rows}

\item{cols}{(integer vector) Index which columns are to be returned in the full matrix. Default \code{cols=NULL} return all the \code{n} columns}

\item{drop}{Either \code{TRUE} or \code{FALSE} to whether return a uni-dimensional vector when data is a matrix with either 1 row or 1 column}

\item{rowcols}{(integer vector) Index which rows/columns are to be returned in the triangular matrix}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing matrix information}
}
\value{
Function 'upperTri' and 'lowerTri' return a vector with the entries above/below (and on if \code{diag=TRUE}) the main diagonal. The output will contain attributes 'n', 'diag', 'uplo'. If the input is a symmetric matrix, \code{upperTri(A, byrow=TRUE)} and \code{lowerTri(A, byrow=FALSE)} will produce the same result. The same is true for \code{upperTri(A, byrow=FALSE)} and \code{lowerTri(A, byrow=TRUE)}

Function 'subsetTri' resizes a packed upper/lower triangular matrix to a triangular matrix containing the rows/columns indexed by \code{rowcols}.

Function 'unpackTri' converts a packed upper/lower triangular matrix into a full \code{n}x\code{n} symmetric matrix. If \code{rows} is not NULL, the result will contain \code{length(rows)} rows. The same is applicable for \code{cols}.

These functions are suitable for symmeric matrices only.
}
\description{
Upper/lower triangular portion from a squared matrix
}
\examples{
  require(SFSI)
  data(wheatHTP)
  
  A = cov2cor2(tcrossprod(M[1:100,]))   # A 100x100 symmetric matrix 
  
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
  
  # Recovering the full matrix from the triangular part
  A0 = unpackTri(B1)
  A0[1:5,1:5]
  # it must be equal to A
  max(A-A0)
  
  # Subsetting a triangular matrix
  index = c(20,30,10,5)
  B2 = subsetTri(B1, index)
  # unpacking this must be equal to A[index,index]
  A0 = unpackTri(B2)
  max(A[index,index] - A0)
  
  # If the diagonal is not included, it is assumed to be 1
  B2 = upperTri(A, diag=FALSE)
  length(B1) - length(B2) # it has fewer elements  
  A0 = unpackTri(B1)
  max(A-A0)
  
}
