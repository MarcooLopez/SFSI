\name{BinaryFiles}
\alias{readBinary}
\alias{saveBinary}
\title{Binary files}
\usage{
saveBinary(X, file = paste0(tempdir(), "/file.bin"), 
           precision.format = c("double","single"), 
           verbose = TRUE)
  
readBinary(file = paste0(tempdir(), "/file.bin"), 
           index.row = NULL, index.col = NULL, 
           verbose = TRUE)
}
\arguments{
\item{X}{(numeric matrix) Data to save}

\item{file}{(character) Name of the binary file to save/read}

\item{precision.format}{(character) Either 'single' or 'double' for single (4 bytes) or double precision (8 bytes), respectively,
that matrix to save will occupy}

\item{index.row}{(integer vector) Which rows are to be read from the file. Default \code{index.row=NULL} will read all the rows}

\item{index.col}{(integer vector) Which columns are to be read from the file. Default \code{index.col=NULL} will read all the columns}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing file information}
}
\value{
Function 'saveBinary' does not return any value but print a description of the file saved. 

Function 'readBinary' returns the data that was read.
}
\description{
Save/read a numeric data as a fortran-formatted binary file at a defined precision (single or double).
}
\examples{
  require(SFSI)

  # Simulate matrix
  X = matrix(rnorm(5000),ncol=5)
  head(X)
  
  # Save matrix
  saveBinary(X,paste0(tempdir(),"/Matrix1.bin"),precision.format="double") # as double-precision
  saveBinary(X,paste0(tempdir(),"/Matrix2.bin"),precision.format="single")  # as single-precision

  # Read the double-precision matrix
  X2 = readBinary(paste0(tempdir(),"/Matrix1.bin"))
  head(X2)
  max(abs(X-X2))     # No loss of precision
  object.size(X2)    # Size of the object

  # Read the single-precision matrix
  X3 = readBinary(paste0(tempdir(),"/Matrix2.bin"))
  head(X3)
  max(abs(X-X3))     # Loss of precision
  object.size(X3)    # But smaller-size object

  # Read specific rows and columns
  index.row = c(2,4,5,8,10)
  index.col = c(1,2,5)
  X2 = readBinary(paste0(tempdir(),"/Matrix1.bin"),index.row=index.row,index.col=index.col)
  X2
  # Equal to: 
  X[index.row,index.col]

}
\author{
Marco Lopez-Cruz (\email{maraloc@gmail.com}) and Gustavo de los Campos
}