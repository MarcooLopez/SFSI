#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>
//#include "SFSI.h"
//#include "utils.c"

// ----------------------------------------------------------
// Extract the upper triangular section of a matrix A
//    IN:  n = dimension of the matrix A
//    IN:  diag = whether the diagonal must be included
//    IN:  byrow = whether concatenate elements by row
//    OUT:  a vector of length
//                   n(n+1)/2   if diag = TRUE
//                   n(n-1)/2   if diag = FALSE
// ----------------------------------------------------------
SEXP R_uppertri(SEXP n_, SEXP A_, SEXP diag_, SEXP byrow_)
{
    long long j, tmp;
    int m;
    int inc1=1;

    int n=INTEGER_VALUE(n_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);

    long long n0 = (long long)n;

    tmp = diag ? n0*(n0+1)/2 : n0*(n0-1)/2;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
    double *out = NUMERIC_POINTER(out_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    tmp = 0;
    if(byrow){  // By column: inc=n
      if(diag){ // Include diagonal: start from A[j,j] for 0<=j<n, copy n-j
        for(j=0; j<n; j++){
          m = n-j;
          F77_NAME(dcopy)(&m, A + n0*j+j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[j,j+1] for 0<=j<n-1, copy n-j-1
        for(j=0; j<n-1; j++){
          m = n-j-1;
          F77_NAME(dcopy)(&m, A + n0*(j+1) + j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }

    }else{  // By column: inc=1
      if(diag){ // Include diagonal: start from A[0,j] for 0<=j<n, copy j+1
        for(j=0; j<n; j++){
          m = j+1;
          F77_NAME(dcopy)(&m, A + n0*j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[0,j] for 1<=j<n, copy j
        for(j=1; j<n; j++){
          m = j;
          F77_NAME(dcopy)(&m, A + n0*j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }
    }

    UNPROTECT(2);

    return(out_);
}

// ----------------------------------------------------------
// Extract the lower triangular section of a matrix A
//    IN:  n = dimension of the matrix A
//    IN:  diag = whether the diagonal must be included
//    IN:  byrow = whether concatenate elements by row
//    OUT:  a vector of length
//                   n(n+1)/2   if diag = TRUE
//                   n(n-1)/2   if diag = FALSE
// ----------------------------------------------------------
SEXP R_lowertri(SEXP n_, SEXP A_, SEXP diag_, SEXP byrow_)
{
    long long j, tmp;
    int m;
    int inc1=1;

    int n=INTEGER_VALUE(n_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);

    long long n0 = (long long)n;

    tmp = diag ? n0*(n0+1)/2 : n0*(n0-1)/2;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
    double *out = NUMERIC_POINTER(out_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    tmp = 0;
    if(byrow){  // By column: inc=n
      if(diag){ // Include diagonal: start from A[j,0] for 0<=j<n, copy j+1 entries
        for(j=0; j<n; j++){
          m = j+1;
          F77_NAME(dcopy)(&m, A + j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[0,j] for 1<=j<n, copy j
        for(j=1; j<n; j++){
          m = j;
          F77_NAME(dcopy)(&m, A + j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }

    }else{  // By column: inc=1
      if(diag){ // Include diagonal: start from A[j,j] for 0<=j<n, copy n-j
        for(j=0; j<n; j++){
          m = n-j;
          F77_NAME(dcopy)(&m, A + n0*j + j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[j+1,j] for 0<=j<n-1, copy n-j-1
        for(j=0; j<n-1; j++){
          m = n-j-1;
          F77_NAME(dcopy)(&m, A + n0*j + j + 1, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }
    }

    UNPROTECT(2);

    return(out_);
}

// ----------------------------------------------------------
// ----------------------------------------------------------

SEXP R_unpacktri(SEXP n_, SEXP A_, SEXP uplo_, SEXP diag_, SEXP byrow_)
{
    long long j;
    int m;
    int inc1=1;

    int n=INTEGER_VALUE(n_);
    int uplo=INTEGER_VALUE(uplo_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);

    long long n0 = (long long)n;

    SEXP out_ = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    double *out = NUMERIC_POINTER(out_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    if(((uplo == 1)&&!byrow) || ((uplo == 2)&&byrow))
    {
      if(diag){
        for(j=0; j<n; j++){
          out[n0*j + j] = A[(j+1)*(j+2)/2 -1];
        }

        for(j=1; j<n; j++){   // Start from 2nd column/row to n
          m = j;
          F77_NAME(dcopy)(&m, A + j*(j+1)/2, &inc1, out + n0*j, &inc1);
          F77_NAME(dcopy)(&m, A + j*(j+1)/2, &inc1, out + j, &n);
        }

      }else{
        for(j=0; j<n; j++){
          out[n0*j + j] = 1;
        }

        for(j=1; j<n; j++){   // Start from 2nd column/row to n
          m = j;
          F77_NAME(dcopy)(&m, A + j*(j-1)/2, &inc1, out + n0*j, &inc1);
          F77_NAME(dcopy)(&m, A + j*(j-1)/2, &inc1, out + j, &n);
        }
      }

    }else{
      if(((uplo == 1)&&byrow) || ((uplo == 2)&&!byrow))
      {
        if(diag){
          for(j=0; j<n; j++){
            out[n0*j + j] = A[n0*j - j*(j-1)/2];
          }

          for(j=0; j<n-1; j++){   // Start from 1st column/row to n-1
            m = n-j-1;
            F77_NAME(dcopy)(&m, A + n0*j - j*(j-1)/2 + 1, &inc1, out + n0*(j+1)+j, &n);
            F77_NAME(dcopy)(&m, A + n0*j - j*(j-1)/2 + 1, &inc1, out + n0*j+j+1, &inc1);
          }

        }else{
          for(j=0; j<n; j++){
            out[n0*j + j] = 1;
          }

          for(j=0; j<n-1; j++){   // Start from 1st column/row to n-1
            m = n-j-1;
            F77_NAME(dcopy)(&m, A + n0*j - j*(j+1)/2, &inc1, out + n0*(j+1)+j, &n);
            F77_NAME(dcopy)(&m, A + n0*j - j*(j+1)/2, &inc1, out + n0*j+j+1, &inc1);
          }
        }

      }else{
        Rprintf(" No action performed ...\n");
      }
    }

    UNPROTECT(2);

    return(out_);
}
