#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>
#include "SFSI.h"
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
SEXP R_upper_tri(SEXP n_, SEXP A_, SEXP diag_, SEXP byrow_)
{
    long long j, tmp;
    int m;
    int inc1=1;

    int n=INTEGER_VALUE(n_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    tmp = diag ? (long long)n*((long long)n+1)/2 : (long long)n*((long long)n-1)/2;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
    double *out = NUMERIC_POINTER(out_);

    tmp = 0;
    if(byrow){  // By column: inc=n
      if(diag){ // Include diagonal: start from A[j,j] for 0<=j<n, copy n-j
        for(j=0; j<n; j++){
          m = n-j;
          F77_NAME(dcopy)(&m, A + n*j+j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[j,j+1] for 0<=j<n-1, copy n-j-1
        for(j=0; j<n-1; j++){
          m = n-j-1;
          F77_NAME(dcopy)(&m, A + n*(j+1) + j, &n, out + tmp, &inc1);
          tmp += m;
        }
      }

    }else{  // By column: inc=1
      if(diag){ // Include diagonal: start from A[0,j] for 0<=j<n, copy j+1
        for(j=0; j<n; j++){
          m = j+1;
          F77_NAME(dcopy)(&m, A + n*j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[0,j] for 1<=j<n, copy j
        for(j=1; j<n; j++){
          m = j;
          F77_NAME(dcopy)(&m, A + n*j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }
    }

    // Set attributes
    setAttrib(out_, mkString("uplo"), mkString("upper"));
    setAttrib(out_, mkString("n"), n_);
    setAttrib(out_, mkString("include.diag"), diag_);
    setAttrib(out_, mkString("byrow"), byrow_);

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
SEXP R_lower_tri(SEXP n_, SEXP A_, SEXP diag_, SEXP byrow_)
{
    long long j, tmp;
    int m;
    int inc1=1;

    int n=INTEGER_VALUE(n_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    tmp = diag ? (long long)n*((long long)n+1)/2 : (long long)n*((long long)n-1)/2;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
    double *out = NUMERIC_POINTER(out_);

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
          F77_NAME(dcopy)(&m, A + n*j + j, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }else{   // Does not include diagonal: start from A[j+1,j] for 0<=j<n-1, copy n-j-1
        for(j=0; j<n-1; j++){
          m = n-j-1;
          F77_NAME(dcopy)(&m, A + n*j + j + 1, &inc1, out + tmp, &inc1);
          tmp += m;
        }
      }
    }

    // Set attributes
    setAttrib(out_, mkString("uplo"), mkString("lower"));
    setAttrib(out_, mkString("n"), n_);
    setAttrib(out_, mkString("include.diag"), diag_);
    setAttrib(out_, mkString("byrow"), byrow_);

    UNPROTECT(2);

    return(out_);
}

// ----------------------------------------------------------
// ----------------------------------------------------------

SEXP R_subset_tri(SEXP nA_, SEXP A_, SEXP index_,
                  SEXP uplo_, SEXP diag_, SEXP byrow_)
{
    int nA=INTEGER_VALUE(nA_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);
    long long ni=Rf_isNull(index_) ? 0 : XLENGTH(index_);

    int uplo=strcmp(CHAR(asChar(uplo_)),"upper")==0 ? 1 : 2;

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(index_=AS_INTEGER(index_));
    int *index=INTEGER_POINTER(index_);

    long long tmp = diag ? ni*(ni+1)/2 : ni*(ni-1)/2;
    SEXP out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
    double *out = NUMERIC_POINTER(out_);

    if(((uplo == 1)&&!byrow) || ((uplo == 2)&&byrow))
    {
      subset_tri1(ni, nA, A, diag, ni, index, out);

    }else{
      if(((uplo == 1)&&byrow) || ((uplo == 2)&&!byrow))
      {
        subset_tri2(ni, nA, A, diag, ni, index, out);

      }else{
        Rprintf(" No action performed ...\n");
      }
    }

    // Set attributes
    setAttrib(out_, mkString("uplo"), uplo_);
    setAttrib(out_, mkString("n"), Rf_ScalarInteger(ni));
    setAttrib(out_, mkString("include.diag"), diag_);
    setAttrib(out_, mkString("byrow"), byrow_);

    UNPROTECT(3);

    return(out_);
}

// ----------------------------------------------------------
// ----------------------------------------------------------

SEXP R_unpack_tri(SEXP nA_, SEXP A_,
                  SEXP irow_, SEXP icol_, SEXP drop_,
                  SEXP uplo_, SEXP diag_, SEXP byrow_)
{

    int nA=INTEGER_VALUE(nA_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);
    int drop=asLogical(drop_);
    long long nirow=Rf_isNull(irow_) ? 0 : XLENGTH(irow_);
    long long nicol=Rf_isNull(icol_) ? 0 : XLENGTH(icol_);

    int uplo=strcmp(CHAR(asChar(uplo_)),"upper")==0 ? 1 : 2;

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(irow_=AS_INTEGER(irow_));
    int *irow=INTEGER_POINTER(irow_);

    PROTECT(icol_=AS_INTEGER(icol_));
    int *icol=INTEGER_POINTER(icol_);

    int nrow = nirow == 0 ? nA : nirow;
    int ncol = nicol == 0 ? nA : nicol;

    SEXP out_;
    if((nrow==1) || (ncol==1))
    {
      if(drop){
        out_ = PROTECT(Rf_allocVector(REALSXP, (long long)nrow*(long long)ncol));
      }else{
        out_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
      }
    }else{
      out_ = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
    }
    double *out = NUMERIC_POINTER(out_);

    int case_set = 1*(nirow>0) + 2*(nicol>0);

    if(((uplo == 1)&&!byrow) || ((uplo == 2)&&byrow))
    {
       unpack_tri1(case_set, nrow, ncol, nA, A, diag, irow, icol, out);

    }else{
      if(((uplo == 1)&&byrow) || ((uplo == 2)&&!byrow))
      {
        unpack_tri2(case_set, nrow, ncol, nA, A, diag, irow, icol, out);

      }else{
        Rprintf(" No action performed ...\n");
      }
    }

    UNPROTECT(4);

    return(out_);
}
