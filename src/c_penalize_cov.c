#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>
#include "SFSI.h"
//#include "utils.c"


// ----------------------------------------------------------
// Penalize a numeric covariance matrix V by adding to the
// diagonal, a penalty lambda
//    lambda:  Numeric value to be added
//        A:   covariance matrix
//    type:    0 = full. For a nxn matrix
//             1 = upper packed triangular
//             2 = lower packed triangular
//             If lower(upper), A must contain the n(n+1)/2
//             entries on or below (above) the diagonal
// ----------------------------------------------------------
SEXP R_penalize_cov(SEXP n_, SEXP A_, SEXP lambda_, SEXP type_, SEXP byrow_)
{
    long long j;

    long long n=(long long)INTEGER_VALUE(n_);
    int nlambda=XLENGTH(lambda_);
    int type=INTEGER_VALUE(type_);
    int byrow=asLogical(byrow_);

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(lambda_=AS_NUMERIC(lambda_));
    double *lambda=NUMERIC_POINTER(lambda_);

    if(type == 0){  // Full matrix
      if(nlambda == n){
        for(j=0; j<n; j++){
          A[n*j + j] += lambda[j];
        }
      }else{
        for(j=0; j<n; j++){
          A[n*j + j] += lambda[0];
        }
      }

    }else{
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
        if(nlambda == n){
          for(j=0; j<n; j++){
            A[(j+1)*(j+2)/2 -1] += lambda[j];
          }
        }else{
          for(j=0; j<n; j++){
            A[(j+1)*(j+2)/2 -1] += lambda[0];
          }
        }

      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
          if(nlambda == n){
            for(j=0; j<n; j++){
              A[n*j - j*(j-1)/2] += lambda[j];
            }
          }else{
            for(j=0; j<n; j++){
              A[n*j - j*(j-1)/2] += lambda[0];
            }
          }

        }else{
          Rprintf(" No action performed ...\n");
        }
      }
    }

    UNPROTECT(2);

    return(R_NilValue);
}
