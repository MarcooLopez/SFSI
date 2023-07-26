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
// Transform a covariance matrix to a correlation matrix
// The correlation between variables i and j is
// r_ij = v_ij/(sqrt(v_ii)*sqrt(v_jj))
// where sqrt(v_ii) is the SD of the variable i.
//
//    A:       Covariance matrix of dimension n
//    type:    0 = full. For a nxn matrix
//             1 = upper packed triangular
//             2 = lower packed triangular
//             If lower(upper), A must contain the n(n+1)/2
//             entries on or below (above) the diagonal
// ----------------------------------------------------------
SEXP R_cov2cor(SEXP n_, SEXP A_, SEXP a_, SEXP type_, SEXP byrow_)
{
    long long i, j;
    long long pos;

    long long n=(long long)INTEGER_VALUE(n_);
    double a=NUMERIC_VALUE(a_);
    int type=INTEGER_VALUE(type_);
    int byrow=asLogical(byrow_);

    double *sd = (double *) R_alloc(n, sizeof(double)); // Standard deviation (diagonal values)

    int nOK=0;

    PROTECT(A_=AS_NUMERIC(A_));
    double *A=NUMERIC_POINTER(A_);

    if(type == 0){  // Full matrix
      for(j=0; j<n; j++){
        sd[j] = sqrt(A[n*j + j]);
        A[n*j + j] = a*1;
        nOK += isfinite(sd[j]);
      }
      for(j=0; j<n-1; j++){
        for(i=j+1; i<n; i++){
          A[n*j + i] = a*A[n*j + i]/(sd[j]*sd[i]);
          A[n*i + j] = a*A[n*i + j]/(sd[j]*sd[i]);
        }
      }

    }else{
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
        for(j=0; j<n; j++){
          sd[j] = sqrt(A[(j+1)*(j+2)/2 -1]);
          A[(j+1)*(j+2)/2 -1] = a*1;
          nOK += isfinite(sd[j]);
        }

        for(j=1; j<n; j++){   // Start from 2nd colum to n
          pos = j*(j+1)/2; // Start of the column j: A[0,j]
          for(i=0; i<j; i++){
            A[pos + i] = a*A[pos + i]/(sd[i]*sd[j]);
          }
        }

      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
          for(j=0; j<n; j++){
            sd[j] = sqrt(A[n*j - j*(j-1)/2]);
            A[n*j - j*(j-1)/2] = a*1;
            nOK += isfinite(sd[j]);
          }

          for(j=0; j<n-1; j++){   // Start from 1st colum to n-1
            pos = n*j - j*(j-1)/2 + 1;  // Next to diagonal entry
            for(i=0; i<n-j-1; i++){
              A[pos + i] = a*A[pos + i]/(sd[j+i+1]*sd[j]);
            }
          }

        }else{
          Rprintf(" No action performed ...\n");
        }
      }
    }

    UNPROTECT(1);

    return(Rf_ScalarInteger(nOK));
}
