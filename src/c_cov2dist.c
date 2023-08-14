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
// Transform a covariance matrix to a distance matrix
// The distance between variables i and j is
// d_ij = a*sqrt(v_ii + v_jj -2*v_ij)
// where v_ii is the variance (diagonal value) of the variable i.
//
//    A:       Covariance matrix of dimension n
//    type:    0 = full. For a nxn matrix
//             1 = upper packed triangular
//             2 = lower packed triangular
//  If lower(upper), A must contain the n(n+1)/2 entries on or
//  below (above) the diagonal
// ----------------------------------------------------------
SEXP R_cov2dist(SEXP n_, SEXP a_, SEXP A_, SEXP type_, SEXP byrow_)
{
    long long i, j;
    long long pos;

    int n=INTEGER_VALUE(n_);
    double a=NUMERIC_VALUE(a_);
    int type=INTEGER_VALUE(type_);
    int byrow=asLogical(byrow_);

    double *v = (double *) R_alloc(n, sizeof(double)); // Variances: Diagonal values

    PROTECT(A_=AS_NUMERIC(A_));
    double *A=NUMERIC_POINTER(A_);

    if(type == 0){  // Full matrix
      for(j=0; j<n; j++){   // Diagonal values
        v[j] = A[n*j + j];
        A[n*j + j] = 0;
      }

      for(j=0; j<n-1; j++){
        for(i=j+1; i<n; i++){
          A[n*j + i] = a*sqrt(v[i]+v[j] -2*A[n*j + i]);
          A[n*i + j] = a*sqrt(v[i]+v[j] -2*A[n*i + j]);
        }
      }

    }else{
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
        for(j=0; j<n; j++){
          v[j] = A[(j+1)*(j+2)/2 -1];
          A[(j+1)*(j+2)/2 -1] = 0;
        }

        for(j=1; j<n; j++){   // Start from 2nd row/colum to n
          pos = j*(j+1)/2; // Start of the row/column j
          for(i=0; i<j; i++){
            A[pos + i] = a*sqrt(v[i]+v[j] -2*A[pos + i]);
          }
        }

      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
          for(j=0; j<n; j++){
            v[j] = A[n*j - j*(j-1)/2];
            A[n*j - j*(j-1)/2] = 0;
          }

          for(j=0; j<n-1; j++){   // Start from 1st colum to n-1
            pos = n*j - j*(j-1)/2 + 1;  // Next to diagonal entry
            for(i=0; i<n-j-1; i++){
              A[pos + i] = a*sqrt(v[j+i+1]+v[j] -2*A[pos + i]);
            }
          }

        }else{
          Rprintf(" No action performed ...\n");
        }
      }
    }

    UNPROTECT(1);

    return(R_NilValue);
}
