#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>
#include "SFSI.h"
//#include "utils.c"

//====================================================================
// Penalize a numeric covariance matrix V by adding to the
// diagonal, a penalty lambda
//    lambda:  Numeric value to be added
//        A:   covariance matrix
//    type:    0 = full. For a nxn matrix
//             1 = upper packed triangular
//             2 = lower packed triangular
//             If lower(upper), A must contain the n(n+1)/2
//             entries on or below (above) the diagonal
//====================================================================
SEXP R_penalize_cov(SEXP a_, SEXP nA_, SEXP A_, SEXP lambda_,
                    SEXP irow_, SEXP icol_, SEXP drop_,
                    SEXP type_, SEXP diag_, SEXP byrow_,
                    SEXP scale_, SEXP inplace_,
                    SEXP unpack_, SEXP verbose_)
{
    long long i, j, tmp;
    double *out;
    int itmp;
    double *Aj;
    int inc1=1;
    SEXP out_=NULL;
    int nprotect=4;

    int nA=INTEGER_VALUE(nA_);
    double a=NUMERIC_VALUE(a_);
    int nlambda=isMatrix(lambda_) ? nrows(lambda_) : XLENGTH(lambda_);
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);
    int scale=asLogical(scale_);
    int inplace=asLogical(inplace_);
    int drop=asLogical(drop_);
    int unpack=asLogical(unpack_);
    int verbose=asLogical(verbose_);
    long long nirow=Rf_isNull(irow_) ? 0 : XLENGTH(irow_);
    long long nicol=Rf_isNull(icol_) ? 0 : XLENGTH(icol_);

    int type=strcmp(CHAR(asChar(type_)),"full")==0 ? 0 : strcmp(CHAR(asChar(type_)),"upper")==0 ? 1 : 2;

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(lambda_=AS_NUMERIC(lambda_));
    double *lambda=NUMERIC_POINTER(lambda_);

    PROTECT(irow_=AS_INTEGER(irow_));
    int *irow=INTEGER_POINTER(irow_);

    PROTECT(icol_=AS_INTEGER(icol_));
    int *icol=INTEGER_POINTER(icol_);

    int nrow = nirow == 0 ? nA : nirow;
    int ncol = nicol == 0 ? nA : nicol;
    int case_set = 1*(nirow>0) + 2*(nicol>0);
    int totri = (type>0) && all_equal(nirow,irow,nicol,icol);

    if(!diag){
      if(totri){
        Rprintf(" Routine can not be implemented for a triangular matrix with no diagonal\n");
        UNPROTECT(nprotect);
        return(R_NilValue);
      }else{
        Rprintf(" Matrix does not contain diagonal data, a diagonal equal to 1 is assumed\n");
      }
    }

    if(totri){
      if(unpack){
        totri = 0;
        if(verbose){
          Rprintf(" Matrix was unpacked to a full matrix\n");
        }
      }
    }else{
      if(verbose & (type>0)){
        Rprintf(" Matrix was unpacked to a full matrix when selecting different rows/cols\n");
      }
    }

    if(inplace){
      if(case_set > 0){
        Rprintf(" Routine can be implemented 'inplace' only when both 'rows' and 'cols' are NULL\n");
        UNPROTECT(nprotect);
        return(R_NilValue);
      }
      out = A;
    }else{
      if(totri){
        tmp = (long long)ncol*((long long)ncol+1)/2;
        out_ = PROTECT(Rf_allocVector(REALSXP, tmp));
      }else{
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
      }
      out=NUMERIC_POINTER(out_);
      nprotect++;
    }

    if(type == 0){  // Full matrix
      if(!inplace){
        //Rprintf(" Making copy of the matrix ...\n");
        switch(case_set){
          case 0: // nirow=0 && nicol=0
            for(j=0; j<ncol; j++){
              F77_NAME(dcopy)(&nrow, A + nA*j, &inc1, out + nrow*j, &inc1);
            }
          break;

          case 1: // nirow>0 && nicol=0
            for(j=0; j<ncol; j++){
              Aj = A + nA*j;
              for(i=0; i<nrow; i++){
                out[nrow*j + i] = Aj[irow[i]];
              }
            }
          break;

          case 2: // nirow=0 && nicol>0
            for(j=0; j<ncol; j++){
               F77_NAME(dcopy)(&nrow, A + nrow*(long long)icol[j], &inc1, out + nrow*j, &inc1);
            }
          break;

          case 3: // nirow>0 && nicol>0
            for(j=0; j<ncol; j++){
              Aj = A + nA*(long long)icol[j];
              for(i=0; i<nrow; i++){
                out[nrow*j + i] = Aj[irow[i]];
              }
            }
          break;
        }
      }

    }else{
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
        if(!inplace){
          //Rprintf(" Making copy of the matrix ...\n");
          if(totri){
            subset_tri1(ncol, nA, A, 1, nicol, icol, out);
          }else{
            unpack_tri1(case_set, nrow, ncol, nA, A, diag, irow, icol, out);
          }
        }

      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
          if(!inplace){
            //Rprintf(" Making copy of the matrix ...\n");
            if(totri){
              subset_tri2(ncol, nA, A, 1, nicol, icol, out);
            }else{
              unpack_tri2(case_set, nrow, ncol, nA, A, diag, irow, icol, out);
            }
          }

        }else{
          Rprintf(" No action performed ...\n");
        }
      }
    }

    if(totri){  // Triangular matrix that remain triangular
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
        if(scale){
          //Rprintf(" Making pointwise multiplication by scalar 'a'...\n");
          for(j=0; j<ncol; j++){
            itmp = j+1;
            F77_NAME(dscal)(&itmp, &a, out + j*(j+1)/2, &inc1);
          }
        }

        //Rprintf(" Making shifting ...\n");
        if(nicol == 0){
          for(j=0; j<ncol; j++){
            out[(j+1)*(j+2)/2 -1] += lambda[nlambda == 1 ? 0 : j];
          }
        }else{
          for(j=0; j<ncol; j++){
            out[(j+1)*(j+2)/2 -1] += lambda[nlambda == 1 ? 0 : icol[j]];
          }
        }

      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
          if(scale){
            //Rprintf(" Making pointwise multiplication by scalar 'a'...\n");
            for(j=0; j<ncol; j++){
              itmp = ncol-j;
              F77_NAME(dscal)(&itmp, &a, out + ncol*j - j*(j-1)/2, &inc1);
            }
          }

          //Rprintf(" Making shifting ...\n");
          if(nicol == 0){
            for(j=0; j<ncol; j++){
              out[ncol*j - j*(j-1)/2] += lambda[nlambda == 1 ? 0 : j];
            }
          }else{
            for(j=0; j<nA; j++){
              out[ncol*j - j*(j-1)/2] += lambda[nlambda == 1 ? 0 : icol[j]];
            }
          }

        }else{
          Rprintf(" No action performed ...\n");
        }
      }

      // Set attributes
      if(!inplace){
        setAttrib(out_, mkString("uplo"), type_);
        setAttrib(out_, mkString("n"), Rf_ScalarInteger(nrow));
        setAttrib(out_, mkString("include.diag"), diag_);
        setAttrib(out_, mkString("byrow"), byrow_);
      }

    }else{ // Full matrix or a triangular unpacked (to full)
      if(scale){
        //Rprintf(" Making pointwise multiplication by scalar 'a'...\n");
        for(j=0; j<ncol; j++){
          F77_NAME(dscal)(&nrow, &a, out + nrow*j, &inc1);
        }
      }

      //Rprintf(" Making shifting ...\n");
      switch(case_set){
        case 0: // nirow=0 && nicol=0
          for(j=0; j<ncol; j++){
            out[nrow*j + j] += lambda[nlambda == 1 ? 0 : j];
          }
        break;

        case 1: // nirow>0 && nicol=0
          for(j=0; j<ncol; j++){
            for(i=0; i<nrow; i++){
              if(j == irow[i]){
                out[nrow*j + i] += lambda[nlambda == 1 ? 0 : j];
              }
            }
          }
        break;

        case 2: // nirow=0 && nicol>0
          for(j=0; j<ncol; j++){
            for(i=0; i<nrow; i++){
              if(icol[j] == i){
                out[nrow*j + i] += lambda[nlambda == 1 ? 0 : icol[j]];
              }
            }
          }
        break;

        case 3: // nirow>0 && nicol>0
          for(j=0; j<ncol; j++){
            for(i=0; i<nrow; i++){
              if(icol[j] == irow[i]){
                out[nrow*j + i] += lambda[nlambda == 1 ? 0 : icol[j]];
              }
            }
          }
        break;
      }
    }

    UNPROTECT(nprotect);

    if(inplace){
      return(R_NilValue);
    }else{
      return(out_);
    }
}

//====================================================================

SEXP R_penalize_cov_kronecker(SEXP na_, SEXP a_, SEXP nA_, SEXP A_, SEXP lambda_,
                              SEXP irow_, SEXP icol_, SEXP drop_,
                              SEXP type_, SEXP diag_, SEXP byrow_)
{
    long long i, j, tmp;
    int *posa, *posA;

    int nA=INTEGER_VALUE(nA_);
    int na=INTEGER_VALUE(na_); // should be equal to nlambda
    int byrow=asLogical(byrow_);
    int diag=asLogical(diag_);
    int drop=asLogical(drop_);
    long long nirow=Rf_isNull(irow_) ? 0 : XLENGTH(irow_);
    long long nicol=Rf_isNull(icol_) ? 0 : XLENGTH(icol_);

    int type=strcmp(CHAR(asChar(type_)),"full")==0 ? 0 : strcmp(CHAR(asChar(type_)),"upper")==0 ? 1 : 2;

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(a_ = AS_NUMERIC(a_));
    double *a = NUMERIC_POINTER(a_);

    PROTECT(lambda_=AS_NUMERIC(lambda_));
    double *lambda=NUMERIC_POINTER(lambda_);

    PROTECT(irow_=AS_INTEGER(irow_));
    int *irow=INTEGER_POINTER(irow_);

    PROTECT(icol_=AS_INTEGER(icol_));
    int *icol=INTEGER_POINTER(icol_);

    long long nmap = (long long)na*nA;
    long long nrow = nirow == 0 ? nmap : nirow;
    long long ncol = nicol == 0 ? nmap : nicol;

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
    double *out=NUMERIC_POINTER(out_);

    int case_set = 1*(nirow>0) + 2*(nicol>0);

    posa=(int *) R_alloc(nmap, sizeof(int)); // For a
    posA=(int *) R_alloc(nmap, sizeof(int)); // For A
    tmp = 0;
    for(i=0; i<na; i++){
      for(j=0; j<nA; j++){
        posa[tmp] = i;
        posA[tmp] = j;
        tmp++;
      }
    }

    //Rprintf(" Making kronecker multiplication by matrix 'a'...\n");
    if(type == 0){  // Full matrix
      make_kronecker_full_full(case_set, nrow, ncol, na, na, a, nA, nA, A,
                               posa, posA, posa, posA, irow, icol, out);

    }else{
      if(((type == 1)&&!byrow) || ((type == 2)&&byrow))
      {
          make_kronecker_tri_full1(case_set, nrow, ncol, nA, A, na, a, diag,
                                   posA, posa, posA, posa, irow, icol, out);


      }else{
        if(((type == 1)&&byrow) || ((type == 2)&&!byrow))
        {
            make_kronecker_tri_full2(case_set, nrow, ncol, nA, A, na, a, diag,
                                     posA, posa, posA, posa, irow, icol, out);

        }else{
          Rprintf(" No action performed ...\n");
        }
      }
      if(!diag){
        Rprintf(" Matrix does not contain diagonal data, a diagonal equal to 1 is assumed\n");
      }
    }

    //Rprintf(" Making shifting ...\n");
    // The output is now a full matrix, thus a common shifting code is applied
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            if(posA[j] == posA[i]){
              out[nrow*j + i] += lambda[na*(long long)posa[j]+(long long)posa[i]];
            }
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            if(posA[j] == posA[irow[i]]){
              out[nrow*j + i] += lambda[na*(long long)posa[j]+(long long)posa[irow[i]]];
            }
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            if(posA[icol[j]] == posA[i]){
              out[nrow*j + i] += lambda[na*(long long)posa[icol[j]]+(long long)posa[i]];
            }
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            if(posA[icol[j]] == posA[irow[i]]){
              out[nrow*j + i] += lambda[na*(long long)posa[icol[j]]+(long long)posa[irow[i]]];
            }
          }
        }
      break;
    }

    UNPROTECT(6);

    return(out_);
}
