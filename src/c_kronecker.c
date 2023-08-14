#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>
#include "SFSI.h"
//#include "utils.c"

//====================================================================
// Make the kronecker product of matrices A and B
//       A:    Numeric matrix of dimension nrowA x ncolA
//       B:    Numeric matrix of dimension nrowB x ncolB
//    type:    Type of matrix B, one of
//             0 = full. For a nxn matrix
//             1 = upper packed triangular
//             2 = lower packed triangular
//             If lower(upper), B is assumed symmetric and
//             must contain n(n+1)/2 or n(n-1)/2 entries
//             depending on whether the diagonal is included
//   totri:    Should the output returned as triangular
//====================================================================
SEXP R_kronecker(SEXP nrowA_, SEXP ncolA_, SEXP A_,
                 SEXP nrowB_, SEXP ncolB_, SEXP B_,
                 SEXP irow_, SEXP icol_, SEXP drop_,
                 SEXP type1_, SEXP byrow1_, SEXP diag1_,
                 SEXP type2_, SEXP byrow2_, SEXP diag2_,
                 SEXP unpack_, SEXP verbose_)
{
    long long tmp, i, j;
    int *posArow, *posBrow, *posAcol, *posBcol;

    int nrowA=INTEGER_VALUE(nrowA_);
    int ncolA=INTEGER_VALUE(ncolA_);
    int nrowB=INTEGER_VALUE(nrowB_);
    int ncolB=INTEGER_VALUE(ncolB_);
    int byrow1=asLogical(byrow1_);
    int byrow2=asLogical(byrow2_);
    int diag1=asLogical(diag1_);
    int diag2=asLogical(diag2_);
    int drop=asLogical(drop_);
    int unpack=asLogical(unpack_);
    int verbose=asLogical(verbose_);
    long long nirow=Rf_isNull(irow_) ? 0 : XLENGTH(irow_);
    long long nicol=Rf_isNull(icol_) ? 0 : XLENGTH(icol_);

    int type1=strcmp(CHAR(asChar(type1_)),"full")==0 ? 0 : strcmp(CHAR(asChar(type1_)),"upper")==0 ? 1 : 2;
    int type2=strcmp(CHAR(asChar(type2_)),"full")==0 ? 0 : strcmp(CHAR(asChar(type2_)),"upper")==0 ? 1 : 2;

    PROTECT(A_ = AS_NUMERIC(A_));
    double *A = NUMERIC_POINTER(A_);

    PROTECT(B_ = AS_NUMERIC(B_));
    double *B = NUMERIC_POINTER(B_);

    PROTECT(irow_=AS_INTEGER(irow_));
    int *irow=INTEGER_POINTER(irow_);

    PROTECT(icol_=AS_INTEGER(icol_));
    int *icol=INTEGER_POINTER(icol_);

    long long nrow = nirow == 0 ? (long long)nrowA*nrowB : nirow;
    long long ncol = nicol == 0 ? (long long)ncolA*ncolB : nicol;

    int upper1 = type1==1 ? 1 : 0;
    int upper2 = type2==1 ? 1 : 0;
    int case_set = 1*(nirow>0) + 2*(nicol>0);   // 0=(nirow=0)&(nicol=0)
    int case_type = 1*(type1>0) + 2*(type2>0);  // 0=(type1=0)&(type2=0)

    int totri = (case_type==3) && all_equal(nirow,irow,nicol,icol);

    if(totri){
      if(unpack){
        totri = 0;
        if(verbose){
          Rprintf(" Matrix was unpacked to a full matrix\n");
        }
      }
    }else{
      if(verbose & (case_type==3)){
        Rprintf(" Matrix was unpacked to a full matrix when selecting different rows/cols\n");
      }
    }

    SEXP out_;
    // If both are triangular, output is triangular too
    // when irow=icol or both NULL
    if(totri){
      tmp = diag1 ? ncol*(ncol+1)/2 : ncol*(ncol-1)/2;
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
    double *out=NUMERIC_POINTER(out_);

    if(!((case_set == 0) && (case_type==0))){
      posAcol=(int *) R_alloc((long long)ncolA*ncolB, sizeof(int));
      posBcol=(int *) R_alloc((long long)ncolA*ncolB, sizeof(int));

      tmp = 0;
      for(i=0; i<ncolA; i++){
        for(j=0; j<ncolB; j++){
          posAcol[tmp] = i;
          posBcol[tmp] = j;
          tmp++;
        }
      }

      if((nrowA==ncolA) && (nrowB==ncolB)){
        // If both are squared matrices, do not repeat the indices
        posArow = posAcol;
        posBrow = posBcol;
      }else{
        posArow=(int *) R_alloc((long long)nrowA*nrowB, sizeof(int));
        posBrow=(int *) R_alloc((long long)nrowA*nrowB, sizeof(int));

        tmp = 0;
        for(i=0; i<nrowA; i++){
          for(j=0; j<nrowB; j++){
            posArow[tmp] = i;
            posBrow[tmp] = j;
            tmp++;
          }
        }
      }
    }

    //Rprintf(" Making kronecker multiplication...\n");
    switch(case_type){
      case 0: // (type1=0)&&(type2=0). Both A and B are full
        make_kronecker_full_full(case_set, nrow, ncol, nrowA, ncolA, A, nrowB, ncolB, B,
                                 posArow, posBrow, posAcol, posBcol, irow, icol, out);
      break;

      case 1: // (type1>0)&&(type2=0). A is triangular, B is full
        if((upper1&&!byrow1)||(!upper1&&byrow1))
        {
            make_kronecker_tri_full1(case_set, nrow, ncol, nrowA, A, nrowB, B, diag1,
                                     posArow, posBrow, posAcol, posBcol, irow, icol, out);

        }else{
          if((upper1&&byrow1)||(!upper1&&!byrow1)){
            make_kronecker_tri_full2(case_set, nrow, ncol, nrowA, A, nrowB, B, diag1,
                                     posArow, posBrow, posAcol, posBcol, irow, icol, out);
          }else{
            Rprintf(" No action performed ...\n");
          }
        }
      break;

      case 2: // (type1=0)&&(type2>0). A is full, B is triangular
        if((upper2&&!byrow2)||(!upper2&&byrow2))
        {
            make_kronecker_tri_full1(case_set, nrow, ncol, nrowB, B, nrowA, A, diag2,
                                     posBrow, posArow, posBcol, posAcol,
                                     irow, icol, out);

        }else{
          if((upper2&&byrow2)||(!upper2&&!byrow2)){
            make_kronecker_tri_full2(case_set, nrow, ncol, nrowB, B, nrowA, A, diag2,
                                     posBrow, posArow, posBcol, posAcol,
                                     irow, icol, out);
          }else{
            Rprintf(" No action performed ...\n");
          }
        }
      break;

      case 3: // (type1>0)&&(type2>0). Both A and B are triangular. Both either with diagonal or not
        if(((upper1&&!byrow1)||(!upper1&&byrow1)) && ((upper2&&!byrow2)||(!upper2&&byrow2)))
        {
           if(totri){
              make_kronecker_tri_tri1_totri(ncol, ncolA, A, ncolB, B, diag1,
                                            posAcol, posBcol, nicol, icol, out);
           }else{
              make_kronecker_tri_tri1(case_set, nrow, ncol, nrowA, A, nrowB, B, diag1,
                                      posArow, posBrow, posAcol, posBcol, irow, icol, out);
           }

        }else{
           if(((upper1&&byrow1)||(!upper1&&!byrow1)) && ((upper2&&byrow2)||(!upper2&&!byrow2)))
           {
             if(totri){
                make_kronecker_tri_tri2_totri(ncol, ncolA, A, ncolB, B, diag1,
                                              posAcol, posBcol, nicol, icol, out);
             }else{
               make_kronecker_tri_tri2(case_set, nrow, ncol, nrowA, A, nrowB, B, diag1,
                                       posArow, posBrow, posAcol, posBcol, irow, icol, out);
             }
           }else{
             Rprintf(" No action performed ...\n");
           }
        }
      break;
    }

    if(totri){
      setAttrib(out_, mkString("uplo"), type1_);
      setAttrib(out_, mkString("n"), Rf_ScalarInteger(nrow));
      setAttrib(out_, mkString("include.diag"), diag1_);
      setAttrib(out_, mkString("byrow"), byrow1_);
    }

    UNPROTECT(5);

    return(out_);
}
