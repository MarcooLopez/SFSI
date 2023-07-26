#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>
#include "SFSI.h"
//#include "utils_all.c"

//====================================================================
// Functions used by the 'readBinary' routine
//====================================================================
void read_integer(FILE *f, long long start, long long ncol, long long nrow,
                 long long n, long long p,
                 int *X, int varsize, int nsetrow, int *setrow,
                 int nsetcol, int *setcol, int *nerror){
  long long i, j;
  int *line;
  size_t out=0;

  fseeko(f, 0, SEEK_END);
  if((nrow*ncol*(long long)varsize + start)==ftello(f))
  {
    fseeko(f, start, SEEK_SET);

    if(nsetrow>0){
     line = (int *) R_alloc(nrow, sizeof(int));
    }

    for(j=0; j<p; j++){
      if(nsetcol > 0){
        fseeko(f, start + nrow*(long long)varsize*((long long)setcol[j]-1), SEEK_SET);
      }

      if(nsetrow>0){
         out+=fread(line, varsize, nrow, f);
         for(i=0; i<n; i++){
           X[n*j + i] = line[setrow[i]-1];
         }
      }else{
         out+=fread(X + n*j, varsize, nrow, f);
      }
    }
    if(out != (nrow*p)){
      Rprintf("  Error: Something went wrong when reading data from file");
      (*nerror)++;
    }
  }else{
    Rprintf("  Error: The function failed to read data from file\n");
    (*nerror)++;
  }
}

//====================================================================

void read_double(FILE *f, long long start, long long ncol, long long nrow,
                 long long n, long long p,
                 double *X, int varsize, int nsetrow, int *setrow,
                 int nsetcol, int *setcol, int *nerror){
  long long i, j;
  size_t out=0;
  double *linedouble=NULL;
  float *linefloat=NULL;

  fseeko(f, 0, SEEK_END);
  if((nrow*ncol*(long long)varsize + start)==ftello(f))
  {
    fseeko(f, start, SEEK_SET);

    if(varsize == sizeof(double)){
      if(nsetrow>0){
        linedouble = (double *) R_alloc(nrow, sizeof(double));
      }
    }else{
      linefloat = (float *) R_alloc(nrow, sizeof(float));
    }

    if(varsize == sizeof(double)){
      for(j=0; j<p; j++){
        if(nsetcol > 0){
          fseeko(f, start + nrow*(long long)varsize*((long long)setcol[j]-1), SEEK_SET);
        }

        if(nsetrow>0){
           out+=fread(linedouble, varsize, nrow, f);
           for(i=0; i<n; i++){
             X[n*j + i] = linedouble[setrow[i]-1];
           }
        }else{
           out+=fread(X + n*j, varsize, nrow, f);
        }
      }
    }else{
      for(j=0; j<p; j++){
        if(nsetcol > 0){
          fseeko(f, start + nrow*(long long)varsize*((long long)setcol[j]-1), SEEK_SET);
        }

        out+=fread(linefloat, varsize, nrow, f);
        if(nsetrow>0){
          for(i=0; i<n; i++){
            X[n*j + i] = (double)linefloat[setrow[i]-1];
          }
        }else{
          for(i=0; i<n; i++){
            X[n*j + i] = (double)linefloat[i];
          }
        }
      }
    }
    if(out!=(nrow*p)){
      Rprintf("  Error: Something went wrong when reading data from file");
      (*nerror)++;
    }
  }else{
    Rprintf("  Error: The function failed to read data from file\n");
    (*nerror)++;
  }
}

//====================================================================
// vartype: 1:integer, 2: logical, 3: double
//====================================================================
SEXP R_writeBinFile(SEXP filename_,
                    SEXP nrow_,
                    SEXP ncol_,
                    SEXP X_,
                    SEXP doubleprecision_)
{
    FILE *f=NULL;
    long long i, j;
    long long n, p;
    int varsize, vartype;
    float valuefloat;
    SEXP list;

    int nrow=INTEGER_VALUE(nrow_);
    int ncol=INTEGER_VALUE(ncol_);
    int doubleprecision=asLogical(doubleprecision_);

    varsize=0;  // Initialize value
    vartype=0;  // Initialize value

    f=fopen(CHAR(STRING_ELT(filename_,0)),"wb");
    fwrite(&nrow, sizeof(int), 1, f);
    fwrite(&ncol, sizeof(int), 1, f);

    n = (long long)nrow;
    p = (long long)ncol;
    if(p == 0){
      p = 1;
    }

    if(TYPEOF(X_) == INTSXP || TYPEOF(X_) == LGLSXP)
    {
      vartype = TYPEOF(X_) ==  INTSXP ? 1 : 2;

      PROTECT(X_=AS_INTEGER(X_));
      int *X=INTEGER_POINTER(X_);   // An integer pointer is also used if is FLOAT
      varsize = sizeof(X[0]);
      fwrite(&vartype, sizeof(int), 1, f);
      fwrite(&varsize, sizeof(int), 1, f);

      for(j=0; j<p; j++){
        fwrite(X + n*j, varsize, n, f);
      }

    }else{
      if(TYPEOF(X_) == REALSXP)
      {
        vartype = 3;
        PROTECT(X_=AS_NUMERIC(X_));
        double *X=NUMERIC_POINTER(X_);

        if(doubleprecision)
        {
          varsize = sizeof(X[0]);
          fwrite(&vartype, sizeof(int), 1, f);
          fwrite(&varsize, sizeof(int), 1, f);

          for(j=0; j<p; j++){
            fwrite(X + n*j, varsize, n, f);
          }
        }else{
          varsize = sizeof(float);
          fwrite(&vartype, sizeof(int), 1, f);
          fwrite(&varsize, sizeof(int), 1, f);

          for(j=0; j<p; j++){
            for(i=0; i<n; i++){
              valuefloat = X[n*j + i];
              fwrite(&valuefloat, varsize, 1, f);
            }
          }
        }

      }else{
        Rprintf("  File can not be saved with the current type format\n");
      }
    }
    fclose(f);

    PROTECT(list = Rf_allocVector(VECSXP, 4));

    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(vartype));
    SET_VECTOR_ELT(list, 3, ScalarInteger(varsize));

    UNPROTECT(2);

    return(list);
}

//====================================================================
//====================================================================

SEXP R_readBinFile(SEXP filename_,
                   SEXP setrow_,
                   SEXP setcol_)
{
    FILE *f=NULL;
    int *setrow, *setcol;
    int varsize, vartype;
    int intval;
    int ismatrix;
    int nerror=0;
    long long nrow, ncol, n, p;
    size_t out;

    SEXP list;
    SEXP X_=NULL;

    int nsetrow=Rf_isNull(setrow_) ? 0 : XLENGTH(setrow_);
    int nsetcol=Rf_isNull(setcol_) ? 0 : XLENGTH(setcol_);

    PROTECT(setrow_=AS_INTEGER(setrow_));
    setrow=INTEGER_POINTER(setrow_);

    PROTECT(setcol_=AS_INTEGER(setcol_));
    setcol=INTEGER_POINTER(setcol_);

    f=fopen(CHAR(STRING_ELT(filename_,0)),"rb");

    out=fread(&intval, sizeof(int), 1, f);
    nrow = (long long)intval;
    out+=fread(&intval, sizeof(int), 1, f);
    ncol = (long long)intval;
    out+=fread(&vartype, sizeof(int), 1, f);
    out+=fread(&varsize, sizeof(int), 1, f);

    ismatrix = ncol > 0 ? 1 : 0;
    if(ncol == 0){
      ncol = 1;
    }

    if(out < 4){
      Rprintf("  Error: The function failed to read data information\n");
      nerror++;
    }

    // Check if any index is larger than n or p
    if(nsetrow > 0){
       intval =  setrow[imax_integer(nsetrow, setrow)];
       //Rprintf("  Max index row=%d",intval);
       if(intval > nrow){
         Rprintf("  Error in reading row %d: file contains only %d rows\n",intval,nrow);
         nerror++;
       }
    }
    if(nsetcol > 0){
       intval = setcol[imax_integer(nsetcol, setcol)];
       //Rprintf("  Max index column=%d",intval);
       if(intval > ncol){
         Rprintf("  Error in reading column %d: file contains only %d columns\n",intval,ncol);
         nerror++;
       }
    }

    n = nsetrow > 0 ? (long long) nsetrow : nrow;
    p = nsetcol > 0 ? (long long) nsetcol : ncol;
    //Rprintf("  To read: n=%d, p=%d, vartype=%d, varsize=%d\n",n,p,vartype,varsize);

    // vartype: 1:integer, 2: logical, 3: double
    if(nerror == 0)
    {
      if(vartype == 1 || vartype == 2)  // INTEGER|LOGICAL
      {
        if(ismatrix){
          if(vartype == 1){
            X_=PROTECT(Rf_allocMatrix(INTSXP, n, p));
          }else{
            X_=PROTECT(Rf_allocMatrix(LGLSXP, n, p));
          }
        }else{
          if(vartype == 1){
            X_=PROTECT(Rf_allocVector(INTSXP, n));
          }else{
            X_=PROTECT(Rf_allocVector(LGLSXP, n));
          }
        }

        int *X=INTEGER_POINTER(X_);
        read_integer(f, 16, ncol, nrow, n, p, X,
                     varsize, nsetrow, setrow, nsetcol, setcol, &nerror);

      }else{
        if(vartype == 3)  // DOUBLE
        {
          if(ismatrix){
            X_=PROTECT(Rf_allocMatrix(REALSXP, n, p));
          }else{
            X_=PROTECT(Rf_allocVector(REALSXP, n));
          }

          double *X=NUMERIC_POINTER(X_);
          read_double(f, 16, ncol, nrow, n, p, X,
                      varsize, nsetrow, setrow, nsetcol, setcol, &nerror);

        }else{
          Rprintf("  Error: File can not be read with the current type format\n");
        }
      }
    }
    fclose(f);

    PROTECT(list = Rf_allocVector(VECSXP, 6));

    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(vartype));
    SET_VECTOR_ELT(list, 3, ScalarInteger(varsize));
    SET_VECTOR_ELT(list, 4, ScalarInteger(nerror));
    SET_VECTOR_ELT(list, 5, X_);

    UNPROTECT(4);
    return(list);
}
