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
//    Append to a vector v, k elements provided in another vector
//    IN: Vector values (length at least k) to append
//    OUT: Vector v will contain the n+k elements
//====================================================================
double* append_to_vector_double(int n, double *v, int k, double *values)
{
  int i;
  v = R_Realloc(v, n+k, double);

  for(i=0; i<k; i++){
    v[n+i] = values[i];
  }

  return(v);
}

int* append_to_vector_integer(int n, int *v, int k, int *values)
{
  int i;
  v = R_Realloc(v, n+k, int);

  for(i=0; i<k; i++){
    v[n+i] = values[i];
  }

  return(v);
}

//====================================================================
//    Subset a vector v1 by selecting k elements
//    provided in the INDEX vector: any of 0,1,2,...,n-1
//    OUT: Vector v2 contain the k subset elements
//====================================================================
void subset_vector_double(double *v1, double *v2, int k, int *index)
{
  int i;
  for(i=0; i<k; i++){
    v2[i] = v1[index[i]];
  }
}

//====================================================================
//    Delete a row or a column in a n x p matrix M
//    irow: any of 0,1,2, ..., n-1, if irow < 0 no row is deleted
//    icol: any of 0,1,2, ..., p-1, if icol < 0 no column is deleted
//====================================================================
void reduce_matrix(int nrow, int ncol, int irow, int icol, double *M)
{
  long long j;
  long long offset;
  long long nnew, pnew;

  long long nrow0 = (long long)nrow;
  long long ncol0 = (long long)ncol;

  nnew = nrow0;
  if(irow > -1){
    nnew = nrow0-1;
    offset = nrow-irow-1;
    for(j=1; j<ncol; j++){
      memmove(M + nnew*j - offset, M + nrow0*j - offset, nnew*sizeof(double));
    }
    if(offset>0){
      memmove(M + nnew*ncol0 - offset,
              M + nrow0*ncol0 - offset,
              offset*sizeof(double));
    }
  }
  pnew = ncol;
  if(icol > -1){
    pnew = ncol0-1;
    offset = ncol0-icol-1;
    if(offset>0){
      memmove(M + nnew*(long long)icol,
              M + nnew*((long long)icol+1),
              nnew*offset*sizeof(double));
    }
  }
}

//====================================================================
//    Delete k elements from a vector of length n
//    k: number of elements to delete
//    INDEX: vector with any of 0,1,2,...,n-1. Must to be sorted ascending
//====================================================================
void reduce_vector_double(int n, double *v, int k, int *index)
{
  int j;
  int pos;
  int chunk;
  int flag;

  pos = index[0];
  j = 0;
  flag=1;
  while(flag){
    while(((j+1)<k) && ((index[j+1]-index[j])==1)){
      j++;
    }
    if((index[j]+1) < n){
      if((j+1) < k){
        chunk=index[j+1]-index[j]-1;
      }else{
        chunk=n-index[j]-1;
        flag=0;
      }
      memmove(v + pos, v + index[j]+1, chunk*sizeof(double));
      pos+=chunk;
      j++;
    }else{
      flag=0;
    }
  }
}

//====================================================================
// Update coefficients in the p x p-1 matrix R
// that has been formed from a p x p upper-triangular
// matrix whose column k (0, 1,...,p-1) was deleted
//====================================================================
void update_deleted_cols(int p, int k, double *R, int nz, double *z)
{
  double a, b, c, s, tau;
  long long i, j;
  long long pos1, pos2;

  long long p0 = (long long)p;

  double eps = DBL_EPSILON;

  for(i=(long long)k; i<p0-1; i++)
  {
    pos1 = p0*i + i;
    a = R[pos1];
    b = R[pos1 + 1];
    if(fabs(b) > eps) //if(b!=0.0f)
    {
      // Compute the rotation
       if(fabs(b)>fabs(a)){
         tau = -a/b;
         s = 1/sqrt(1 + tau*tau);
         c = s * tau;
       }else{
         tau = -b/a;
         c = 1/sqrt(1 + tau*tau);
         s = c * tau;
       }

       // update r and z
       R[pos1] = c*a - s*b;
       R[pos1 + 1] = s*a + c*b;

       for(j=i+1; j<p0-1; j++)
       {
         pos2 = p0*j + i;
         a = R[pos2];
         b = R[pos2 + 1];
         R[pos2] = c*a - s*b;
         R[pos2 + 1] = s*a + c*b;
       }
       for(j=0; j<nz; j++)
       {
         pos2 = p0*j + i;
         a = z[pos2];
         b = z[pos2 + 1];
         z[pos2] = c*a - s*b;
         z[pos2 + 1] = s*a + c*b;
       }
    }
  }
}

//====================================================================
// Drop row and columm k from a p x p upper-triangular matrix R
// and update resulting (p-1)x(p-1) matrix
// k: any of 0,1,...,p-1
//====================================================================
void downdate_chol(int p, int k, double *R, int nz, double *z)
{
  reduce_matrix(p, p, -1, k, R);

  update_deleted_cols(p, k, R, nz, z);

  reduce_matrix(p, p-1, p-1, -1, R);
}

//====================================================================

SEXP R_lars(SEXP XtX_, SEXP Xty_,
          SEXP eps_, SEXP dfmax_, SEXP scale_,
          SEXP sd_, SEXP isLASSO_,
          SEXP filename_,
          SEXP doubleprecision_, SEXP verbose_)
{
    double *XtX, *Xty;
    double *lambda, *sd, *w, *rhs, *Sign, *a;
    double A, gamhat, zmin;
    int i, j, k;
    int *df, *activeignores;
    int *active, *inactive, *im, *itmp;
    int info, nsteps, inew, nnew, nR;
    int inc1=1;
    int varsize, vartype;
    double Cmax;
    double *B;
    double *covar, *R;
    double *b, *b0, *z, *tmp;
    double value;
    float valuefloat;
    int nprotect=6;
    FILE *f=NULL;
    SEXP list, B_=NULL, lambda_=NULL, df_=NULL;

    int p=XLENGTH(Xty_);
    int dfmax=INTEGER_VALUE(dfmax_);
    int verbose=asLogical(verbose_);
    int scale=asLogical(scale_);
    int isLASSO=asLogical(isLASSO_);
    double eps=NUMERIC_VALUE(eps_);
    int doubleprecision=asLogical(doubleprecision_);
    int save=!Rf_isNull(filename_);

    PROTECT(XtX_=AS_NUMERIC(XtX_));
    XtX=NUMERIC_POINTER(XtX_);

    PROTECT(Xty_=AS_NUMERIC(Xty_));
    Xty=NUMERIC_POINTER(Xty_);

    PROTECT(sd_=AS_NUMERIC(sd_));
    sd=NUMERIC_POINTER(sd_);

    nsteps = 1.25*dfmax;

    lambda=(double *) R_Calloc(nsteps, double);
    df=(int *) R_Calloc(nsteps, int);

    // Allocate memory for B. Allocated memory is set to zero (as in calloc)
    if(save){
      varsize = doubleprecision ? sizeof(double) : sizeof(float);
      vartype = 3;
      f=fopen(CHAR(STRING_ELT(filename_,0)),"wb");
      fwrite(&p, sizeof(int), 1, f);
      fwrite(&nsteps, sizeof(int), 1, f);
      fwrite(&vartype, sizeof(int), 1, f);
      fwrite(&varsize, sizeof(int), 1, f);
    }else{
      B = (double *) R_Calloc(p*nsteps, double);
    }

    rhs=(double *) R_alloc(p, sizeof(double));
    im=(int *) R_alloc(p, sizeof(int));
    covar=(double *) R_alloc(p, sizeof(double));
    a=(double *) R_alloc(p, sizeof(double));
    Sign=(double *) R_alloc(p, sizeof(double));
    inactive=(int *) R_alloc(p, sizeof(int));
    active=(int *) R_alloc(p, sizeof(int));
    activeignores=(int *) R_alloc(p, sizeof(int));
    w =(double *) R_alloc(p, sizeof(double));
    z=(double *) R_alloc(p, sizeof(double));
    b=(double *) R_alloc(p, sizeof(double));
    b0=(double *) R_alloc(p, sizeof(double));  // Output
    R =(double *) R_alloc(p*p, sizeof(double));
    itmp=(int *) R_Calloc(0, int);
    tmp=(double *) R_alloc(p, sizeof(double));

    memset(b, 0, p*sizeof(double));        // Initialize all coefficients to zero
    memset(df, 0, nsteps*sizeof(int));

    // Initialize first column of B to zero
    if(save){
      if(doubleprecision){
        fwrite(b, varsize, p, f);
      }else{  // Cast to float one by one
        valuefloat = 0;
        for(j=0; j<p; j++){
          fwrite(&valuefloat, varsize, 1, f);
        }
      }
    }else{
      memset(B, 0, p*sizeof(double));
    }

    for(j=0; j<p; j++){
      inactive[j] = j;
      im[j] = j;
      z[j] = 1;
    }

    F77_NAME(dcopy)(&p, Xty, &inc1, rhs, &inc1); //  rhs <- Xty
    int nactive = 0;
    int ninactive = p;
    int ncovar = p;
    int nignores = 0;
    int ndrops = 0;
    int extrasteps = (0.01*p)<1 ? 1 : 0.01*p;  // Enlarge by 0.01xp the output's memory

    k=0;
    while(nactive<dfmax && nactive<(p-nignores))
    {
      subset_vector_double(rhs, covar, ninactive, inactive);
      ncovar=ninactive;

      Cmax=fabs(covar[F77_NAME(idamax)(&ninactive, covar, &inc1)-1]);
      if(Cmax < eps*100){
        if(verbose){
          Rprintf(" Max absolute correlation is zero. Exiting...\n");
        }
        break;
      }
      lambda[k++] = Cmax;

      if(ndrops == 0){
        nnew = 0;
        for(i=0; i<ncovar; i++){
          //if(fabs(covar[i]) >= Cmax-eps){
          if(fabs(fabs(covar[i])-Cmax) <= eps){
            itmp=append_to_vector_integer(nnew++, itmp, 1, &i);
          }
        }
        reduce_vector_double(ncovar, covar, nnew, itmp);
        ncovar-=nnew;

        for(i=0; i<nnew; i++){
          inew=inactive[itmp[i]];
          append_to_sorted_vector_integer(nactive+nignores, activeignores, 1, &inew);
          update_chol(p, XtX, nactive, R, inew, active, &eps, tmp, &info);

          if(info==0){
            active[nactive] = inew;
            Sign[nactive] = sign(rhs[inew]);
            nactive++;
            if(verbose){
              Rprintf("--------------------------------------------------------------\n");
              Rprintf(" Step=%5d Lambda=%1.8f nActive=%5d  IN=%5d\n",k,Cmax,nactive,inew+1);
            }
          }else{
            reduce_matrix(nactive+1, nactive+1, nactive, nactive, R);
            nignores++;
            if(verbose){
              Rprintf("  Feature %5d is collinear, dropped for good\n",inew+1);
            }
          }
        }
      }
      memcpy(w, Sign, nactive*sizeof(double));
      backsolvet(nactive, R, w);
      backsolve(nactive, R, w);

      A = 1/sqrt(F77_NAME(ddot)(&nactive, w, &inc1, Sign, &inc1));
      F77_NAME(dscal)(&nactive, &A, w, &inc1);

      // Get inactive subjects
      memcpy(inactive, im, p*sizeof(int));
      reduce_vector_integer(p, inactive, nactive+nignores, activeignores);
      ninactive = p-nactive-nignores;

      if(nactive >= (p-nignores)){
        gamhat = Cmax/A;
      }else{
        // a <- XtX[-c(active,ignores),active] %*% w
        matrix_vector_product_subset(p,p,XtX,w,a,ninactive,inactive,
                                     nactive,active,0,tmp);
        gamhat = Cmax/A;

        for(i=0; i<ncovar; i++){
          value=(Cmax-covar[i])/(A-a[i]);
          if(value>eps && value<gamhat){
            gamhat = value;
          }

          value=(Cmax+covar[i])/(A+a[i]);
          if(value>eps && value<gamhat){
            gamhat = value;
          }
        }
      }

      ndrops = 0;
      if(isLASSO){
        zmin = gamhat;
        for(i=0; i<nactive; i++){ // z <- -B[,j]/w
          tmp[i] = -1*b[active[i]]/w[i];
          if(tmp[i]>eps && tmp[i]<zmin){
            zmin = tmp[i];
          }
        }
        if(zmin < gamhat){
          gamhat = zmin;
          for(i=0; i<nactive; i++){
            if(fabs(tmp[i]-zmin) < eps){
              itmp=append_to_vector_integer(ndrops++, itmp, 1, &i);
            }
          }
        }
      }

      for(i=0; i<nactive; i++){
        b[active[i]] += gamhat*w[i];
      }

      // Update covariances: rhs <- rhs - gamhat*XtX[,active]%*%w
      matrix_vector_product_subset(p,p,XtX,w,a,0,NULL,nactive,active,0,tmp);
      value = -1*gamhat;
      F77_NAME(daxpy)(&p, &value, a, &inc1, rhs, &inc1);

      if(isLASSO && ndrops>0){
        nR = nactive;
        for(i=0; i<ndrops; i++){
          if(verbose){
            Rprintf("--------------------------------------------------------------\n");
            Rprintf(" Step=%5d Lambda=%1.8f nActive=%5d OUT=%5d\n",
                    k+1,Cmax,nactive-i-1,active[itmp[ndrops-i-1]]+1);
          }
          downdate_chol(nR, itmp[ndrops-i-1], R, 1, z);

          b[active[itmp[i]]] = 0;
          for(j=0; j<nactive+nignores; j++){
            if(activeignores[j]==active[itmp[i]]){
              reduce_vector_integer(nactive+nignores, activeignores, 1, &j);
              break;
            }
          }

          // Return dropped active to inactive
          append_to_sorted_vector_integer(ninactive++, inactive, 1, active+itmp[i]);
          nR--;
        }
        reduce_vector_integer(nactive, active, ndrops, itmp);
        reduce_vector_double(nactive, Sign, ndrops, itmp);
        nactive-=ndrops;
      }

      F77_NAME(dcopy)(&p, b, &inc1, b0, &inc1);

      if(scale){
        for(j=0; j<p; j++){
          b0[j] = b0[j]/sd[j];
        }
      }

      df[k]=0;
      for(j=0; j<p; j++){
        if(fabs(b0[j])>0) df[k]++;
      }

      if(save){
        if(doubleprecision){
          fwrite(b0, varsize, p, f);
        }else{  // Cast to float one by one
          for(j=0; j<p; j++){
            valuefloat = b0[j];
            fwrite(&valuefloat, varsize, 1, f);
          }
        }
      }else{
        F77_NAME(dcopy)(&p, b0, &inc1, B+(long long)p*(long long)k, &inc1);
      }

      if(k == nsteps-1){  // Enlarge memory for outputs if needed
        nsteps+=extrasteps;
        lambda = R_Realloc(lambda, nsteps, double);
        df = R_Realloc(df, nsteps, int);
        if(!save){
          B = R_Realloc(B, p*nsteps, double);
        }
      }
    }

    if(dfmax<p){ // Get the next max correlation
      subset_vector_double(rhs, covar, ninactive, inactive);
      lambda[k] = fabs(covar[F77_NAME(idamax)(&ninactive, covar, &inc1)-1]);
    }else{
      lambda[k] = fabs(rhs[F77_NAME(idamax)(&p, rhs, &inc1)-1]);
    }
    k++;

    lambda_=PROTECT(Rf_allocVector(REALSXP, k));
    memcpy(NUMERIC_POINTER(lambda_), lambda, k*sizeof(double));

    df_=PROTECT(Rf_allocVector(INTSXP, k));
    memcpy(INTEGER_POINTER(df_), df, k*sizeof(int));

    if(save){
      fseek(f, 4, SEEK_SET); // Save the final number of solutions
      fwrite(&k, 4, 1, f);
      fclose(f);
      B_ = R_NilValue;
    }else{
      B_ = PROTECT(Rf_allocMatrix(REALSXP, p, k));
      memcpy(NUMERIC_POINTER(B_), B, p*k*sizeof(double));
      R_Free(B);
      nprotect++;
    }

    R_Free(lambda);
    R_Free(df);
    R_Free(itmp);

    // Creating a list with 4 vector elements:
    PROTECT(list = Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(list, 0, B_);
    SET_VECTOR_ELT(list, 1, lambda_);
    SET_VECTOR_ELT(list, 2, df_);

    UNPROTECT(nprotect);

    return(list);
}

//====================================================================
// Calculate betas for a set of decreasing lambdas: l[0]>l[1]>...>l[k]
// The procedure is iteratively done until convergence, i.e., the
// difference between two consecutive solutions is smaller than a threshold.
// It uses the covariance matrix XtX and the covariance vector Xty
// XtX is assumed to be scaled with diagonal equal to one
//
//     p:         Number of beta parameters (dim of XtX)
//     XtX:       Crossprod matrix X'X
//     Xty:       Crossprod vector X'y
//     q:         Number of lambdas
//     lambda:    Vector of lambdas for which the betas will be calculated
//     alpha:     Alpha value in the Elastic-net problem
//     tol:       Maximum value between two consecutive solutions for convergence
//     maxiter:   Number of iterations to run before the updating stops
//
// -----------------------------------------------------------------------------
// The fitted values yHatNotj[i] excluding the contribution from X[,j] is
//                     yHatNotj[i] = sum{k!=j}X[i,k]b[k]
// and can be written for the whole vector yHatNotj = {yHatNotj[i]} as
//                     yHatNotj = Xb - X[,j]b[j]
//
// Fitting the partial residual r[i] = y[i] - yHatNotj[i] to X[i,j] has
// an OLS estimator equal to
//             bOLS[j] = (1/n)sum{i=1:n}[X[i,j](y[i] - yHatNotj[i])]
// which can be written as
//             bOLS[j] = X[,j]'y - XtyHatNotj[j]
// where
//                XtyHatNotj[j] = X[,j]'Xb - X[,j]'X[,j]b[j]
//
// After soft-thresholding the OLS estimator for variable j:
//                       bNew[j] <- S(bOLS[j], L1)/(1+L2)
// the terms XtyHatNotj[k] are updated for all k (different from j) if
// there is a change in delta = bNew[j]-b[j]:
//         XtyHatNotj[k] <- XtyHatNotj[k] + X[,k]'X[,k]*(bNew[j]-b[j])
//
// This will replace the contribution of the current value b[j] by the
// contribution of the updated value bNew[j] in XtyHatNotj
//====================================================================
SEXP R_updatebeta(SEXP XtX_, SEXP Xty_,
                   SEXP lambda_, SEXP alpha_, SEXP tol_,
                   SEXP maxiter_, SEXP dfmax_, SEXP scale_,
                   SEXP sd_, SEXP filename_,
                   SEXP doubleprecision_, SEXP verbose_)
{
    double *lambda, *sd;
    double L1, L2, error;
    long long j;
    int k, iter;
    int *df;
    int varsize, vartype;
    int inc1=1;
    double delta, bNew;
    double *B;
    double *b, *b0, *XtyHatNotj;
    float valuefloat;
    int nprotect=7;
    FILE *f=NULL;
    SEXP list, lambda2_=NULL, df_=NULL, B_=NULL;

    int p=XLENGTH(Xty_);
    int nlambda=XLENGTH(lambda_);
    int maxiter=INTEGER_VALUE(maxiter_);
    int dfmax=INTEGER_VALUE(dfmax_);
    int verbose=asLogical(verbose_);
    int scale=asLogical(scale_);
    double alpha=NUMERIC_VALUE(alpha_);
    double tol=NUMERIC_VALUE(tol_);
    int doubleprecision=asLogical(doubleprecision_);
    int save=!Rf_isNull(filename_);

    long long p0 = (long long)p;

    PROTECT(lambda_=AS_NUMERIC(lambda_));
    lambda=NUMERIC_POINTER(lambda_);

    PROTECT(sd_=AS_NUMERIC(sd_));
    sd=NUMERIC_POINTER(sd_);

    PROTECT(XtX_=AS_NUMERIC(XtX_));
    double *XtX=NUMERIC_POINTER(XtX_);

    PROTECT(Xty_=AS_NUMERIC(Xty_));
    double *Xty=NUMERIC_POINTER(Xty_);

    df=(int *) R_alloc(nlambda, sizeof(int));
    b=(double *) R_alloc(p, sizeof(double));   // Current b[j] values
    b0=(double *) R_alloc(p, sizeof(double));  // Output
    XtyHatNotj=(double *) R_alloc(p, sizeof(double)); // XtyHatNotj = {XtX[,j]b - XtX[j,j]b[j]}

    memset(b, 0, sizeof(double)*p);           // Initialize all b[j] to zero
    memset(XtyHatNotj, 0, sizeof(double)*p);  // Since all b[j] are initially 0, all XtyHatNotj[j] are so
    for(k=0; k<nlambda; k++) df[k] = p;

    // Allocate memory for B. Allocated memory is set to zero (as in calloc)
    if(save){
      varsize = doubleprecision ? sizeof(double) : sizeof(float);
      vartype = 3;
      f=fopen(CHAR(STRING_ELT(filename_,0)),"wb");
      fwrite(&p, sizeof(int), 1, f);
      fwrite(&nlambda, sizeof(int), 1, f);
      fwrite(&vartype, sizeof(int), 1, f);
      fwrite(&varsize, sizeof(int), 1, f);
    }else{
      B = (double *) R_Calloc(p*nlambda, double);
    }

    for(k=0; k<nlambda; k++)
    {
        L1 = alpha*lambda[k];
        L2 = (1-alpha)*lambda[k];
        iter = 0;          // Set iter < maxiter to enter the WHILE
        error = tol + 1.0;   // Set a error > tol to enter the WHILE
        while(iter<maxiter && error>tol)
        {
            iter++;
            //Rprintf(" iter=%d\n",iter);
            error = 0;
            for(j=0; j<p0; j++)
            {
                // varj = XtX[p0*j + j];  // Variance of predictor j
                // bOLS = (Xty[j] - XtyHat[j])/varj;
                bNew = soft_threshold(Xty[j] - XtyHatNotj[j], L1)/(1+L2);

                delta = bNew-b[j];
                //Rprintf(" j=%d. XtyHatNotj=%1.8f bOLS=%1.8f  bNew=%1.8f  delta=%f\n",j+1,XtyHatNotj[j],Xty[j]-XtyHatNotj[j], bNew, delta);
                if(fabs(delta)>0){ // Update only if there is a change
                  // Update XtyHatNotj for all k!=j: XtyHatNotj[k] <- XtyHatNotj[k] + XtX[k,k]*(bNew[j]-b[j])
                  F77_NAME(daxpy)(&p, &delta, XtX + p0*j, &inc1, XtyHatNotj, &inc1);

                  XtyHatNotj[j] -= delta; // Except for k=j. delta*varj
                  if(fabs(delta)>error){
                    error = fabs(delta);
                  }
                  b[j] = bNew;
                }
            }
        }
        if(verbose){
            Rprintf(" lambda[%d]=%1.8f  nIters=%5d  Error=%G\n",k+1,lambda[k],iter,error);
            if(error>tol){
              Rprintf(" Warning: The process did not converge after %d iterations for lambda[%d]=%f\n",maxiter,k+1,lambda[k]);
            }
        }

        F77_NAME(dcopy)(&p, b, &inc1, b0, &inc1);

        if(scale){
          for(j=0; j<p0; j++){
            b0[j] = b0[j]/sd[j];
          }
        }

        df[k] = 0;
        for(j=0; j<p0; j++){
          if(fabs(b0[j])>0) df[k]++;
        }

        if(save){
          if(doubleprecision){
            fwrite(b0, varsize, p, f);
          }else{  // Cast to float one by one
            for(j=0; j<p0; j++){
              valuefloat = b0[j];
              fwrite(&valuefloat, varsize, 1, f);
            }
          }
        }else{
          F77_NAME(dcopy)(&p, b0, &inc1, B + p*k, &inc1);
        }

        if(dfmax<p && df[k]>=dfmax){
          break;
        }
    }

    lambda2_=PROTECT(Rf_allocVector(REALSXP, k));
    memcpy(NUMERIC_POINTER(lambda2_), lambda, k*sizeof(double));

    df_=PROTECT(Rf_allocVector(INTSXP, k));
    memcpy(INTEGER_POINTER(df_), df, k*sizeof(int));

    if(save){
      fseek(f, 4, SEEK_SET); // Save the final number of solutions
      fwrite(&k, 4, 1, f);
      fclose(f);
      B_ = R_NilValue;
    }else{
      B_ = PROTECT(Rf_allocMatrix(REALSXP, p, k));
      memcpy(NUMERIC_POINTER(B_), B, p*k*sizeof(double));
      nprotect++;
    }

    PROTECT(list = Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(list, 0, B_);
    SET_VECTOR_ELT(list, 1, lambda2_);
    SET_VECTOR_ELT(list, 2, df_);

    UNPROTECT(nprotect);

    return(list);
}
