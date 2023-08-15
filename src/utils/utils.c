#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <Rconfig.h>
#ifndef FCONE
# define FCONE
#endif

//====================================================================
// Return 1 if the all elements in set1 are
// equal to those in set2, of if both are NULL
// Return 0 otherwise
//====================================================================
int all_equal(long long n1, int *set1, long long n2, int *set2){
  int out;
  long long i;
  if((n1 + n2) == 0){
    out = 1;
  }else{
    if(n1 == n2){
      out = 1;
      for(i=0; i<n1; i++){
        if(set1[i] != set2[i]){
          out = 0;
          break;
        }
      }
    }else{
      out = 0;
    }
  }
  return(out);
}

//====================================================================
// Return the position of the maximum element of 'x'
//   x: an integer vector of length 'n'
//====================================================================
int imax_integer(int n, int *x){
  int index=0;
  int i;
  for(i=0; i<n; i++){
    if(x[i] > x[index]){
      index = i;
    }
  }
  return(index);
}

//====================================================================
// Soft-thresholding operator for Coordinate Descend Algorithm
//====================================================================
double soft_threshold(double z, double penalty)
{
  if(penalty < fabs(z)){
    return(sign(z)*(fabs(z)-penalty));
  }else{
    return(0);
  }
}

//====================================================================
// Matrix-vector multiplication y = aAx or y=aA'x
// First ncol (or nrow if transpose=1) elements from 'x' are referenced
// OUT: Length of y must be the at least nrow
//      or ncol if transpose
// The level-two BLAS routine 'dgemv' is used
//====================================================================
void matrix_vector_product(int nrow,
                           int ncol,
                           double *a,
                           double *A,
                           double *x,
                           int incx,
                           double *y,
                           int transpose)
{
  char *trans;
  double beta=0;
  int inc1=1;

  trans = transpose ? "T" : "N";

  F77_NAME(dgemv)(trans, &nrow, &ncol, a, A, &nrow, x, &incx,
                  &beta, y, &inc1 FCONE);
}

//====================================================================
// Matrix-matrix multiplication C = AB
// IN: nrow(B) must be the same as ncol(A)
// OUT: Length of C must be the at least nrow(A)*ncol(B)
// The level-3 BLAS routine 'dgemm' is used
//====================================================================
void matrix_matrix_product(int nrowA, int ncolA, int ncolB, double *A,
                           double *B, double *C)
{
  char *trans = "N";
  double beta=0;
  double one=1;

  F77_NAME(dgemm)(trans,trans,&nrowA,&ncolB,&ncolA,&one,A,&nrowA,
                  B,&ncolA,&beta,C,&nrowA FCONE FCONE);
}

//====================================================================
// crossprod C = A'B
// nrow(A) must be the same as nrow(B)
// OUT: Matrix C must contain at least ncolA x ncolB elements
// The level-two BLAS routine 'dgemm' is used
//====================================================================
void crossproduct(int nrow, int ncolA, int ncolB,
                  double *A, double *B, double *C)
{
  char *trans1 = "T";
  char *trans2 = "N";
  double beta=0;
  double one=1;

  F77_NAME(dgemm)(trans1,trans2,&ncolA,&ncolB,&nrow,&one,A,&nrow,
                  B,&nrow,&beta,C,&ncolA FCONE FCONE);
}

//====================================================================
// tcrossprod C = A B'
// ncol(A) must be the same as ncol(B)
// OUT: C must contain at least nrowA x nrowB elements
// If BLAS, the level-two BLAS routine 'dgemm' is used
//====================================================================
void tcrossproduct(int nrowA, int ncol, int nrowB,
                   double *A, double *B, double *C)
{
  char *trans1 = "N";
  char *trans2 = "T";
  double beta=0;
  double one=1;

  F77_NAME(dgemm)(trans1,trans2,&nrowA,&nrowB,&ncol,&one,A,&nrowA,
                  B,&nrowB,&beta,C,&nrowA FCONE FCONE);
}

//====================================================================
// Matrix-vector multiplication y = Ax or y=A'x
// First ncol (or nrow if transpose=1) elements from 'x' are referenced
// IN: Matrix A is squared matrix of dimension n=ncol=nrow
//     upperlower: =1 if A is upper or =2 if A is lower triangular
// OUT: Length of y must be the at least nrow(A)
//      or ncol(A) if transpose
// The level-two BLAS routine 'dtrmv' is used
//====================================================================
void matrix_vector_product_tri(int n,
                               double *A,
                               double *x,
                               int incx,
                               double *y,
                               int transpose,
                               int upperlower)
{
  int inc1=1;
  char *trans, *uplo;
  char *diag = "N";

  trans = transpose ? "T" : "N";
  uplo = upperlower==1 ? "U" : "L";

  F77_NAME(dcopy)(&n, x, &incx, y, &inc1); // copy x to y

  F77_NAME(dtrmv)(uplo, trans, diag, &n, A, &n, y, &inc1 FCONE FCONE FCONE);
}

//====================================================================
// tcrossprod C = A B'
// A,B must be upper triangular squared matrices of dimension n
// OUT: C must contain at least n x n elements
//====================================================================
void tcrossproduct_tri(int n, double *A, double *B, double *C)
{
   long long i;

   // Loop over all rows of B
   for(i=0; i<n; i++){
     // Obtain the i column of C: C[,i] = A*B[i,]'
     matrix_vector_product_tri(n, A,
                               B + i, n,
                               C + n*i,
                               0, 1);
   }
}

//====================================================================
// crossprod C = A'D B,  with D = {d} being a diagonal matrix
// of dimension LENGTH(d) = nrow(A) = nrow(B)
// nrow(A) must be the same as nrow(B)
// OUT: Matrix C must contain at least ncolA x ncolB elements
// WORK is a working space of dimension of at least nrowA
//====================================================================
void crossproduct_scale(int nrowAB, int ncolA, int ncolB,
                        double *A, double *d, double *B,
                        double *C, double *work)
{
  long long i, j;
  double one=1;

  // Loop over all columns of B
  for(j=0; j<ncolB; j++){
    // Obtain the j column of B scaled by d: w =  B[,j]*d
    for(i=0; i<nrowAB; i++){
      work[i] = d[i]*B[nrowAB*j + i];
    }

    // Obtain the j colum of C: C[,j] = A'w
    matrix_vector_product(nrowAB, ncolA, &one, A, work, 1,
                          C + ncolA*j, 1);
  }
}

//====================================================================
// tcrossprod C = A D B' with D = {d} being a diagonal matrix
// of dimension LENGTH(d) = ncol(A) = ncol(B)
// ncol(A) must be the same as ncol(B)
// OUT: Matrix B must contain at least nrowA x nrowB elements
// WORK is a working space of dimension of at least ncol(A)
//====================================================================
void tcrossproduct_scale(int nrowA, int ncolAB, int nrowB,
                         double *A, double *d, double *B,
                         double *C, double *work)
{
  long long i, j;
  double one=1;

  // Loop over all rows of B
  for(i=0; i<nrowB; i++){
    // Obtain the i row of B scaled by d: w =  B[i,]*d
    for(j=0; j<ncolAB; j++){
      work[j] = d[j]*B[nrowB*j + i];
    }

    // Obtain the i column of C: C[,i] = Aw
    matrix_vector_product(nrowA, ncolAB, &one, A, work, 1,
                          C + nrowA*i, 0);
  }
}

//====================================================================
//   Subset from an A matrix by selecting n elements
//   provided in INDEX vector being any of:
//     0,1,2,...,ncol(A) of row k      if margin=1.
//     0,1,2,...,nrow(A) of column k   if margin=2 or.
//   If n>0 and index=NULL, the first n entries of k row (if margin=1)
//   or column (if margin=2) are retrieved.
//   Output is stored in first n elements of x vector
//====================================================================
void slice_matrix(int nrow, double *A, double *x,
                  int n, int *index,
                  int k,
                  int margin)
{
  long long i;
  int inc1=1;
  //int one=1;

  if(margin == 1){
    if(index == NULL){  // First n entries of row k
      F77_NAME(dcopy)(&n, A+k, &nrow, x, &inc1);
    }else{
      for(i=0; i<n; i++){
        x[i] = A[nrow*(long long)index[i] + (long long)k];
      }
    }

  }else{
    if(index == NULL){  // First n entries of column k
      F77_NAME(dcopy)(&n, A + nrow*(long long)k, &inc1, x, &inc1);
    }else{
      for(i=0; i<n; i++){
        x[i] = A[nrow*(long long)k + (long long)index[i]];
      }
    }
  }
}

//====================================================================

void resize_matrix(int nrow, int ncol,
                   int nrownew, int ncolnew,
                   double *A)
{
  long long j, offset;

  int m = ncolnew < ncol ? ncolnew : ncol; // min(pnew,p)

  if(nrownew < nrow){
    for(j=1; j<m; j++){
      //Rprintf("    Copying %d elements from %d to %d columns\n", m, n*j, nnew*j);
      memmove(A + nrownew*j, A + nrow*j, (long long)nrownew*sizeof(double));
    }
  }

  if(nrownew > nrow){
    offset = nrownew - nrow;
    if(ncol > 1){
      //  Shifting elements starting at last column and filling with zeros
      for(j=m-1; j>=1; j--){
        memmove(A + nrownew*j, A + nrow*j, (long long)nrow*sizeof(double));
        memset(A + nrownew*j + nrow, 0, offset*sizeof(double)); // set 0 extra rows
      }
    }
    memset(A + nrow, 0, offset*sizeof(double));
  }
  if(ncolnew > ncol){
    memset(A + nrownew*(long long)ncol, 0, (ncolnew-ncol)*(long long)nrownew*sizeof(double));
  }
}

//====================================================================
// Solves the system equations
//       Ax =  b  (backsolve)
//       A'x = b  (backsolvet)
// for A being a upper (or right) n x n triangular matrix
// vector b can be passed and solution x is wrriten on it
// The level-1 BLAS routine 'dtrsv' is used
//====================================================================
void backsolvet(int n, double *A, double *b)
{
  char *uplo = "U";
  char *trans = "T";
  char *diag = "N";
  int inc1=1;

  F77_NAME(dtrsv)(uplo, trans, diag, &n, A, &n, b, &inc1 FCONE FCONE FCONE);
}

void backsolve(int n, double *A, double *b)
{
  char *uplo = "U";
  char *trans = "N";
  char *diag = "N";
  int inc1=1;

  F77_NAME(dtrsv)(uplo, trans, diag, &n, A, &n, b, &inc1 FCONE FCONE FCONE);
}

//====================================================================
// Update the Cholesky factorization A = R'R where R is
// an upper-triangular matrix.
// The current decomposition A[index,index] = R R'
// is updated to when adding the new element k to
// correspond to A[(index,k),(index,k)] = R R'
// where the new R is of dimensions length(index)+1
//====================================================================
void update_chol(int n, double *A, int nR, double *R,
                 int k, int *index, double *eps,
                 double *work, int *info)
{
  if(nR==0){
    R[0] = sqrt(A[(long long)n*(long long)k + (long long)k]);
    *info = 0;
  }else{
    int inc1=1;
    double rpp;
    //Rprintf("    Obtaining the kth vector: v=A[k,index] ...\n");
    slice_matrix(n,A,work,nR,index,k,1);

    //Rprintf("    Obtaining backsolve solution Lv=r...\n");
    backsolvet(nR, R, work);
    rpp = A[(long long)n*(long long)k + (long long)k] - F77_NAME(ddot)(&nR, work, &inc1, work, &inc1);

    if(rpp <= (*eps)){
      rpp = *eps;
      *info = 1;
    }else{
      rpp=sqrt(rpp);
      *info = 0;
    }
    //Rprintf("    Enlarging R ...\n");
    resize_matrix(nR, nR, nR+1, nR+1, R);

    //Rprintf("    Copying new values to enlarged R ...\n");
    F77_NAME(dcopy)(&nR, work, &inc1, R + ((long long)nR+1)*(long long)nR, &inc1);
    R[((long long)nR+1)*(long long)nR + (long long)nR] = rpp;
  }
}

//====================================================================
//====================================================================

void invert_matrix(int n, double *A, double *Ainv, double *eps, double *work)
{
  int k, info;

  // Get the Cholesky factorization stepwise A = R'R
  info = 0;
  for(k=0; k<n; k++){
    update_chol(n, A, k, Ainv, k, NULL, eps, work, &info);
    if(info != 0){
      Rprintf(" Row/column %d is collinear. No inverse is calculated\n",k+1);
      break;
    }
  }

  if(info == 0){
    long long j;

    // Invert R matrix by columwise backsolve: Rinv[,j]
    memset(work, 0, (long long)n*(long long)n*sizeof(double));
    for(j=0; j<n; j++){
      work[(long long)n*j + j] = 1;
      backsolve(n, Ainv, work + (long long)n*j);
    }
    // Get the inverse: Ainv <- Rinv*(Rinv)'
    tcrossproduct_tri(n, work, work, Ainv);
    //tcrossproduct(n, n, n, work, work, Ainv);
  }
}

//====================================================================
// Matrix-vector multiplication y = A[irow,icol]x or y=A[irow,icol]'x
// First nicol (or nirow if transpose=1) elements from 'x' are referenced
// OUT: Length of y must be the at least nrow
//      or ncol if transpose
//====================================================================
void matrix_vector_product_subset(int nrow,
                                  int ncol,
                                  double *A,
                                  double *x,
                                  double *y,
                                  int nirow, int *irow,
                                  int nicol, int *icol,
                                  int transpose,
                                  double *work)
{

  int i, j;
  int inc1=1;
  double one=1;

  int task = 1*(nirow>0) + 2*(nicol>0);

  switch(task){
    case 0: // nirow==0 && nicol==0
      matrix_vector_product(nrow,ncol,&one,A,x,inc1,y,transpose);
    break;

    case 1: // nirow>0 && nicol==0
      if(transpose){
        for(j=0; j<ncol; j++){
          slice_matrix(nrow,A,work,nirow,irow,j,2);  // Get the jth rowwise indexed column of A
          y[j] = F77_NAME(ddot)(&nirow,work,&inc1,x,&inc1);
        }
      }else{
        for(i=0; i<nirow; i++){
          y[i] = F77_NAME(ddot)(&ncol,A+irow[i],&nrow,x,&inc1);
        }
      }
    break;

    case 2: // nirow==0 && nicol>0
      if(transpose){
        for(j=0; j<nicol; j++){
          y[j] = F77_NAME(ddot)(&nrow,A + nrow*(long long)icol[j],&inc1,x,&inc1);
        }
      }else{
        for(i=0; i<nrow; i++){
          slice_matrix(nrow,A,work,nicol,icol,i,1);  // Get the ith indexed row of A
          y[i] = F77_NAME(ddot)(&nicol,work,&inc1,x,&inc1);
        }
      }
    break;

    case 3: // nirow>0 && nicol>0
      if(transpose){
        for(j=0; j<nicol; j++){
          slice_matrix(nrow,A,work,nirow,irow,icol[j],2); // Get the jth rowwise indexed column of A
          y[j] = F77_NAME(ddot)(&nirow,work,&inc1,x,&inc1);
        }
      }else{
        for(i=0; i<nirow; i++){
          slice_matrix(nrow,A,work,nicol,icol,irow[i],1); // Get the ith rowwise indexed row of A
          y[i] = F77_NAME(ddot)(&nicol,work,&inc1,x,&inc1);
        }
      }
    break;
  }
}

//====================================================================
//    Append to a vector v, k elements provided in another
//    vector and keep ordered values: v[0]<v[1]<...<v[n+k-1]
//        values: vector whose length at least k
//        k: number of elements to be aappend to vector v
//====================================================================
void append_to_sorted_vector_integer(int n,
                                     int *v,
                                     int k,
                                     int *values)
{
  int i, j, start;

  if(n==0){
    v[0] = values[0];
    start=1;
  }else{
    start=0;
  }

  for(i=start; i<k; i++){
    j=0;
    while(v[j]<values[i]){
      j++;
      if(j==(n+i)){
        break;
      }
    }
    if(j < (n+i)){
      memmove(v + j+1, v + j, (n+i-j)*sizeof(int));
    }
    v[j]=values[i];
  }
}

//====================================================================
//    Delete k elements from a vector of length n
//    k: number of elements to delete
//    INDEX: vector with any of 0,1,2,...,n-1. Must to be sorted ascending
//====================================================================
void reduce_vector_integer(int n,
                           int *v,
                           int k,
                           int *index)
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
      memmove(v + pos, v + index[j]+1, chunk*sizeof(int));
      pos+=chunk;
      j++;
    }else{
      flag=0;
    }
  }
}

//====================================================================
//    Retrieve the element A[i,j] from a matrix A that is stored
//    as a triangular (uplo='upper' or uplo='lower') matrix stacked
//    row-wise (byrow=TRUE) or column-wise (byrow=FALSE).
//    It is assumed that A is symmetric, therefore A[i,j] = A[j,i]
//
//    'get_tri_diag1' and 'get_tri_nodiag1' works for cases
//    (uplo='upper')&(byrow=FALSE) or (uplo='lower')&(byrow=TRUE)
//
//    'get_tri_diag2' and 'get_tri_nodiag2' works for cases
//    (uplo='upper')&(byrow=TRUE) or (uplo='lower')&(byrow=FALSE)
//====================================================================

double ij_tri_diag1(double *A, long long i, long long j)
{
  if(i <= j){
    return(A[j*(j+1)/2 + i]);
  }else{
    return(A[i*(i+1)/2 + j]);
  }
}

//====================================================================

double ij_tri_nodiag1(double *A, long long i, long long j)
{
  if(i < j){
    return(A[j*(j-1)/2 + i]);
  }else{
    if(i > j){
      return(A[i*(i-1)/2 + j]);
    }else{
      return(1);
    }
  }
}

//====================================================================

double ij_tri_diag2(int nA, double *A, long long i, long long j)
{
  if(i >= j){
    return(A[nA*j - j*(j-1)/2 + i - j]);
  }else{
    return(A[nA*i - i*(i-1)/2 + j - i]);
  }
}

//====================================================================

double ij_tri_nodiag2(int nA, double *A, long long i, long long j)
{
  if(i > j){
    return(A[nA*j - j*(j+1)/2 + i - j - 1]);
  }else{
    if(i < j){
      return(A[nA*i - i*(i+1)/2 + j - i - 1]);
    }else{
      return(1);
    }
  }
}

//====================================================================
// The folowing five functions are used to create a kronecker
// between different combinations of full/triangular matrices.
// Four cases of setting rows/columns are considered:
//     case_set=0 : nirow=0 && nicol=0
//     case_set=1 : nirow>0 && nicol=0
//     case_set=2 : nirow=0 && nicol>0
//     case_set=3 : nirow>0 && nicol>0
//====================================================================
void make_kronecker_full_full(int case_set, long long nrow, long long ncol,
                              int nrowA, int ncolA, double *A,
                              int nrowB, int ncolB, double *B,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out)
{
  long long i, j, k;
  int inc1=1;
  double a;
  double *Aj, *Bj;

  switch(case_set){
    case 0: // nirow=0 && nicol=0
      for(j=0; j<ncolA; j++){
        for(i=0; i<nrowA; i++){
          a = A[nrowA*j+i];
          for(k=0; k<ncolB; k++){
            F77_NAME(dcopy)(&nrowB, B+nrowB*k, &inc1, out+nrowA*nrowB*(ncolB*j+k) + nrowB*i, &inc1);
            F77_NAME(dscal)(&nrowB, &a, out+nrowA*nrowB*(ncolB*j+k) + nrowB*i, &inc1);
          }
        }
      }
    break;

    case 1: // nirow>0 && nicol=0
      for(j=0; j<ncol; j++){
        Aj = A + nrowA*(long long)posAcol[j];
        Bj = B + nrowB*(long long)posBcol[j];
        for(i=0; i<nrow; i++){
          out[nrow*j + i] = Aj[posArow[irow[i]]]*Bj[posBrow[irow[i]]];
        }
      }
    break;

    case 2: // nirow=0 && nicol>0
      for(j=0; j<ncol; j++){
        Aj = A + nrowA*(long long)posAcol[icol[j]];
        Bj = B + nrowB*(long long)posBcol[icol[j]];
        for(i=0; i<nrow; i++){
          out[nrow*j + i] = Aj[posArow[i]]*Bj[posBrow[i]];
        }
      }
    break;

    case 3: // nirow>0 && nicol>0
      for(j=0; j<ncol; j++){
        Aj = A + nrowA*(long long)posAcol[icol[j]];
        Bj = B + nrowB*(long long)posBcol[icol[j]];
        for(i=0; i<nrow; i++){
          out[nrow*j + i] = Aj[posArow[irow[i]]]*Bj[posBrow[irow[i]]];
        }
      }
    break;
  }
}

//====================================================================

void make_kronecker_tri_full1(int case_set, long long nrow, long long ncol,
                              int nrowA, double *A, int nrowB, double *B, int diag,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out)
{
  long long i, j;
  double *Bj;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, posArow[i], posAcol[j])*Bj[posBrow[i]];
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, posArow[irow[i]], posAcol[j])*Bj[posBrow[irow[i]]];
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, posArow[i], posAcol[icol[j]])*Bj[posBrow[i]];
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, posArow[irow[i]], posAcol[icol[j]])*Bj[posBrow[irow[i]]];
          }
        }
      break;
    }

  }else{ // No diagonal
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, posArow[i], posAcol[j])*Bj[posBrow[i]];
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, posArow[irow[i]], posAcol[j])*Bj[posBrow[irow[i]]];
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, posArow[i], posAcol[icol[j]])*Bj[posBrow[i]];
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, posArow[irow[i]], posAcol[icol[j]])*Bj[posBrow[irow[i]]];
          }
        }
      break;
    }
  }
}

//====================================================================

void make_kronecker_tri_full2(int case_set, long long nrow, long long ncol,
                              int nrowA, double *A, int nrowB, double *B, int diag,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out)
{
  long long i, j;
  double *Bj;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA, A, posArow[i], posAcol[j])*Bj[posBrow[i]];
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA, A, posArow[irow[i]], posAcol[j])*Bj[posBrow[irow[i]]];
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA, A, posArow[i], posAcol[icol[j]])*Bj[posBrow[i]];
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA, A, posArow[irow[i]], posAcol[icol[j]])*Bj[posBrow[irow[i]]];
          }
        }
      break;
    }

  }else{ // No diagonal
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA, A, posArow[i], posAcol[j])*Bj[posBrow[i]];
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[j];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA, A, posArow[irow[i]], posAcol[j])*Bj[posBrow[irow[i]]];
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA, A, posArow[i], posAcol[icol[j]])*Bj[posBrow[i]];
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          Bj = B + nrowB*(long long)posBcol[icol[j]];
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA, A, posArow[irow[i]], posAcol[icol[j]])*Bj[posBrow[irow[i]]];
          }
        }
      break;
    }
  }
}

//====================================================================

void make_kronecker_tri_tri1_totri(long long n,
                                   int nA, double *A, int nB, double *B,
                                   int diag, int *posA, int *posB,
                                   int nindex, int *index, double *out)
{
  long long i, j, tmp;

  tmp = 0;
  if(nindex == 0){ // No subsetting
    if(diag){
      for(j=0; j<n; j++){
        for(i=0; i<j+1; i++){
          out[tmp++] = ij_tri_diag1(A,posA[i],posA[j])*ij_tri_diag1(B,posB[i],posB[j]);
        }
      }
    }else{
      for(j=0; j<n; j++){
        for(i=0; i<j; i++){
          out[tmp++] = ij_tri_nodiag1(A,posA[i],posA[j])*ij_tri_nodiag1(B,posB[i],posB[j]);
        }
      }
    }

  }else{
    if(diag){
      for(j=0; j<n; j++){
        for(i=0; i<j+1; i++){
          out[tmp++] = ij_tri_diag1(A,posA[index[i]],posA[index[j]])*ij_tri_diag1(B,posB[index[i]],posB[index[j]]);
        }
      }
    }else{
      for(j=0; j<n; j++){
        for(i=0; i<j; i++){
          out[tmp++] = ij_tri_nodiag1(A,posA[index[i]],posA[index[j]])*ij_tri_nodiag1(B,posB[index[i]],posB[index[j]]);
        }
      }
    }
  }
}

void make_kronecker_tri_tri1(int case_set, long long nrow, long long ncol,
                             int nrowA, double *A, int nrowB, double *B, int diag,
                             int *posArow, int *posBrow, int *posAcol, int *posBcol,
                             int *irow, int *icol, double *out)
{
  long long i, j;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A,posArow[i],posAcol[j])*ij_tri_diag1(B,posBrow[i],posBcol[j]);
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A,posArow[irow[i]],posAcol[j])*ij_tri_diag1(B,posBrow[irow[i]],posBcol[j]);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A,posArow[i],posAcol[icol[j]])*ij_tri_diag1(B,posBrow[i],posBcol[icol[j]]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A,posArow[irow[i]],posAcol[icol[j]])*ij_tri_diag1(B,posBrow[irow[i]],posBcol[icol[j]]);
          }
        }
      break;
    }

  }else{ // No diagonal
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A,posArow[i],posAcol[j])*ij_tri_nodiag1(B,posBrow[i],posBcol[j]);
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A,posArow[irow[i]],posAcol[j])*ij_tri_nodiag1(B,posBrow[irow[i]],posBcol[j]);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A,posArow[i],posAcol[icol[j]])*ij_tri_nodiag1(B,posBrow[i],posBcol[icol[j]]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A,posArow[irow[i]],posAcol[icol[j]])*ij_tri_nodiag1(B,posBrow[irow[i]],posBcol[icol[j]]);
          }
        }
      break;
    }
  }
}

//====================================================================

void make_kronecker_tri_tri2_totri(long long n,
                                   int nA, double *A, int nB, double *B,
                                   int diag, int *posA, int *posB,
                                   int nindex, int *index, double *out)
{
  long long i, j, tmp;

  tmp = 0;
  if(nindex == 0){ // No subsetting
    if(diag){
      for(j=0; j<n; j++){
        for(i=j; i<n; i++){
          out[tmp++] = ij_tri_diag2(nA,A,posA[i],posA[j])*ij_tri_diag2(nB,B,posB[i],posB[j]);
        }
      }
    }else{
      for(j=0; j<n; j++){
        for(i=j+1; i<n; i++){
          out[tmp++] = ij_tri_nodiag2(nA,A,posA[i],posA[j])*ij_tri_nodiag2(nB,B,posB[i],posB[j]);
        }
      }
    }

  }else{
    if(diag){
      for(j=0; j<n; j++){
        for(i=j; i<n; i++){
          out[tmp++] = ij_tri_diag2(nA,A,posA[index[i]],posA[index[j]])*ij_tri_diag2(nB,B,posB[index[i]],posB[index[j]]);
        }
      }
    }else{
      for(j=0; j<n; j++){
        for(i=j+1; i<n; i++){
          out[tmp++] = ij_tri_nodiag2(nA,A,posA[index[i]],posA[index[j]])*ij_tri_nodiag2(nB,B,posB[index[i]],posB[index[j]]);
        }
      }
    }
  }
}

void make_kronecker_tri_tri2(int case_set, long long nrow, long long ncol,
                             int nrowA, double *A, int nrowB, double *B, int diag,
                             int *posArow, int *posBrow, int *posAcol, int *posBcol,
                             int *irow, int *icol, double *out)
{
  long long i, j;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA,A,posArow[i],posAcol[j])*ij_tri_diag2(nrowB,B,posBrow[i],posBcol[j]);
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA, A,posArow[irow[i]],posAcol[j])*ij_tri_diag2(nrowB, B,posBrow[irow[i]],posBcol[j]);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA,A,posArow[i],posAcol[icol[j]])*ij_tri_diag2(nrowB,B,posBrow[i],posBcol[icol[j]]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nrowA,A,posArow[irow[i]],posAcol[icol[j]])*ij_tri_diag2(nrowB,B,posBrow[irow[i]],posBcol[icol[j]]);
          }
        }
      break;
    }

  }else{ // No diagonal
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA,A,posArow[i],posAcol[j])*ij_tri_nodiag2(nrowB,B,posBrow[i],posBcol[j]);
          }
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA,A,posArow[irow[i]],posAcol[j])*ij_tri_nodiag2(nrowB,B,posBrow[irow[i]],posBcol[j]);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA,A,posArow[i],posAcol[icol[j]])*ij_tri_nodiag2(nrowB,B,posBrow[i],posBcol[icol[j]]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nrowA,A,posArow[irow[i]],posAcol[icol[j]])*ij_tri_nodiag2(nrowB,B,posBrow[irow[i]],posBcol[icol[j]]);
          }
        }
      break;
    }
  }
}

//====================================================================
//====================================================================

void subset_tri1(int n, int nA, double *A, int diag,
                 int nindex, int *index, double *out)
{
  long long tmp, i, j;

  tmp = 0;
  if(nindex == 0)
  {
    if(diag){
      for(j=0; j<n; j++){
        for(i=0; i<j+1; i++){
          out[tmp++] = ij_tri_diag1(A, i, j);
        }
      }
    }else{ // without diagonal
      for(j=1; j<n; j++){
        for(i=0; i<j; i++){
          out[tmp++] = ij_tri_nodiag1(A, i, j);
        }
      }
    }

  }else{
    if(diag){
      for(j=0; j<n; j++){
        for(i=0; i<j+1; i++){
          out[tmp++] = ij_tri_diag1(A, index[i], index[j]);
        }
      }
    }else{ // without diagonal
      for(j=1; j<n; j++){
        for(i=0; i<j; i++){
          out[tmp++] = ij_tri_nodiag1(A, index[i], index[j]);
        }
      }
    }
  }
}

//====================================================================

void subset_tri2(int n, int nA, double *A, int diag,
                 int nindex, int *index, double *out)
{
  long long tmp, i, j;

  tmp = 0;
  if(nindex == 0)
  {
    if(diag){
      for(j=0; j<n; j++){
        for(i=j; i<n; i++){
          out[tmp++] = ij_tri_diag2(nA, A, i, j);
        }
      }

    }else{ // without diagonal
      for(j=0; j<n-1; j++){
        for(i=j+1; i<n; i++){
          out[tmp++] = ij_tri_nodiag2(nA, A, i, j);
        }
      }
    }

  }else{
    if(diag){
      for(j=0; j<n; j++){
        for(i=j; i<n; i++){
          out[tmp++] = ij_tri_diag2(nA, A, index[i], index[j]);
        }
      }
    }else{ // without diagonal
      for(j=0; j<n-1; j++){
        for(i=j+1; i<n; i++){
          out[tmp++] = ij_tri_nodiag2(nA, A, index[i], index[j]);
        }
      }
    }
  }
}

//====================================================================
//====================================================================

void unpack_tri1(int case_set, long long nrow, long long ncol,
                 int nA, double *A, int diag,
                 int *irow, int *icol, double *out)
{
  long long i, j;
  int tmp;
  int inc1=1;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<nA; j++){
          out[nA*j + j] = A[(j+1)*(j+2)/2 -1];
        }

        for(j=1; j<nA; j++){   // Start from 2nd column/row to n
          tmp = j;
          F77_NAME(dcopy)(&tmp, A + j*(j+1)/2, &inc1, out + nA*j, &inc1);
          F77_NAME(dcopy)(&tmp, A + j*(j+1)/2, &inc1, out + j, &nA);
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, irow[i], j);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, i, icol[j]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag1(A, irow[i], icol[j]);
          }
        }
      break;
    }

  }else{ // without diagonal but is set to 1
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<nA; j++){
          out[nA*j + j] = 1;
        }

        for(j=1; j<nA; j++){   // Start from 2nd column/row to n
          tmp = j;
          F77_NAME(dcopy)(&tmp, A + j*(j-1)/2, &inc1, out + nA*j, &inc1);
          F77_NAME(dcopy)(&tmp, A + j*(j-1)/2, &inc1, out + j, &nA);
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, irow[i], j);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, i, icol[j]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag1(A, irow[i], icol[j]);
          }
        }
      break;
    }
  }
}

//====================================================================

void unpack_tri2(int case_set, long long nrow, long long ncol,
                 int nA, double *A, int diag,
                 int *irow, int *icol, double *out)
{
  long long i, j;
  int tmp;
  int inc1=1;

  if(diag){
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<nA; j++){
          out[nA*j + j] = A[nA*j - j*(j-1)/2];
        }

        for(j=0; j<nA-1; j++){   // Start from 1st column/row to n-1
          tmp = nA-j-1;
          F77_NAME(dcopy)(&tmp, A + nA*j - j*(j-1)/2 + 1, &inc1, out + nA*(j+1) + j, &nA);
          F77_NAME(dcopy)(&tmp, A + nA*j - j*(j-1)/2 + 1, &inc1, out + nA*j + j + 1, &inc1);
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nA, A, irow[i], j);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nA, A, i, icol[j]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_diag2(nA, A, irow[i], icol[j]);
          }
        }
      break;
    }

  }else{ // without diagonal but is set to 1
    switch(case_set){
      case 0: // nirow=0 && nicol=0
        for(j=0; j<nA; j++){
          out[nA*j + j] = 1;
        }

        for(j=0; j<nA-1; j++){   // Start from 1st column/row to n-1
          tmp = nA-j-1;
          F77_NAME(dcopy)(&tmp, A + nA*j - j*(j+1)/2, &inc1, out + nA*(j+1) + j, &nA);
          F77_NAME(dcopy)(&tmp, A + nA*j - j*(j+1)/2, &inc1, out + nA*j + j + 1, &inc1);
        }
      break;

      case 1: // nirow>0 && nicol=0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nA, A, irow[i], j);
          }
        }
      break;

      case 2: // nirow=0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nA, A, i, icol[j]);
          }
        }
      break;

      case 3: // nirow>0 && nicol>0
        for(j=0; j<ncol; j++){
          for(i=0; i<nrow; i++){
            out[nrow*j + i] = ij_tri_nodiag2(nA, A, irow[i], icol[j]);
          }
        }
      break;
    }
  }
}
