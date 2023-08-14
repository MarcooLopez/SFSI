

double soft_threshold(double z, double penalty);

int imax_integer(int n, int *x);

int all_equal(long long n1, int *set1, long long n2, int *set2);

void matrix_vector_product(int nrow,
                           int ncol,
                           double *a,
                           double *A,
                           double *x,
                           int incx,
                           double *y,
                           int transpose);

//void matrix_matrix_product(int nrowA, int ncolA, int ncolB, double *A,
//                           double *B, double *C)

void crossproduct(int nrow, int ncolA, int ncolB,
                 double *A, double *B, double *C);

void tcrossproduct(int nrowA, int ncol, int nrowB,
                  double *A, double *B, double *C);

void tcrossproduct_tri(int n, double *A, double *B, double *C);

void crossproduct_scale(int nrow, int ncolA, int ncolB,
                  double *A, double *d, double *B,
                  double *C, double *work);

void tcrossproduct_scale(int nrowA, int ncol, int nrowB,
                   double *A, double *d, double *B,
                   double *C, double *work);

void slice_matrix(int nrow, double *A, double *x,
                 int n, int *index,
                 int k,
                 int margin);

void resize_matrix(int nrow, int ncol, int nrownew, int ncolnew, double *A);

void backsolvet(int n, double *A, double *b);

void backsolve(int n, double *A, double *b);

void update_chol(int n, double *A, int nR, double *R,
                 int k, int *index, double *eps,
                 double *work, int *info);

void invert_matrix(int n, double *A, double *Ainv, double *eps, double *work);

void matrix_vector_product_subset(int nrow,
                                 int ncol,
                                 double *A,
                                 double *x,
                                 double *y,
                                 int nirow, int *irow,
                                 int nicol, int *icol,
                                 int transpose,
                                 double *work);

void append_to_sorted_vector_integer(int n,
                                    int *v,
                                    int k,
                                    int *values);

void reduce_vector_integer(int n,
                           int *v,
                           int k,
                           int *index);

double ij_tri_diag1(double *A, long long i, long long j);

double ij_tri_nodiag1(double *A, long long i, long long j);

double ij_tri_diag2(int nA, double *A, long long i, long long j);

double ij_tri_nodiag2(int nA, double *A, long long i, long long j);

void make_kronecker_full_full(int case_set, long long nrow, long long ncol,
                              int nrowA, int ncolA, double *A,
                              int nrowB, int ncolB, double *B,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out);

void make_kronecker_tri_full1(int case_set, long long nrow, long long ncol,
                              int nrowA, double *A, int nrowB, double *B, int diag,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out);

void make_kronecker_tri_full2(int case_set, long long nrow, long long ncol,
                              int nrowA, double *A, int nrowB, double *B, int diag,
                              int *posArow, int *posBrow, int *posAcol, int *posBcol,
                              int *irow, int *icol, double *out);

void make_kronecker_tri_tri1_totri(long long n,
                                   int nA, double *A, int nB, double *B,
                                   int diag, int *posA, int *posB,
                                   int nindex, int *index, double *out);

void make_kronecker_tri_tri1(int case_set, long long nrow, long long ncol,
                             int nrowA, double *A, int nrowB, double *B, int diag,
                             int *posArow, int *posBrow, int *posAcol, int *posBcol,
                             int *irow, int *icol, double *out);

void make_kronecker_tri_tri2_totri(long long n,
                                   int nA, double *A, int nB, double *B,
                                   int diag, int *posA, int *posB,
                                   int nindex, int *index, double *out);

void make_kronecker_tri_tri2(int case_set, long long nrow, long long ncol,
                             int nrowA, double *A, int nrowB, double *B, int diag,
                             int *posArow, int *posBrow, int *posAcol, int *posBcol,
                             int *irow, int *icol, double *out);

void subset_tri1(int n, int nA, double *A, int diag,
                int nindex, int *index, double *out);

void subset_tri2(int n, int nA, double *A, int diag,
                int nindex, int *index, double *out);

void unpack_tri1(int case_set, long long nrow, long long ncol,
                int nA, double *A, int diag,
                int *irow, int *icol, double *out);

void unpack_tri2(int case_set, long long nrow, long long ncol,
                int nA, double *A, int diag,
                int *irow, int *icol, double *out);
