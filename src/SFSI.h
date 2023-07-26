

double soft_threshold(double z, double penalty);

int imax_integer(int n, int *x);

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
