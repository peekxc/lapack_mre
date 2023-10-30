#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

//definition for lapack / blas routines
#ifdef __cplusplus
extern "C" {
#endif
  int dgemm_( char *transa, char *transb, const int *m, const int *n, const int *k,
            double *alpha, const double *a, const int *lda, const double *b,
            const int *ldb, double *beta, double *c, const int *ldc);
  int dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
#ifdef __cplusplus
}
#endif


void mk_diag(double* diag, double* mat, int n, const int m) {
    int i;
    for(i =0; i < m; ++i)
        mat[i*n+i] = diag[i];
}

int main() {
    int i, j;
    double m[16] = { //rowmajor
        0, 1, 2, 3,
        1, 2, 6, 4,
        2, 6, 7, 5,
        3, 4, 5, 6 };
    int n = 4;
    int info = 0;
    char jobz = 'V';
    char uplo = 'U';
    int lwork = n * n;

//print initial matrix
    printf ("\nInitial Matrix (A):\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j){
            //printf("[%5.1f] ", m[i][j] );
            printf("[%6.1f] ", m[i*n+j]);
        }
        printf("\n\t");
    }
//find eigenvalues and eigenvectors of a symmetric matrix
    double* x = (double*) malloc( n*sizeof(double));
    double* u = (double*) malloc( n*n*sizeof(double));
    double* work = (double*) malloc( n*n*sizeof(double));
    memcpy (u, m, n*n*sizeof(double));

    dsyev_(&jobz,   //include eigenvectors? (V=yes)
           &uplo,   //Upper triangle of A is stored
           &n,      //order of the matrix
           u,       //** input matrix, returns orthonormal eigenvectors
           &n,      //leading dimension of array A
           x,       //** returns eigenvalues (as ROWS)
           work,    //work space
           &lwork,  //length of work; LWORK >= max(1,3*n-1)
           &info);  //**returns info on results
                    //== 0, success; 
                    // < 0, -ith argument had issue;
                    // > 0, failed in convergence
    if ( info != 0 ){
        printf("error: %d", info);
        return 0;
    }
//create all the vectors/matrices we need
//  diagonal matrix of eigenvalues...
    double* d = (double*) malloc( n*n*sizeof(double));
    mk_diag( x, d, n, n);

    printf("\nD :\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", d[i*n+j]);
        printf("\n\t"); 
    }

    printf("\nUt :\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", u[i*n+j]);
        printf("\n\t"); 
    }
    char _tran = 'T', ntran = 'N';
    double one = 1.0;
    double* I = (double*) malloc( n*n*sizeof(double));
    dgemm_(   &ntran, //const char BLAS_TRANSPOSE  ~(op) on A (Cblas[Trans/NoTrans/Conj/ConjTrans])
              &_tran, //const char BLAS_TRANSPOSE  ~(op) on B (Cblas[Trans/NoTrans/Conj/ConjTrans])
               &n,    //const int  M                ~number of rows in A
               &n,    //const int  N                ~number of cols in B 
               &n,    //const int  K                ~number of cols in A
               &one,  //const double alpha          ~scalar to multiply after computation
                u,    //const double *A             ~matrix A
               &n,    //const int lda               ~stride of A (inc for inner loop)
                u,    //const double *B             ~matrix B
               &n,    //const int ldb               ~stride of B (inc for inner loop)
               &one,  //const double beta           ~scalar to multiply C by
                I,    //      double *C             ~matrix C -- and output
               &n );  //const int ldc               ~stride of C (inc for inner loop)

    printf("\nU * Ut == I :\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", I[i*n+j]);
        printf("\n\t"); 
    }

//---------------------------------------------------------------------------------------
    double* ud = (double*) malloc( n*n*sizeof(double));
    dgemm_(   &ntran, //const char BLAS_TRANSPOSE  ~(op) on A (Cblas[Trans/NoTrans/Conj/ConjTrans])
              &ntran, //const char BLAS_TRANSPOSE  ~(op) on B (Cblas[Trans/NoTrans/Conj/ConjTrans])
               &n,    //const int  M                ~number of rows in A
               &n,    //const int  N                ~number of cols in B
               &n,    //const int  K                ~number of cols in A
               &one,  //const double alpha          ~scalar to multiply after computation
                u,    //const double *A             ~matrix A
               &n,    //const int lda               ~stride of A (inc for inner loop)
                d,    //const double *B             ~matrix B
               &n,    //const int ldb               ~stride of B (inc for inner loop)
               &one,  //const double beta           ~scalar to multiply C by
                ud,   //      double *C             ~matrix C -- and output
               &n );  //const int ldc               ~stride of C (inc for inner loop)

    printf("\n U * D :\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", ud[i*n+j]);
        printf("\n\t"); 
    }

    double* udut = (double*) malloc( n*n*sizeof(double));
    //PERFORMS: C = b*(C)+a*(A*B)
    dgemm_(   &ntran, //const char BLAS_TRANSPOSE  ~(op) on A (Cblas[Trans/NoTrans/Conj/ConjTrans])
              &_tran, //const char BLAS_TRANSPOSE  ~(op) on B (Cblas[Trans/NoTrans/Conj/ConjTrans])
               &n,    //const int  M                ~number of rows in A
               &n,    //const int  N                ~number of cols in B
               &n,    //const int  K                ~number of cols in A
               &one,  //const double alpha          ~scalar to multiply after computation
                ud,    //const double *A             ~matrix A
               &n,    //const int lda               ~stride of A (inc for inner loop)
                u,    //const double *B             ~matrix B
               &n,    //const int ldb               ~stride of B (inc for inner loop)
               &one,  //const double beta           ~scalar to multiply C by
               udut,  //      double *C             ~matrix C -- and output
               &n );  //const int ldc               ~stride of C (inc for inner loop)

    printf("\n(U * D) * U^t == A :\n\t");
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j)
            printf("[%6.5f] ", udut[i*n+j]);
        printf("\n\t"); 
    }
    putchar ('\n');
    return 0;
}