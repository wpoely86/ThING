#ifndef LAPACK_H
#define LAPACK_H

extern "C" {

   void dcopy_(int *n,double *x,int *incx,double *y,int *incy);
   void dscal_(int *n,double *alpha,double *x,int *incx);
   void dgemm_(char *transA,char *transB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
   double dlansy_(char *norm, char *uplo, int *N, double *A, int *lda, double *work);
   void dsygvd_(int *itype, char *jobz, char *uplo, int *N, double *A, int *lda, double *B, int *ldb, double *W, double *work, int *lwork, int *iwork, int *liwork, int *info);
   void dsytrf_(char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info);
   void dsycon_(char *uplo, int *n, double *A, int *lda, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);

}

#endif
