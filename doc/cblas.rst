.. index::
   single: Low-level CBLAS
   single: CBLAS, Low-level interface
   single: BLAS, Low-level C interface
   single: Basic Linear Algebra Subroutines (BLAS)

.. _chap_cblas:

*****************
GSL CBLAS Library
*****************

.. include:: include.rst

The prototypes for the low-level |cblas| functions are declared in
the file :file:`gsl_cblas.h`.  For the definition of the functions
consult the documentation available from Netlib (:ref:`see BLAS References and
Further Reading <sec_blas-references>`).

Level 1
=======

.. function:: float cblas_sdsdot (const int N, const float alpha, const float * x, const int incx, const float * y, const int incy)

.. function:: double cblas_dsdot (const int N, const float * x, const int incx, const float * y, const int incy)

.. function:: float cblas_sdot (const int N, const float * x, const int incx, const float * y, const int incy)

.. function:: double cblas_ddot (const int N, const double * x, const int incx, const double * y, const int incy)

.. function:: void cblas_cdotu_sub (const int N, const void * x, const int incx, const void * y, const int incy, void * dotu)

.. function:: void cblas_cdotc_sub (const int N, const void * x, const int incx, const void * y, const int incy, void * dotc)

.. function:: void cblas_zdotu_sub (const int N, const void * x, const int incx, const void * y, const int incy, void * dotu)

.. function:: void cblas_zdotc_sub (const int N, const void * x, const int incx, const void * y, const int incy, void * dotc)

.. function:: float cblas_snrm2 (const int N, const float * x, const int incx)

.. function:: float cblas_sasum (const int N, const float * x, const int incx)

.. function:: double cblas_dnrm2 (const int N, const double * x, const int incx)

.. function:: double cblas_dasum (const int N, const double * x, const int incx)

.. function:: float cblas_scnrm2 (const int N, const void * x, const int incx)

.. function:: float cblas_scasum (const int N, const void * x, const int incx)

.. function:: double cblas_dznrm2 (const int N, const void * x, const int incx)

.. function:: double cblas_dzasum (const int N, const void * x, const int incx)

.. function:: CBLAS_INDEX cblas_isamax (const int N, const float * x, const int incx)

.. function:: CBLAS_INDEX cblas_idamax (const int N, const double * x, const int incx)

.. function:: CBLAS_INDEX cblas_icamax (const int N, const void * x, const int incx)

.. function:: CBLAS_INDEX cblas_izamax (const int N, const void * x, const int incx)

.. function:: void cblas_sswap (const int N, float * x, const int incx, float * y, const int incy)

.. function:: void cblas_scopy (const int N, const float * x, const int incx, float * y, const int incy)

.. function:: void cblas_saxpy (const int N, const float alpha, const float * x, const int incx, float * y, const int incy)

.. function:: void cblas_dswap (const int N, double * x, const int incx, double * y, const int incy)

.. function:: void cblas_dcopy (const int N, const double * x, const int incx, double * y, const int incy)

.. function:: void cblas_daxpy (const int N, const double alpha, const double * x, const int incx, double * y, const int incy)

.. function:: void cblas_cswap (const int N, void * x, const int incx, void * y, const int incy)

.. function:: void cblas_ccopy (const int N, const void * x, const int incx, void * y, const int incy)

.. function:: void cblas_caxpy (const int N, const void * alpha, const void * x, const int incx, void * y, const int incy)

.. function:: void cblas_zswap (const int N, void * x, const int incx, void * y, const int incy)

.. function:: void cblas_zcopy (const int N, const void * x, const int incx, void * y, const int incy)

.. function:: void cblas_zaxpy (const int N, const void * alpha, const void * x, const int incx, void * y, const int incy)

.. function:: void cblas_srotg (float * a, float * b, float * c, float * s)

.. function:: void cblas_srotmg (float * d1, float * d2, float * b1, const float b2, float * P)

.. function:: void cblas_srot (const int N, float * x, const int incx, float * y, const int incy, const float c, const float s)

.. function:: void cblas_srotm (const int N, float * x, const int incx, float * y, const int incy, const float * P)

.. function:: void cblas_drotg (double * a, double * b, double * c, double * s)

.. function:: void cblas_drotmg (double * d1, double * d2, double * b1, const double b2, double * P)

.. function:: void cblas_drot (const int N, double * x, const int incx, double * y, const int incy, const double c, const double s)

.. function:: void cblas_drotm (const int N, double * x, const int incx, double * y, const int incy, const double * P)

.. function:: void cblas_sscal (const int N, const float alpha, float * x, const int incx)

.. function:: void cblas_dscal (const int N, const double alpha, double * x, const int incx)

.. function:: void cblas_cscal (const int N, const void * alpha, void * x, const int incx)

.. function:: void cblas_zscal (const int N, const void * alpha, void * x, const int incx)

.. function:: void cblas_csscal (const int N, const float alpha, void * x, const int incx)

.. function:: void cblas_zdscal (const int N, const double alpha, void * x, const int incx)

Level 2
=======

.. function:: void cblas_sgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const float alpha, const float * A, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)

.. function:: void cblas_sgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const int KL, const int KU, const float alpha, const float * A, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)

.. function:: void cblas_strmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const float * A, const int lda, float * x, const int incx)

.. function:: void cblas_stbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const float * A, const int lda, float * x, const int incx)

.. function:: void cblas_stpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const float * Ap, float * x, const int incx)

.. function:: void cblas_strsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const float * A, const int lda, float * x, const int incx)

.. function:: void cblas_stbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const float * A, const int lda, float * x, const int incx)

.. function:: void cblas_stpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const float * Ap, float * x, const int incx)

.. function:: void cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double * A, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)

.. function:: void cblas_dgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const int KL, const int KU, const double alpha, const double * A, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)

.. function:: void cblas_dtrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const double * A, const int lda, double * x, const int incx)

.. function:: void cblas_dtbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const double * A, const int lda, double * x, const int incx)

.. function:: void cblas_dtpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const double * Ap, double * x, const int incx)

.. function:: void cblas_dtrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const double * A, const int lda, double * x, const int incx)

.. function:: void cblas_dtbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const double * A, const int lda, double * x, const int incx)

.. function:: void cblas_dtpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const double * Ap, double * x, const int incx)

.. function:: void cblas_cgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_cgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const int KL, const int KU, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_ctrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ctbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ctpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * Ap, void * x, const int incx)

.. function:: void cblas_ctrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ctbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ctpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * Ap, void * x, const int incx)

.. function:: void cblas_zgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_zgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const int M, const int N, const int KL, const int KU, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_ztrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ztbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ztpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * Ap, void * x, const int incx)

.. function:: void cblas_ztrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ztbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const int K, const void * A, const int lda, void * x, const int incx)

.. function:: void cblas_ztpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int N, const void * Ap, void * x, const int incx)

.. function:: void cblas_ssymv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * A, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)

.. function:: void cblas_ssbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const int K, const float alpha, const float * A, const int lda, const float * x, const int incx, const float beta, float * y, const int incy)

.. function:: void cblas_sspmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * Ap, const float * x, const int incx, const float beta, float * y, const int incy)

.. function:: void cblas_sger (const enum CBLAS_ORDER order, const int M, const int N, const float alpha, const float * x, const int incx, const float * y, const int incy, float * A, const int lda)

.. function:: void cblas_ssyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * x, const int incx, float * A, const int lda)

.. function:: void cblas_sspr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * x, const int incx, float * Ap)

.. function:: void cblas_ssyr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * x, const int incx, const float * y, const int incy, float * A, const int lda)

.. function:: void cblas_sspr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const float * x, const int incx, const float * y, const int incy, float * A)

.. function:: void cblas_dsymv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * A, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)

.. function:: void cblas_dsbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const int K, const double alpha, const double * A, const int lda, const double * x, const int incx, const double beta, double * y, const int incy)

.. function:: void cblas_dspmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * Ap, const double * x, const int incx, const double beta, double * y, const int incy)

.. function:: void cblas_dger (const enum CBLAS_ORDER order, const int M, const int N, const double alpha, const double * x, const int incx, const double * y, const int incy, double * A, const int lda)

.. function:: void cblas_dsyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * x, const int incx, double * A, const int lda)

.. function:: void cblas_dspr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * x, const int incx, double * Ap)

.. function:: void cblas_dsyr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * x, const int incx, const double * y, const int incy, double * A, const int lda)

.. function:: void cblas_dspr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double * x, const int incx, const double * y, const int incy, double * A)

.. function:: void cblas_chemv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_chbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const int K, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_chpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * Ap, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_cgeru (const enum CBLAS_ORDER order, const int M, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_cgerc (const enum CBLAS_ORDER order, const int M, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_cher (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const void * x, const int incx, void * A, const int lda)

.. function:: void cblas_chpr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const float alpha, const void * x, const int incx, void * A)

.. function:: void cblas_cher2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_chpr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * Ap)

.. function:: void cblas_zhemv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_zhbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const int K, const void * alpha, const void * A, const int lda, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_zhpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * Ap, const void * x, const int incx, const void * beta, void * y, const int incy)

.. function:: void cblas_zgeru (const enum CBLAS_ORDER order, const int M, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_zgerc (const enum CBLAS_ORDER order, const int M, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_zher (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const void * x, const int incx, void * A, const int lda)

.. function:: void cblas_zhpr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const void * x, const int incx, void * A)

.. function:: void cblas_zher2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * A, const int lda)

.. function:: void cblas_zhpr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const void * alpha, const void * x, const int incx, const void * y, const int incy, void * Ap)

Level 3
=======

.. function:: void cblas_sgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, const float * A, const int lda, const float * B, const int ldb, const float beta, float * C, const int ldc)

.. function:: void cblas_ssymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const float alpha, const float * A, const int lda, const float * B, const int ldb, const float beta, float * C, const int ldc)

.. function:: void cblas_ssyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const float alpha, const float * A, const int lda, const float beta, float * C, const int ldc)

.. function:: void cblas_ssyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const float alpha, const float * A, const int lda, const float * B, const int ldb, const float beta, float * C, const int ldc)

.. function:: void cblas_strmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float * A, const int lda, float * B, const int ldb)

.. function:: void cblas_strsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const float alpha, const float * A, const int lda, float * B, const int ldb)

.. function:: void cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double * A, const int lda, const double * B, const int ldb, const double beta, double * C, const int ldc)

.. function:: void cblas_dsymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const double alpha, const double * A, const int lda, const double * B, const int ldb, const double beta, double * C, const int ldc)

.. function:: void cblas_dsyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const double alpha, const double * A, const int lda, const double beta, double * C, const int ldc)

.. function:: void cblas_dsyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const double alpha, const double * A, const int lda, const double * B, const int ldb, const double beta, double * C, const int ldc)

.. function:: void cblas_dtrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double * A, const int lda, double * B, const int ldb)

.. function:: void cblas_dtrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double * A, const int lda, double * B, const int ldb)

.. function:: void cblas_cgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_csymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_csyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * beta, void * C, const int ldc)

.. function:: void cblas_csyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_ctrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const void * alpha, const void * A, const int lda, void * B, const int ldb)

.. function:: void cblas_ctrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const void * alpha, const void * A, const int lda, void * B, const int ldb)

.. function:: void cblas_zgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_zsymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_zsyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * beta, void * C, const int ldc)

.. function:: void cblas_zsyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_ztrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const void * alpha, const void * A, const int lda, void * B, const int ldb)

.. function:: void cblas_ztrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const int M, const int N, const void * alpha, const void * A, const int lda, void * B, const int ldb)

.. function:: void cblas_chemm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_cherk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const float alpha, const void * A, const int lda, const float beta, void * C, const int ldc)

.. function:: void cblas_cher2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const float beta, void * C, const int ldc)

.. function:: void cblas_zhemm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const void * beta, void * C, const int ldc)

.. function:: void cblas_zherk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const double alpha, const void * A, const int lda, const double beta, void * C, const int ldc)

.. function:: void cblas_zher2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const int N, const int K, const void * alpha, const void * A, const int lda, const void * B, const int ldb, const double beta, void * C, const int ldc)

.. function:: void cblas_xerbla (int p, const char * rout, const char * form, ...)

Examples
========

The following program computes the product of two matrices using the
Level-3 |blas| function SGEMM,

.. only:: not texinfo

   .. math::

      \left(
        \begin{array}{ccc}
          0.11 & 0.12 & 0.13\\
          0.21 & 0.22 & 0.23
        \end{array}
      \right)
      \left(
        \begin{array}{cc}
          1011 & 1012\\
          1021 & 1022\\
          1031 & 1032
        \end{array}
      \right)
      =
      \left(
        \begin{array}{cc}
          367.76 & 368.12\\
          674.06 & 674.72
        \end{array}
      \right)

.. only:: texinfo

   ::

      [ 0.11 0.12 0.13 ]  [ 1011 1012 ]     [ 367.76 368.12 ]
      [ 0.21 0.22 0.23 ]  [ 1021 1022 ]  =  [ 674.06 674.72 ]
                          [ 1031 1032 ]

The matrices are stored in row major order but could be stored in column
major order if the first argument of the call to :func:`cblas_sgemm` was
changed to :code:`CblasColMajor`.

.. include:: examples/cblas.c
   :code:

To compile the program use the following command line::

  $ gcc -Wall demo.c -lgslcblas

There is no need to link with the main library :code:`-lgsl` in this
case as the |cblas| library is an independent unit. Here is the output
from the program,

.. include:: examples/cblas.txt
   :code:
