.. index::
   single: linear algebra, BLAS
   single: matrix, operations
   single: vector, operations
   single: BLAS
   single: CBLAS
   single: Basic Linear Algebra Subroutines (BLAS)

.. _chap_blas-support:

************
BLAS Support
************

The Basic Linear Algebra Subprograms (BLAS) define a set of fundamental
operations on vectors and matrices which can be used to create optimized
higher-level linear algebra functionality.

The library provides a low-level layer which corresponds directly to the
C-language BLAS standard, referred to here as "CBLAS", and a
higher-level interface for operations on GSL vectors and matrices.
Users who are interested in simple operations on GSL vector and matrix
objects should use the high-level layer described
in this chapter.  The functions are declared in the file
:file:`gsl_blas.h` and should satisfy the needs of most users.  

Note that GSL matrices are implemented using dense-storage so the
interface only includes the corresponding dense-storage BLAS
functions.  The full BLAS functionality for band-format and
packed-format matrices is available through the low-level CBLAS
interface.  Similarly, GSL vectors are restricted to positive strides,
whereas the low-level CBLAS interface supports negative
strides as specified in the BLAS standard [#f1]_.
 
The interface for the :code:`gsl_cblas` layer is specified in the file
:file:`gsl_cblas.h`.  This interface corresponds to the BLAS Technical
Forum's standard for the C interface to legacy BLAS
implementations. Users who have access to other conforming CBLAS
implementations can use these in place of the version provided by the
library.  Note that users who have only a Fortran BLAS library can
use a CBLAS conformant wrapper to convert it into a CBLAS
library.  A reference CBLAS wrapper for legacy Fortran
implementations exists as part of the CBLAS standard and can
be obtained from Netlib.  The complete set of CBLAS functions is
listed in an :ref:`appendix <chap_cblas>`.

There are three levels of BLAS operations,

=========== ===============================================================
**Level 1** Vector operations, e.g. :math:`y = \alpha x + y`
**Level 2** Matrix-vector operations, e.g. :math:`y = \alpha A x + \beta y`
**Level 3** Matrix-matrix operations, e.g. :math:`C = \alpha A B + C`
=========== ===============================================================

Each routine has a name which specifies the operation, the type of
matrices involved and their precisions.  Some of the most common
operations and their names are given below,

======== =====================================
**DOT**  scalar product, :math:`x^T y`
**AXPY** vector sum, :math:`\alpha x + y`
**MV**   matrix-vector product, :math:`A x`
**SV**   matrix-vector solve, :math:`inv(A) x`
**MM**   matrix-matrix product, :math:`A B`
**SM**   matrix-matrix solve, :math:`inv(A) B`
======== =====================================

The types of matrices are,

====== =================
**GE** general
**GB** general band
**SY** symmetric
**SB** symmetric band
**SP** symmetric packed
**HE** hermitian
**HB** hermitian band
**HP** hermitian packed
**TR** triangular 
**TB** triangular band
**TP** triangular packed
====== =================

Each operation is defined for four precisions,

===== ==============
**S** single real
**D** double real
**C** single complex
**Z** double complex
===== ==============

Thus, for example, the name SGEMM stands for "single-precision
general matrix-matrix multiply" and ZGEMM stands for
"double-precision complex matrix-matrix multiply".

Note that the vector and matrix arguments to BLAS functions must not
be aliased, as the results are undefined when the underlying arrays
overlap (:ref:`aliasing-of-arrays`).

GSL BLAS Interface
==================

GSL provides dense vector and matrix objects, based on the relevant
built-in types.  The library provides an interface to the BLAS
operations which apply to these objects.  The interface to this
functionality is given in the file :file:`gsl_blas.h`.

.. CblasNoTrans, CblasTrans, CblasConjTrans
.. CblasUpper, CblasLower
.. CblasNonUnit, CblasUnit
.. CblasLeft, CblasRight

Level 1
-------

.. index::
   single: DOT, Level-1 BLAS

.. function:: int gsl_blas_sdsdot (float alpha, const gsl_vector_float * x, const gsl_vector_float * y, float * result)

   This function computes the sum :math:`\alpha + x^T y` for the vectors
   :data:`x` and :data:`y`, returning the result in :data:`result`.

.. function:: int gsl_blas_sdot (const gsl_vector_float * x, const gsl_vector_float * y, float * result)
              int gsl_blas_dsdot (const gsl_vector_float * x, const gsl_vector_float * y, double * result)
              int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)

   These functions compute the scalar product :math:`x^T y` for the vectors
   :data:`x` and :data:`y`, returning the result in :data:`result`.

.. function:: int gsl_blas_cdotu (const gsl_vector_complex_float * x, const gsl_vector_complex_float * y, gsl_complex_float * dotu)
              int gsl_blas_zdotu (const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_complex * dotu)

   These functions compute the complex scalar product :math:`x^T y` for the
   vectors :data:`x` and :data:`y`, returning the result in :data:`dotu`

.. function:: int gsl_blas_cdotc (const gsl_vector_complex_float * x, const gsl_vector_complex_float * y, gsl_complex_float * dotc)
              int gsl_blas_zdotc (const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_complex * dotc)

   These functions compute the complex conjugate scalar product :math:`x^H y`
   for the vectors :data:`x` and :data:`y`, returning the result in
   :data:`dotc`

.. index::
   single: NRM2, Level-1 BLAS

.. function:: float gsl_blas_snrm2 (const gsl_vector_float * x)
              double gsl_blas_dnrm2 (const gsl_vector * x)

   These functions compute the Euclidean norm 
   :math:`||x||_2 = \sqrt{\sum x_i^2}` of the vector :data:`x`.

.. function:: float gsl_blas_scnrm2 (const gsl_vector_complex_float * x)
              double gsl_blas_dznrm2 (const gsl_vector_complex * x)

   These functions compute the Euclidean norm of the complex vector :data:`x`,

   .. math:: ||x||_2 = \sqrt{\sum (\Re(x_i)^2 + \Im(x_i)^2)}.

.. index::
   single: ASUM, Level-1 BLAS

.. function:: float gsl_blas_sasum (const gsl_vector_float * x)
              double gsl_blas_dasum (const gsl_vector * x)

   These functions compute the absolute sum :math:`\sum |x_i|` of the
   elements of the vector :data:`x`.

.. function:: float gsl_blas_scasum (const gsl_vector_complex_float * x)
              double gsl_blas_dzasum (const gsl_vector_complex * x)

   These functions compute the sum of the magnitudes of the real and
   imaginary parts of the complex vector :data:`x`, 
   :math:`\sum \left( |\Re(x_i)| + |\Im(x_i)| \right)`.

.. index::
   single: AMAX, Level-1 BLAS

.. function:: CBLAS_INDEX_t gsl_blas_isamax (const gsl_vector_float * x)
              CBLAS_INDEX_t gsl_blas_idamax (const gsl_vector * x)
              CBLAS_INDEX_t gsl_blas_icamax (const gsl_vector_complex_float * x)
              CBLAS_INDEX_t gsl_blas_izamax (const gsl_vector_complex * x)

   These functions return the index of the largest element of the vector
   :data:`x`. The largest element is determined by its absolute magnitude for
   real vectors and by the sum of the magnitudes of the real and imaginary
   parts :math:`|\Re(x_i)| + |\Im(x_i)|` for complex vectors.  If the
   largest value occurs several times then the index of the first
   occurrence is returned.

.. index::
   single: SWAP, Level-1 BLAS

.. function:: int gsl_blas_sswap (gsl_vector_float * x, gsl_vector_float * y)
              int gsl_blas_dswap (gsl_vector * x, gsl_vector * y)
              int gsl_blas_cswap (gsl_vector_complex_float * x, gsl_vector_complex_float * y)
              int gsl_blas_zswap (gsl_vector_complex * x, gsl_vector_complex * y)

   These functions exchange the elements of the vectors :data:`x` and :data:`y`.

.. index::
   single: COPY, Level-1 BLAS

.. function:: int gsl_blas_scopy (const gsl_vector_float * x, gsl_vector_float * y)
              int gsl_blas_dcopy (const gsl_vector * x, gsl_vector * y)
              int gsl_blas_ccopy (const gsl_vector_complex_float * x, gsl_vector_complex_float * y)
              int gsl_blas_zcopy (const gsl_vector_complex * x, gsl_vector_complex * y)

   These functions copy the elements of the vector :data:`x` into the vector
   :data:`y`.

.. index::
   single: AXPY, Level-1 BLAS
   single: DAXPY, Level-1 BLAS
   single: SAXPY, Level-1 BLAS

.. function:: int gsl_blas_saxpy (float alpha, const gsl_vector_float * x, gsl_vector_float * y)
              int gsl_blas_daxpy (double alpha, const gsl_vector * x, gsl_vector * y)
              int gsl_blas_caxpy (const gsl_complex_float alpha, const gsl_vector_complex_float * x, gsl_vector_complex_float * y)
              int gsl_blas_zaxpy (const gsl_complex alpha, const gsl_vector_complex * x, gsl_vector_complex * y)

   These functions compute the sum :math:`y = \alpha x + y` for the vectors
   :data:`x` and :data:`y`.

.. index::
   single: SCAL, Level-1 BLAS

.. function:: void gsl_blas_sscal (float alpha, gsl_vector_float * x)
              void gsl_blas_dscal (double alpha, gsl_vector * x)
              void gsl_blas_cscal (const gsl_complex_float alpha, gsl_vector_complex_float * x)
              void gsl_blas_zscal (const gsl_complex alpha, gsl_vector_complex * x)
              void gsl_blas_csscal (float alpha, gsl_vector_complex_float * x)
              void gsl_blas_zdscal (double alpha, gsl_vector_complex * x)

   These functions rescale the vector :data:`x` by the multiplicative factor
   :data:`alpha`.

.. index::
   single: ROTG, Level-1 BLAS
   single: Givens Rotation, BLAS

.. function:: int gsl_blas_srotg (float a[], float b[], float c[], float s[])
              int gsl_blas_drotg (double a[], double b[], double c[], double s[])

   These functions compute a Givens rotation :math:`(c,s)` which zeroes the
   vector :math:`(a,b)`,

   .. only:: not texinfo

      .. math::

         \left(
         \begin{matrix}
            c & s \\
           -s & c
         \end{matrix}
         \right)
         \left(
         \begin{matrix}
           a \\
           b
         \end{matrix}
         \right)
         =
         \left(
         \begin{matrix}
           r' \\
           0
         \end{matrix}
         \right)

   .. only:: texinfo

      ::

         [  c  s ] [ a ] = [ r ]
         [ -s  c ] [ b ]   [ 0 ]

   The variables :data:`a` and :data:`b` are overwritten by the routine.

.. function:: int gsl_blas_srot (gsl_vector_float * x, gsl_vector_float * y, float c, float s)
              int gsl_blas_drot (gsl_vector * x, gsl_vector * y, const double c, const double s)

   These functions apply a Givens rotation :math:`(x', y') = (c x + s y, -s x + c y)`
   to the vectors :data:`x`, :data:`y`.

.. index::
   single: Modified Givens Rotation, BLAS
   single: Givens Rotation, Modified, BLAS

.. function:: int gsl_blas_srotmg (float d1[], float d2[], float b1[], float b2, float P[])
              int gsl_blas_drotmg (double d1[], double d2[], double b1[], double b2, double P[])

   These functions compute a modified Givens transformation.  The modified
   Givens transformation is defined in the original Level-1 BLAS
   specification, given in the references.

.. function:: int gsl_blas_srotm (gsl_vector_float * x, gsl_vector_float * y, const float P[])
              int gsl_blas_drotm (gsl_vector * x, gsl_vector * y, const double P[])

   These functions apply a modified Givens transformation.  

Level 2
-------

.. index::
   single: GEMV, Level-2 BLAS

.. function:: int gsl_blas_sgemv (CBLAS_TRANSPOSE_t TransA, float alpha, const gsl_matrix_float * A, const gsl_vector_float * x, float beta, gsl_vector_float * y)
              int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
              int gsl_blas_cgemv (CBLAS_TRANSPOSE_t TransA, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_vector_complex_float * x, const gsl_complex_float beta, gsl_vector_complex_float * y)
              int gsl_blas_zgemv (CBLAS_TRANSPOSE_t TransA, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_vector_complex * x, const gsl_complex beta, gsl_vector_complex * y)

   These functions compute the matrix-vector product and sum :math:`y = \alpha op(A) x + \beta y`,
   where :math:`op(A) = A`, :math:`A^T`, :math:`A^H` for :data:`TransA` = :code:`CblasNoTrans`,
   :code:`CblasTrans`, :code:`CblasConjTrans`.

.. index::
   single: TRMV, Level-2 BLAS

.. function:: int gsl_blas_strmv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_float * A, gsl_vector_float * x)
              int gsl_blas_dtrmv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix * A, gsl_vector * x)
              int gsl_blas_ctrmv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_complex_float * A, gsl_vector_complex_float * x)
              int gsl_blas_ztrmv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_complex * A, gsl_vector_complex * x)

   These functions compute the matrix-vector product 
   :math:`x = op(A) x` for the triangular matrix :data:`A`, where
   :math:`op(A) = A`, :math:`A^T`, :math:`A^H` for :data:`TransA` =
   :code:`CblasNoTrans`, :code:`CblasTrans`, :code:`CblasConjTrans`.  When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle of :data:`A` is
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   of :data:`A` is used.  If :data:`Diag` is :code:`CblasNonUnit` then the
   diagonal of the matrix is used, but if :data:`Diag` is :code:`CblasUnit`
   then the diagonal elements of the matrix :data:`A` are taken as unity and
   are not referenced.

.. index::
   single: TRSV, Level-2 BLAS

.. function:: int gsl_blas_strsv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_float * A, gsl_vector_float * x)
              int gsl_blas_dtrsv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix * A, gsl_vector * x)
              int gsl_blas_ctrsv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_complex_float * A, gsl_vector_complex_float * x)
              int gsl_blas_ztrsv (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_matrix_complex * A, gsl_vector_complex * x)

   These functions compute :math:`inv(op(A)) x` for :data:`x`, where
   :math:`op(A) = A`, :math:`A^T`, :math:`A^H` for :data:`TransA` =
   :code:`CblasNoTrans`, :code:`CblasTrans`, :code:`CblasConjTrans`.  When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle of :data:`A` is
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   of :data:`A` is used.  If :data:`Diag` is :code:`CblasNonUnit` then the
   diagonal of the matrix is used, but if :data:`Diag` is :code:`CblasUnit`
   then the diagonal elements of the matrix :data:`A` are taken as unity and
   are not referenced.

.. index::
   single: SYMV, Level-2 BLAS

.. function:: int gsl_blas_ssymv (CBLAS_UPLO_t Uplo, float alpha, const gsl_matrix_float * A, const gsl_vector_float * x, float beta, gsl_vector_float * y)
              int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)

   These functions compute the matrix-vector product and sum :math:`y = \alpha A x + \beta y`
   for the symmetric matrix :data:`A`.  Since the
   matrix :data:`A` is symmetric only its upper half or lower half need to be
   stored.  When :data:`Uplo` is :code:`CblasUpper` then the upper triangle
   and diagonal of :data:`A` are used, and when :data:`Uplo` is
   :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.

.. index::
   single: HEMV, Level-2 BLAS

.. function:: int gsl_blas_chemv (CBLAS_UPLO_t Uplo, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_vector_complex_float * x, const gsl_complex_float beta, gsl_vector_complex_float * y)
              int gsl_blas_zhemv (CBLAS_UPLO_t Uplo, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_vector_complex * x, const gsl_complex beta, gsl_vector_complex * y)

   These functions compute the matrix-vector product and sum :math:`y = \alpha A x + \beta y`
   for the hermitian matrix :data:`A`.  Since the
   matrix :data:`A` is hermitian only its upper half or lower half need to be
   stored.  When :data:`Uplo` is :code:`CblasUpper` then the upper triangle
   and diagonal of :data:`A` are used, and when :data:`Uplo` is
   :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.  The imaginary elements of the diagonal are automatically assumed
   to be zero and are not referenced.

.. index::
   single: GER, Level-2 BLAS
   single: GERU, Level-2 BLAS

.. function:: int gsl_blas_sger (float alpha, const gsl_vector_float * x, const gsl_vector_float * y, gsl_matrix_float * A)
              int gsl_blas_dger (double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A)
              int gsl_blas_cgeru (const gsl_complex_float alpha, const gsl_vector_complex_float * x, const gsl_vector_complex_float * y, gsl_matrix_complex_float * A)
              int gsl_blas_zgeru (const gsl_complex alpha, const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_matrix_complex * A)

   These functions compute the rank-1 update :math:`A = \alpha x y^T + A` of
   the matrix :data:`A`.

.. index::
   single: GERC, Level-2 BLAS

.. function:: int gsl_blas_cgerc (const gsl_complex_float alpha, const gsl_vector_complex_float * x, const gsl_vector_complex_float * y, gsl_matrix_complex_float * A)
              int gsl_blas_zgerc (const gsl_complex alpha, const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_matrix_complex * A)

   These functions compute the conjugate rank-1 update :math:`A = \alpha x y^H + A`
   of the matrix :data:`A`.

.. index::
   single: SYR, Level-2 BLAS

.. function:: int gsl_blas_ssyr (CBLAS_UPLO_t Uplo, float alpha, const gsl_vector_float * x, gsl_matrix_float * A)
              int gsl_blas_dsyr (CBLAS_UPLO_t Uplo, double alpha, const gsl_vector * x, gsl_matrix * A)

   These functions compute the symmetric rank-1 update :math:`A = \alpha x x^T + A`
   of the symmetric matrix :data:`A`.  Since the matrix :data:`A` is
   symmetric only its upper half or lower half need to be stored.  When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle and diagonal of
   :data:`A` are used, and when :data:`Uplo` is :code:`CblasLower` then the
   lower triangle and diagonal of :data:`A` are used.

.. index::
   single: HER, Level-2 BLAS

.. function:: int gsl_blas_cher (CBLAS_UPLO_t Uplo, float alpha, const gsl_vector_complex_float * x, gsl_matrix_complex_float * A)
              int gsl_blas_zher (CBLAS_UPLO_t Uplo, double alpha, const gsl_vector_complex * x, gsl_matrix_complex * A)

   These functions compute the hermitian rank-1 update :math:`A = \alpha x x^H + A`
   of the hermitian matrix :data:`A`.  Since the matrix :data:`A` is
   hermitian only its upper half or lower half need to be stored.  When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle and diagonal of
   :data:`A` are used, and when :data:`Uplo` is :code:`CblasLower` then the
   lower triangle and diagonal of :data:`A` are used.  The imaginary elements
   of the diagonal are automatically set to zero.

.. index::
   single: SYR2, Level-2 BLAS

.. function:: int gsl_blas_ssyr2 (CBLAS_UPLO_t Uplo, float alpha, const gsl_vector_float * x, const gsl_vector_float * y, gsl_matrix_float * A)
              int gsl_blas_dsyr2 (CBLAS_UPLO_t Uplo, double alpha, const gsl_vector * x, const gsl_vector * y, gsl_matrix * A)

   These functions compute the symmetric rank-2 update :math:`A = \alpha x y^T + \alpha y x^T + A`
   of the symmetric matrix :data:`A`.  Since the
   matrix :data:`A` is symmetric only its upper half or lower half need to be
   stored.  When :data:`Uplo` is :code:`CblasUpper` then the upper triangle
   and diagonal of :data:`A` are used, and when :data:`Uplo` is
   :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.

.. index::
   single: HER2, Level-2 BLAS

.. function:: int gsl_blas_cher2 (CBLAS_UPLO_t Uplo, const gsl_complex_float alpha, const gsl_vector_complex_float * x, const gsl_vector_complex_float * y, gsl_matrix_complex_float * A)
              int gsl_blas_zher2 (CBLAS_UPLO_t Uplo, const gsl_complex alpha, const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_matrix_complex * A)

   These functions compute the hermitian rank-2 update :math:`A = \alpha x y^H + \alpha^* y x^H + A`
   of the hermitian matrix :data:`A`.  Since the
   matrix :data:`A` is hermitian only its upper half or lower half need to be
   stored.  When :data:`Uplo` is :code:`CblasUpper` then the upper triangle
   and diagonal of :data:`A` are used, and when :data:`Uplo` is
   :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.  The imaginary elements of the diagonal are automatically set to zero.

Level 3
-------

.. index::
   single: GEMM, Level-3 BLAS

.. function:: int gsl_blas_sgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, float alpha, const gsl_matrix_float * A, const gsl_matrix_float * B, float beta, gsl_matrix_float * C)
              int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
              int gsl_blas_cgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, const gsl_complex_float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, const gsl_complex beta, gsl_matrix_complex * C)

   These functions compute the matrix-matrix product and sum :math:`C = \alpha op(A) op(B) + \beta C`
   where :math:`op(A) = A`, :math:`A^T`,
   :math:`A^H` for :data:`TransA` = :code:`CblasNoTrans`, :code:`CblasTrans`,
   :code:`CblasConjTrans` and similarly for the parameter :data:`TransB`.

.. index::
   single: SYMM, Level-3 BLAS

.. function:: int gsl_blas_ssymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, float alpha, const gsl_matrix_float * A, const gsl_matrix_float * B, float beta, gsl_matrix_float * C)
              int gsl_blas_dsymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
              int gsl_blas_csymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, const gsl_complex_float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zsymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, const gsl_complex beta, gsl_matrix_complex * C)

   These functions compute the matrix-matrix product and sum :math:`C = \alpha A B + \beta C`
   for :data:`Side` is :code:`CblasLeft` and :math:`C = \alpha B A + \beta C`
   for :data:`Side` is :code:`CblasRight`, where the
   matrix :data:`A` is symmetric.  When :data:`Uplo` is :code:`CblasUpper` then
   the upper triangle and diagonal of :data:`A` are used, and when :data:`Uplo`
   is :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.

.. index::
   single: HEMM, Level-3 BLAS

.. function:: int gsl_blas_chemm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, const gsl_complex_float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zhemm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, const gsl_complex beta, gsl_matrix_complex * C)

   These functions compute the matrix-matrix product and sum :math:`C = \alpha A B + \beta C`
   for :data:`Side` is :code:`CblasLeft` and :math:`C = \alpha B A + \beta C`
   for :data:`Side` is :code:`CblasRight`, where the
   matrix :data:`A` is hermitian.  When :data:`Uplo` is :code:`CblasUpper` then
   the upper triangle and diagonal of :data:`A` are used, and when :data:`Uplo`
   is :code:`CblasLower` then the lower triangle and diagonal of :data:`A` are
   used.  The imaginary elements of the diagonal are automatically set to
   zero.

.. index::
   single: TRMM, Level-3 BLAS

.. function:: int gsl_blas_strmm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, float alpha, const gsl_matrix_float * A, gsl_matrix_float * B)
              int gsl_blas_dtrmm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, double alpha, const gsl_matrix * A, gsl_matrix * B)
              int gsl_blas_ctrmm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, gsl_matrix_complex_float * B)
              int gsl_blas_ztrmm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_complex alpha, const gsl_matrix_complex * A, gsl_matrix_complex * B)

   These functions compute the matrix-matrix product :math:`B = \alpha op(A) B`
   for :data:`Side` is :code:`CblasLeft` and :math:`B = \alpha B op(A)` for
   :data:`Side` is :code:`CblasRight`.  The matrix :data:`A` is triangular and
   :math:`op(A) = A`, :math:`A^T`, :math:`A^H` for :data:`TransA` =
   :code:`CblasNoTrans`, :code:`CblasTrans`, :code:`CblasConjTrans`. When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle of :data:`A` is
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   of :data:`A` is used.  If :data:`Diag` is :code:`CblasNonUnit` then the
   diagonal of :data:`A` is used, but if :data:`Diag` is :code:`CblasUnit` then
   the diagonal elements of the matrix :data:`A` are taken as unity and are
   not referenced.

.. index::
   single: TRSM, Level-3 BLAS

.. function:: int gsl_blas_strsm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, float alpha, const gsl_matrix_float * A, gsl_matrix_float * B)
              int gsl_blas_dtrsm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, double alpha, const gsl_matrix * A, gsl_matrix * B)
              int gsl_blas_ctrsm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, gsl_matrix_complex_float * B)
              int gsl_blas_ztrsm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const gsl_complex alpha, const gsl_matrix_complex * A, gsl_matrix_complex * B)

   These functions compute the inverse-matrix matrix product 
   :math:`B = \alpha op(inv(A))B` for :data:`Side` is 
   :code:`CblasLeft` and :math:`B = \alpha B op(inv(A))` for
   :data:`Side` is :code:`CblasRight`.  The matrix :data:`A` is triangular and
   :math:`op(A) = A`, :math:`A^T`, :math:`A^H` for :data:`TransA` =
   :code:`CblasNoTrans`, :code:`CblasTrans`, :code:`CblasConjTrans`. When
   :data:`Uplo` is :code:`CblasUpper` then the upper triangle of :data:`A` is
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   of :data:`A` is used.  If :data:`Diag` is :code:`CblasNonUnit` then the
   diagonal of :data:`A` is used, but if :data:`Diag` is :code:`CblasUnit` then
   the diagonal elements of the matrix :data:`A` are taken as unity and are
   not referenced.

.. index::
   single: SYRK, Level-3 BLAS

.. function:: int gsl_blas_ssyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, float alpha, const gsl_matrix_float * A, float beta, gsl_matrix_float * C)
              int gsl_blas_dsyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix * A, double beta, gsl_matrix * C)
              int gsl_blas_csyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_complex_float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zsyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_complex beta, gsl_matrix_complex * C)

   These functions compute a rank-k update of the symmetric matrix :data:`C`,
   :math:`C = \alpha A A^T + \beta C` when :data:`Trans` is
   :code:`CblasNoTrans` and :math:`C = \alpha A^T A + \beta C` when
   :data:`Trans` is :code:`CblasTrans`.  Since the matrix :data:`C` is symmetric
   only its upper half or lower half need to be stored.  When :data:`Uplo` is
   :code:`CblasUpper` then the upper triangle and diagonal of :data:`C` are
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   and diagonal of :data:`C` are used.

.. index::
   single: HERK, Level-3 BLAS

.. function:: int gsl_blas_cherk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, float alpha, const gsl_matrix_complex_float * A, float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zherk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix_complex * A, double beta, gsl_matrix_complex * C)

   These functions compute a rank-k update of the hermitian matrix :data:`C`,
   :math:`C = \alpha A A^H + \beta C` when :data:`Trans` is
   :code:`CblasNoTrans` and :math:`C = \alpha A^H A + \beta C` when
   :data:`Trans` is :code:`CblasConjTrans`.  Since the matrix :data:`C` is hermitian
   only its upper half or lower half need to be stored.  When :data:`Uplo` is
   :code:`CblasUpper` then the upper triangle and diagonal of :data:`C` are
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   and diagonal of :data:`C` are used.  The imaginary elements of the
   diagonal are automatically set to zero.

.. index::
   single: SYR2K, Level-3 BLAS

.. function:: int gsl_blas_ssyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, float alpha, const gsl_matrix_float * A, const gsl_matrix_float * B, float beta, gsl_matrix_float * C)
              int gsl_blas_dsyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
              int gsl_blas_csyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, const gsl_complex_float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zsyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, const gsl_complex beta, gsl_matrix_complex * C)

   These functions compute a rank-2k update of the symmetric matrix :data:`C`,
   :math:`C = \alpha A B^T + \alpha B A^T + \beta C` when :data:`Trans` is
   :code:`CblasNoTrans` and :math:`C = \alpha A^T B + \alpha B^T A + \beta C` when
   :data:`Trans` is :code:`CblasTrans`.  Since the matrix :data:`C` is symmetric
   only its upper half or lower half need to be stored.  When :data:`Uplo` is
   :code:`CblasUpper` then the upper triangle and diagonal of :data:`C` are
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   and diagonal of :data:`C` are used.

.. index::
   single: HER2K, Level-3 BLAS

.. function:: int gsl_blas_cher2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, float beta, gsl_matrix_complex_float * C)
              int gsl_blas_zher2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, double beta, gsl_matrix_complex * C)

   These functions compute a rank-2k update of the hermitian matrix :data:`C`,
   :math:`C = \alpha A B^H + \alpha^* B A^H + \beta C` when :data:`Trans` is
   :code:`CblasNoTrans` and :math:`C = \alpha A^H B + \alpha^* B^H A + \beta C` when
   :data:`Trans` is :code:`CblasConjTrans`.  Since the matrix :data:`C` is hermitian
   only its upper half or lower half need to be stored.  When :data:`Uplo` is
   :code:`CblasUpper` then the upper triangle and diagonal of :data:`C` are
   used, and when :data:`Uplo` is :code:`CblasLower` then the lower triangle
   and diagonal of :data:`C` are used.  The imaginary elements of the
   diagonal are automatically set to zero.

Examples
========

The following program computes the product of two matrices using the
Level-3 BLAS function DGEMM,

.. only:: not texinfo

   .. math::

      \left(
      \begin{matrix}
        0.11&0.12&0.13 \\
        0.21&0.22&0.23
      \end{matrix}
      \right)
      \left(
      \begin{matrix}
        1011&1012 \\
        1021&1022 \\
        1031&1031
      \end{matrix}
      \right)
      =
      \left(
      \begin{matrix}
        367.76&368.12 \\
        674.06&674.72
      \end{matrix}
      \right)

.. only:: texinfo

   ::

      [ 0.11 0.12 0.13 ]  [ 1011 1012 ]     [ 367.76 368.12 ]
      [ 0.21 0.22 0.23 ]  [ 1021 1022 ]  =  [ 674.06 674.72 ]
                          [ 1031 1032 ]

The matrices are stored in row major order, according to the C convention 
for arrays.

.. include:: examples/blas.c
   :code:

Here is the output from the program,

.. include:: examples/blas.txt
   :code:

.. _sec_blas-references:

References and Further Reading
==============================

Information on the BLAS standards, including both the legacy and
updated interface standards, is available online from the BLAS
Homepage and BLAS Technical Forum web-site.

* BLAS Homepage, http://www.netlib.org/blas/

* BLAS Technical Forum, http://www.netlib.org/blas/blast-forum/

The following papers contain the specifications for Level 1, Level 2 and
Level 3 BLAS.

* C. Lawson, R. Hanson, D. Kincaid, F. Krogh, "Basic Linear Algebra
  Subprograms for Fortran Usage", ACM Transactions on Mathematical
  Software, Vol.: 5 (1979), Pages 308--325.

* J.J. Dongarra, J. DuCroz, S. Hammarling, R. Hanson, "An Extended Set of
  Fortran Basic Linear Algebra Subprograms", ACM Transactions on
  Mathematical Software, Vol.: 14, No.: 1 (1988), Pages 1--32.

* J.J. Dongarra, I. Duff, J. DuCroz, S. Hammarling, "A Set of
  Level 3 Basic Linear Algebra Subprograms", ACM Transactions on
  Mathematical Software, Vol.: 16 (1990), Pages 1--28.

Postscript versions of the latter two papers are available from
http://www.netlib.org/blas/. A CBLAS wrapper for Fortran BLAS
libraries is available from the same location.

.. rubric:: Footnotes

.. [#f1] In the low-level CBLAS interface, a negative stride accesses the vector elements
         in reverse order, i.e. the :math:`i`-th element is given by
         :math:`(N-i)*|incx|` for :math:`incx < 0`.
