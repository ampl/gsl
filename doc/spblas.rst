.. index::
   single: sparse BLAS
   single: BLAS, sparse

*******************
Sparse BLAS Support
*******************

.. include:: include.rst

The Sparse Basic Linear Algebra Subprograms (|blas|) define a set of
fundamental operations on vectors and sparse matrices which can be used
to create optimized higher-level linear algebra functionality.
GSL supports a limited number of BLAS operations for sparse matrices.

The header file :file:`gsl_spblas.h` contains the prototypes for the
sparse BLAS functions and related declarations.

.. index::
   single: sparse matrices, BLAS operations

Sparse BLAS operations
======================

.. function:: int gsl_spblas_dgemv (const CBLAS_TRANSPOSE_t TransA, const double alpha, const gsl_spmatrix * A, const gsl_vector * x, const double beta, gsl_vector * y)

   This function computes the matrix-vector product and sum
   :math:`y \leftarrow \alpha op(A) x + \beta y`, where
   :math:`op(A) = A, A^T` for :data:`TransA` = :code:`CblasNoTrans`,
   :code:`CblasTrans`. In-place computations are not supported, so
   :data:`x` and :data:`y` must be distinct vectors.
   The matrix :data:`A` may be in triplet or compressed format.

.. function:: int gsl_spblas_dgemm (const double alpha, const gsl_spmatrix * A, const gsl_spmatrix * B, gsl_spmatrix * C)

   This function computes the sparse matrix-matrix product
   :math:`C = \alpha A B`. The matrices must be in compressed format.

.. index::
   single: sparse BLAS, references

References and Further Reading
==============================

The algorithms used by these functions are described in the
following sources:

* Davis, T. A., Direct Methods for Sparse Linear Systems, SIAM, 2006.

* CSparse software library, https://www.cise.ufl.edu/research/sparse/CSparse
