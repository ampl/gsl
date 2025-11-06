/* linalg/rqrc.c
 * 
 * Copyright (C) 2020 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/*
 * this module contains routines for the QR factorization of a matrix
 * using the recursive Level 3 BLAS algorithm of Elmroth and Gustavson with
 * additional modifications courtesy of Julien Langou.
 */

static int unpack_Q1(gsl_matrix_complex * Q);
static int unpack_Q2(const gsl_matrix_complex * QR, const gsl_matrix_complex * T, gsl_matrix_complex * Q);
static int aux_ULH(const gsl_matrix_complex * L, gsl_matrix_complex * U);
static int aux_mLU(gsl_matrix_complex * A);
static int aux_ApUBH(const gsl_matrix_complex * U, const gsl_matrix_complex * B, gsl_matrix_complex * A);

/*
gsl_linalg_complex_QR_decomp_r()
  QR decomposition using Level 3 BLAS recursive algorithm of:
  
Elmroth, E. and Gustavson, F.G., 2000. Applying recursion to serial and parallel
  QR factorization leads to better performance. IBM Journal of Research and Development,
  44(4), pp.605-624.

Inputs: A - matrix to be factored, M-by-N with M >= N
        T - N-by-N upper triangular factor of block reflector

Return: success/error

Notes:
1) on output, upper triangle of A contains R; elements below the diagonal
are columns of V. The matrix Q is

Q = I - V T V^H

where T is upper triangular. Note that diag(T) = tau
*/

int
gsl_linalg_complex_QR_decomp_r (gsl_matrix_complex * A, gsl_matrix_complex * T)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (T->size1 != T->size2)
    {
      GSL_ERROR ("T matrix must be square", GSL_ENOTSQR);
    }
  else if (T->size1 != N)
    {
      GSL_ERROR ("T matrix does not match dimensions of A", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* base case, compute householder transform for single column matrix */

      gsl_complex * T00 = gsl_matrix_complex_ptr(T, 0, 0);
      gsl_vector_complex_view v = gsl_matrix_complex_column(A, 0);

      *T00 = gsl_linalg_complex_householder_transform(&v.vector);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * partition matrices:
       *
       *       N1  N2              N1  N2
       * N1 [ A11 A12 ] and  N1 [ T11 T12 ]
       * M2 [ A21 A22 ]      N2 [  0  T22 ]
       */
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;
      const size_t M2 = M - N1;

      gsl_matrix_complex_view A11 = gsl_matrix_complex_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_complex_view A12 = gsl_matrix_complex_submatrix(A, 0, N1, N1, N2);
      gsl_matrix_complex_view A21 = gsl_matrix_complex_submatrix(A, N1, 0, M2, N1);
      gsl_matrix_complex_view A22 = gsl_matrix_complex_submatrix(A, N1, N1, M2, N2);

      gsl_matrix_complex_view T11 = gsl_matrix_complex_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_complex_view T12 = gsl_matrix_complex_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_complex_view T22 = gsl_matrix_complex_submatrix(T, N1, N1, N2, N2);

      gsl_matrix_complex_view m;

      /*
       * Eq. 2: recursively factor
       *
       * [ A11 ] = Q1 [ R11 ]
       * [ A21 ]      [  0  ]
       *                                          N1
       * Note: Q1 = I - V1 T11 V1^H, where V1 = [ V11 ] N1
       *                                        [ V21 ] M2
       */
      m = gsl_matrix_complex_submatrix(A, 0, 0, M, N1);
      status = gsl_linalg_complex_QR_decomp_r(&m.matrix, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       * [ R12 ] := Q1^H [ A12 ] = [ A12 ] - [ V11 W ]
       * [ A22 ]         [ A22 ]   [ A22 ]   [ V21 W ]
       *
       * where W = T11^H (V11^H A12 + V21^H A22), and using T12 as temporary storage
       */
      gsl_matrix_complex_memcpy(&T12.matrix, &A12.matrix);                                            /* W := A12 */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasConjTrans, CblasUnit, GSL_COMPLEX_ONE,
                     &A11.matrix, &T12.matrix);                                                       /* W := V11^H * A12 */
      gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, &A21.matrix, &A22.matrix,
                     GSL_COMPLEX_ONE, &T12.matrix);                                                   /* W := W + V21^H * A22 */
      gsl_blas_ztrmm(CblasLeft, CblasUpper, CblasConjTrans, CblasNonUnit, GSL_COMPLEX_ONE,
                     &T11.matrix, &T12.matrix);                                                       /* W := T11^H * W */
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_NEGONE, &A21.matrix, &T12.matrix,
                     GSL_COMPLEX_ONE, &A22.matrix);                                                   /* A22 = A22 - V21 * W */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_ONE, &A11.matrix,
                     &T12.matrix);                                                                    /* tmp = V11 * W */
      gsl_matrix_complex_sub(&A12.matrix, &T12.matrix);                                               /* R12 := A12 - V11 * W */

      /*
       * Eq. 4: recursively factor
       *
       * A22 = Q2~ R22
       *
       *              N1 M2
       * Note: Q2 = [ I   0  ] N1
       *            [ 0  Q2~ ] M2
       */
      status = gsl_linalg_complex_QR_decomp_r(&A22.matrix, &T22.matrix);
      if (status)
        return status;

      /*
       * Eq. 13: update T12 := -T11 * V1^H * V2 * T22
       *
       * where:
       *
       *        N1                N2
       * V1 = [ V11 ] N1   V2 = [  0  ] N1
       *      [ V21 ] N2        [ V22 ] N2
       *      [ V31 ] M-N       [ V32 ] M-N
       *
       * Note: V1^H V2 = V21^H V22 + V31^H V32
       * Also, V11, V22 are unit lower triangular
       */

      m = gsl_matrix_complex_submatrix(&A21.matrix, 0, 0, N2, N1);                                 /* V21 */
      gsl_matrix_complex_conjtrans_memcpy(&T12.matrix, &m.matrix);                                 /* T12 := V21^H */

      m = gsl_matrix_complex_submatrix(A, N1, N1, N2, N2);                                         /* V22 */
      gsl_blas_ztrmm(CblasRight, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_ONE,
                     &m.matrix, &T12.matrix);                                                      /* T12 := V21^H * V22 */

      if (M > N)
        {
          gsl_matrix_complex_view V31 = gsl_matrix_complex_submatrix(A, N, 0, M - N, N1);
          gsl_matrix_complex_view V32 = gsl_matrix_complex_submatrix(A, N, N1, M - N, N2);

          gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, &V31.matrix, &V32.matrix,
                         GSL_COMPLEX_ONE, &T12.matrix);                                            /* T12 := T12 + V31^T * V32 */
        }

      gsl_blas_ztrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_NEGONE,
                     &T11.matrix, &T12.matrix);                                                    /* T12 := -T11 * T12 */
      gsl_blas_ztrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_ONE,
                     &T22.matrix, &T12.matrix);                                                    /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
}

/* Solves the square system A x = b for x using the QR factorisation,
 *
 *  R x = Q^H b
 *
 * where Q = I - V T V^H
 */

int
gsl_linalg_complex_QR_solve_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                               const gsl_vector_complex * b, gsl_vector_complex * x)
{
  const size_t N = QR->size2;

  if (QR->size1 != N)
    {
      GSL_ERROR ("QR matrix must be square", GSL_ENOTSQR);
    }
  else if (T->size1 != QR->size1 || T->size2 != QR->size2)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (N != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (N != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute Q^H b = [I - V T^H V^H] b */

      /* x := V^H b */
      gsl_vector_complex_memcpy(x, b);
      gsl_blas_ztrmv(CblasLower, CblasConjTrans, CblasUnit, QR, x);

      /* x = T^H * x */
      gsl_blas_ztrmv(CblasUpper, CblasConjTrans, CblasNonUnit, T, x);

      /* x = V * x */
      gsl_blas_ztrmv(CblasLower, CblasNoTrans, CblasUnit, QR, x);

      /* x = b - V * x */
      for (i = 0; i < N; ++i)
        {
          gsl_complex * xi = gsl_vector_complex_ptr(x, i);
          gsl_complex bi = gsl_vector_complex_get(b, i);
          GSL_REAL(*xi) = GSL_REAL(bi) - GSL_REAL(*xi);
          GSL_IMAG(*xi) = GSL_IMAG(bi) - GSL_IMAG(*xi);
        }

      /* Solve R x = Q^H b, storing x in-place */
      gsl_blas_ztrsv (CblasUpper, CblasNoTrans, CblasNonUnit, QR, x);

      return GSL_SUCCESS;
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   A x = b 
 *  
 * for M >= N using the QR factorization A = Q R. 
 *
 * Inputs: QR   - [R; V] matrix, M-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         b    - right hand side, size M
 *         x    - (output) solution, size M
 *                x(1:N) = least squares solution vector
 *                x(N+1:M) = vector whose norm equals ||b - Ax||
 *         work - workspace, size N
 */

int
gsl_linalg_complex_QR_lssolve_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                                 const gsl_vector_complex * b, gsl_vector_complex * x, gsl_vector_complex * work)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M < N)
    {
      GSL_ERROR ("QR matrix must have M >= N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (M != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (M != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (N != work->size)
    {
      GSL_ERROR ("matrix size must match work size", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view R = gsl_matrix_complex_const_submatrix (QR, 0, 0, N, N);
      gsl_vector_complex_view x1 = gsl_vector_complex_subvector(x, 0, N);

      /* compute x = Q^H b */
      gsl_vector_complex_memcpy(x, b);
      gsl_linalg_complex_QR_QHvec_r (QR, T, x, work);

      /* Solve R x = Q^H b */
      gsl_blas_ztrsv (CblasUpper, CblasNoTrans, CblasNonUnit, &R.matrix, &x1.vector);

      return GSL_SUCCESS;
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   A x = B
 *  
 * for M >= N using the QR factorization A = Q R. 
 *
 * Inputs: QR   - [R; V] matrix, M-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         B    - right hand sides, size M-by-nrhs
 *         X    - (output) solution, size M-by-nrhs
 *                x(1:N,k)   = least squares solution vector for B(:,k)
 *                x(N+1:M,k) = vector whose norm equals || B(:,k) - A X(1:N,k) ||
 *         work - workspace, size N-by-nrhs
 */

int
gsl_linalg_complex_QR_lssolvem_r (const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                                  const gsl_matrix_complex * B, gsl_matrix_complex * X, gsl_matrix_complex * work)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;
  const size_t nrhs = B->size2;

  if (M < N)
    {
      GSL_ERROR ("QR matrix must have M >= N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (B->size1 != M)
    {
      GSL_ERROR ("matrix size must match B size", GSL_EBADLEN);
    }
  else if (X->size1 != M || X->size2 != nrhs)
    {
      GSL_ERROR ("solution matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (work->size1 != N || work->size2 != nrhs)
    {
      GSL_ERROR ("work matrix has wrong dimensions", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view R = gsl_matrix_complex_const_submatrix (QR, 0, 0, N, N);
      gsl_matrix_complex_view X1 = gsl_matrix_complex_submatrix(X, 0, 0, N, nrhs);

      /* compute X = Q^H B */
      gsl_matrix_complex_memcpy(X, B);
      gsl_linalg_complex_QR_QHmat_r (QR, T, X, work);

      /* Solve R X = Q^H B */
      gsl_blas_ztrsm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
                      GSL_COMPLEX_ONE, &R.matrix, &X1.matrix);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_complex_QR_unpack_r()
  Unpack matrices Q and R

Inputs: QR - packed QR format, M-by-N
        T  - block reflector matrix, N-by-N
        Q  - (output) Q matrix, M-by-M
        R  - (output) R matrix, N-by-N

Return: success/error

Notes:
1) Implementation provided by Julien Langou
*/

int
gsl_linalg_complex_QR_unpack_r(const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                               gsl_matrix_complex * Q, gsl_matrix_complex * R)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (Q->size1 != M || Q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M-by-M", GSL_EBADLEN);
    }
  else if (R->size1 != N || R->size2 != N)
    {
      GSL_ERROR ("R matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view RV = gsl_matrix_complex_const_submatrix(QR, 0, 0, N, N);
      gsl_matrix_complex_view Q1 = gsl_matrix_complex_submatrix(Q, 0, 0, M, N);
      gsl_matrix_complex_view m;

      /*
       * set Q1 = [ T ]
       *          [ V ]
       */
      m = gsl_matrix_complex_submatrix(Q, 0, 0, N, N);
      gsl_matrix_complex_tricpy(CblasUpper, CblasNonUnit, &m.matrix, T);
      gsl_matrix_complex_tricpy(CblasLower, CblasUnit, &m.matrix, &RV.matrix);

      if (M > N)
        {
          gsl_matrix_complex_const_view tmp = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);
          m = gsl_matrix_complex_submatrix(Q, N, 0, M - N, N);
          gsl_matrix_complex_memcpy(&m.matrix, &tmp.matrix);
        }

      unpack_Q1(&Q1.matrix);

      if (M > N)
        {
          gsl_matrix_complex_view Q2 = gsl_matrix_complex_submatrix(Q, 0, N, M, M - N);
          unpack_Q2(QR, T, &Q2.matrix);
        }

      /* copy R */
      gsl_matrix_complex_tricpy(CblasUpper, CblasNonUnit, R, &RV.matrix);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_QR_QHvec_r()
  Apply M-by-M Q^H to the M-by-1 vector b

Inputs: QR   - [R; V] matrix encoded by gsl_linalg_complex_QR_decomp_r
        T    - block reflector matrix
        b    - M-by-1 vector replaced by Q^H b on output
        work - workspace, length N

Notes:
1) Q^H b = (I - V T^H V^H) b
         = b - V T^H [ V1^H V2^H ] [ b1 ]
                                   [ b2 ]
         = b - V T^H [ V1^H b1 + V2^H b2 ]
         = [ b1 ] - [ V1 w ]
           [ b2 ]   [ V2 w ]

where w = T^H ( V1^H b1 + V2^H b2 )
*/

int
gsl_linalg_complex_QR_QHvec_r(const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                              gsl_vector_complex * b, gsl_vector_complex * work)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (b->size != M)
    {
      GSL_ERROR ("b vector must have length M", GSL_EBADLEN);
    }
  else if (work->size != N)
    {
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view V1 = gsl_matrix_complex_const_submatrix(QR, 0, 0, N, N);
      gsl_vector_complex_view b1 = gsl_vector_complex_subvector(b, 0, N);
      gsl_vector_complex_view b2;

      /* work := V1^H b1 */
      gsl_vector_complex_memcpy(work, &b1.vector);
      gsl_blas_ztrmv(CblasLower, CblasConjTrans, CblasUnit, &V1.matrix, work);

      if (M > N)
        {
          gsl_matrix_complex_const_view V2 = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);

          /* work = work + V2^H b2 */
          b2 = gsl_vector_complex_subvector(b, N, M - N);
          gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, &V2.matrix, &b2.vector, GSL_COMPLEX_ONE, work);
        }

      /* work = T^H * work */
      gsl_blas_ztrmv(CblasUpper, CblasConjTrans, CblasNonUnit, T, work);

      if (M > N)
        {
          /* b2 = b2 - V2 * work */
          gsl_matrix_complex_const_view V2 = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);
          gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_NEGONE, &V2.matrix, work, GSL_COMPLEX_ONE, &b2.vector);
        }

      /* b1 = b1 - V1 * work */
      gsl_blas_ztrmv(CblasLower, CblasNoTrans, CblasUnit, &V1.matrix, work);
      gsl_vector_complex_sub(&b1.vector, work);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_QR_QHmat_r()
  Apply M-by-M Q^H to the M-by-K matrix B

Inputs: QR   - [R; V] matrix encoded by gsl_linalg_complex_QR_decomp_r
        T    - block reflector matrix
        B    - M-by-K matrix replaced by Q^H B on output
        work - N-by-K workspace

Notes:
1) Provided by Julien Langou
*/

int
gsl_linalg_complex_QR_QHmat_r(const gsl_matrix_complex * QR, const gsl_matrix_complex * T,
                              gsl_matrix_complex * B, gsl_matrix_complex * work)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;
  const size_t K = B->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (B->size1 != M)
    {
      GSL_ERROR ("B matrix must have M rows", GSL_EBADLEN);
    }
  else if (work->size1 != N || work->size2 != K)
    {
      GSL_ERROR ("workspace must be N-by-K", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view V1 = gsl_matrix_complex_const_submatrix(QR, 0, 0, N, N);
      gsl_matrix_complex_view B1 = gsl_matrix_complex_submatrix(B, 0, 0, N, K);
      gsl_matrix_complex_view B2;

      /* work := V1^H B1 */
      gsl_matrix_complex_memcpy(work, &B1.matrix);
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasConjTrans, CblasUnit,
                     GSL_COMPLEX_ONE, &V1.matrix, work);

      if (M > N)
        {
          gsl_matrix_complex_const_view V2 = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);

          /* work = work + V2^H B2 */
          B2 = gsl_matrix_complex_submatrix(B, N, 0, M - N, K);
          gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE,
                         &V2.matrix, &B2.matrix, GSL_COMPLEX_ONE, work);
        }

      /* work = T^H * work */
      gsl_blas_ztrmm(CblasLeft, CblasUpper, CblasConjTrans, CblasNonUnit, GSL_COMPLEX_ONE, T, work);

      if (M > N)
        {
          /* B2 = B2 - V2 * work */
          gsl_matrix_complex_const_view V2 = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);
          gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_NEGONE,
                         &V2.matrix, work, GSL_COMPLEX_ONE, &B2.matrix);
        }

      /* B1 = B1 - V1 * work */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_ONE, &V1.matrix, work);
      gsl_matrix_complex_sub(&B1.matrix, work);

      return GSL_SUCCESS;
    }
}

/*
unpack_Q1()
  Compute Q_1

Inputs: Q  - on input, contains T in upper triangle and V in lower trapezoid
             on output, contains Q_1
             M-by-N

Return: success/error

Notes:
1)                                                         N
Q1 = [ Q1 Q2 ] [ I_n ] = (I - V T V^H) [ I; 0 ] = [ I - V1 T V1^H ] N
               [  0  ]                            [   - V2 T V1^H ] M - N
*/

static int
unpack_Q1(gsl_matrix_complex * Q)
{
  int status;
  const size_t M = Q->size1;
  const size_t N = Q->size2;
  gsl_matrix_complex_view Q1 = gsl_matrix_complex_submatrix(Q, 0, 0, N, N);
  gsl_vector_complex_view diag = gsl_matrix_complex_diagonal(&Q1.matrix);

  /* Q1 := T V1^H */
  status = aux_ULH(&Q1.matrix, &Q1.matrix);
  if (status)
    return status;

  if (M > N)
    {
      /* compute Q2 := - V2 T V1^H */
      gsl_matrix_complex_view V2 = gsl_matrix_complex_submatrix(Q, N, 0, M - N, N);
      gsl_blas_ztrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_NEGONE, &Q1.matrix, &V2.matrix);
    }

  /* Q1 := - V1 T V1^H */
  status = aux_mLU(&Q1.matrix);
  if (status)
    return status;

  /* Q1 := I - V1 T V1^H */
  gsl_vector_complex_add_constant(&diag.vector, GSL_COMPLEX_ONE);

  return GSL_SUCCESS;
}

/*
unpack_Q2()
  Compute Q_2

Inputs: QR - [R; V] from QR_decomp_r, M-by-N
        T  - upper triangular T factor, N-by-N
        Q  - (output) Q_2 factor, M-by-(M-N)

Return: success/error

Notes:                         N   M-N
1) Since Q = I - V T V^H = M [ Q1  Q2  ], we have

           M-N
Q2 = Q [    0    ] N
       [ I_{M-N} ] M-N

So, Q2 = Q [ 0; I ] = (I - V T V^H) [ 0; I ] = [   - V1 T V2^H ]
                                               [ I - V2 T V2^H ]
*/

static int
unpack_Q2(const gsl_matrix_complex * QR, const gsl_matrix_complex * T, gsl_matrix_complex * Q)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M <= N)
    {
      GSL_ERROR ("M must be > N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (Q->size1 != M || Q->size2 != (M - N))
    {
      GSL_ERROR ("Q matrix must be M-by-(M-N)", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_complex_const_view V1 = gsl_matrix_complex_const_submatrix(QR, 0, 0, N, N);
      gsl_matrix_complex_const_view V2 = gsl_matrix_complex_const_submatrix(QR, N, 0, M - N, N);
      gsl_matrix_complex_view Q1 = gsl_matrix_complex_submatrix(Q, 0, 0, N, M - N);
      gsl_matrix_complex_view Q2 = gsl_matrix_complex_submatrix(Q, N, 0, M - N, M - N);
      gsl_vector_complex_view diag = gsl_matrix_complex_diagonal(&Q2.matrix);

      /* Q1 := V2^H */
      gsl_matrix_complex_conjtrans_memcpy(&Q1.matrix, &V2.matrix);

      /* Q1 := - T V2^H */
      gsl_blas_ztrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_NEGONE, T, &Q1.matrix);

      /* Q2 := - V2 T V2^H */
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, &V2.matrix, &Q1.matrix, GSL_COMPLEX_ZERO, &Q2.matrix);

      /* Q2 := I - V2 T V2^H */
      gsl_vector_complex_add_constant(&diag.vector, GSL_COMPLEX_ONE);

      /* Q1 := - V1 T V2^H */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_ONE, &V1.matrix, &Q1.matrix);

      return GSL_SUCCESS;
    }
}

/* U := U L^H for triangular matrices L and U; L is unit lower triangular */
static int
aux_ULH(const gsl_matrix_complex * L, gsl_matrix_complex * U)
{
  const size_t N = L->size1;

  if (N != L->size2)
    {
      GSL_ERROR ("L matrix must be square", GSL_ENOTSQR);
    }
  else if (U->size1 != N || U->size2 != N)
    {
      GSL_ERROR ("U matrix must be same size as L", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* nothing to do */
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_complex_const_view L11 = gsl_matrix_complex_const_submatrix(L, 0, 0, N1, N1);
      gsl_matrix_complex_const_view L21 = gsl_matrix_complex_const_submatrix(L, N1, 0, N2, N1);
      gsl_matrix_complex_const_view L22 = gsl_matrix_complex_const_submatrix(L, N1, N1, N2, N2);

      gsl_matrix_complex_view U11 = gsl_matrix_complex_submatrix(U, 0, 0, N1, N1);
      gsl_matrix_complex_view U12 = gsl_matrix_complex_submatrix(U, 0, N1, N1, N2);
      gsl_matrix_complex_view U22 = gsl_matrix_complex_submatrix(U, N1, N1, N2, N2);

      /* U12 = U12 * L22^H */
      gsl_blas_ztrmm(CblasRight, CblasLower, CblasConjTrans, CblasUnit, GSL_COMPLEX_ONE, &L22.matrix, &U12.matrix);

      /* U12 = U12 + U11 * L21^H */
      status = aux_ApUBH(&U11.matrix, &L21.matrix, &U12.matrix);
      if (status)
        return status;

      status = aux_ULH(&L11.matrix, &U11.matrix);
      if (status)
        return status;

      status = aux_ULH(&L22.matrix, &U22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/* store -L*U in A */
static int
aux_mLU(gsl_matrix_complex * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N == 1)
    {
      gsl_complex *A00 = gsl_matrix_complex_ptr(A, 0, 0);
      GSL_REAL(*A00) = -GSL_REAL(*A00);
      GSL_IMAG(*A00) = -GSL_IMAG(*A00);
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_complex_view A11 = gsl_matrix_complex_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_complex_view A12 = gsl_matrix_complex_submatrix(A, 0, N1, N1, N2);
      gsl_matrix_complex_view A21 = gsl_matrix_complex_submatrix(A, N1, 0, N2, N1);
      gsl_matrix_complex_view A22 = gsl_matrix_complex_submatrix(A, N1, N1, N2, N2);

      /* A22 = - L22 U22 */
      status = aux_mLU(&A22.matrix);
      if (status)
        return status;

      /* A22 = A22 - L21 U12 */
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_NEGONE, &A21.matrix, &A12.matrix, GSL_COMPLEX_ONE, &A22.matrix);

      /* A12 - -L11 U12 */
      gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, GSL_COMPLEX_NEGONE, &A11.matrix, &A12.matrix);

      /* A21 = -L21 U11 */
      gsl_blas_ztrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, GSL_COMPLEX_NEGONE, &A11.matrix, &A21.matrix);

      /* A11 = - L11 U11 */
      status = aux_mLU(&A11.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/* A := A + U B^H where U is upper triangular */
static int
aux_ApUBH(const gsl_matrix_complex * U, const gsl_matrix_complex * B, gsl_matrix_complex * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (U->size1 != M || U->size2 != M)
    {
      GSL_ERROR ("U matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (B->size1 != N || B->size2 != M)
    {
      GSL_ERROR ("B matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (M == 1 && N == 1)
    {
      gsl_complex *aptr = gsl_matrix_complex_ptr(A, 0, 0);
      const gsl_complex U00 = gsl_matrix_complex_get(U, 0, 0);
      const gsl_complex B00_conj = gsl_complex_conjugate(gsl_matrix_complex_get(B, 0, 0));
      const gsl_complex prod = gsl_complex_mul(U00, B00_conj);
      GSL_REAL(*aptr) += GSL_REAL(prod);
      GSL_IMAG(*aptr) += GSL_IMAG(prod);
      return GSL_SUCCESS;
    }
  else if (M == 1)
    {
      gsl_complex U00 = gsl_matrix_complex_get(U, 0, 0);
      size_t i;

      for (i = 0; i < N; ++i)
        {
          gsl_complex * ai = gsl_matrix_complex_ptr(A, 0, i);
          gsl_complex bi = gsl_matrix_complex_get(B, i, 0);
          gsl_complex prod = gsl_complex_mul(U00, gsl_complex_conjugate(bi));
          GSL_REAL(*ai) += GSL_REAL(prod);
          GSL_IMAG(*ai) += GSL_IMAG(prod);
        }

      return GSL_SUCCESS;
    }
  else if (N == 1)
    {
      /*
       * partition:
       *
       *         M1  M2
       * B = 1 [ B11 B12 ]
       *
       *       1        M1  M2       1
       * M1 [ A11 ] + [ U11 U12 ] [ B11^H ] M1
       * M2 [ A21 ]   [  0  U22 ] [ B12^H ] M2
       */
      int status;
      const size_t M1 = M / 2;
      const size_t M2 = M - M1;

      gsl_matrix_complex_view A11 = gsl_matrix_complex_submatrix(A, 0, 0, M1, 1);
      gsl_matrix_complex_view A21 = gsl_matrix_complex_submatrix(A, M1, 0, M2, 1);

      gsl_matrix_complex_const_view U11 = gsl_matrix_complex_const_submatrix(U, 0, 0, M1, M1);
      gsl_matrix_complex_const_view U12 = gsl_matrix_complex_const_submatrix(U, 0, M1, M1, M2);
      gsl_matrix_complex_const_view U22 = gsl_matrix_complex_const_submatrix(U, M1, M1, M2, M2);

      gsl_matrix_complex_const_view B11 = gsl_matrix_complex_const_submatrix(B, 0, 0, 1, M1);
      gsl_matrix_complex_const_view B12 = gsl_matrix_complex_const_submatrix(B, 0, M1, 1, M2);

      /* A11 += U12 * B12^H */
      gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, &U12.matrix, &B12.matrix, GSL_COMPLEX_ONE, &A11.matrix);

      /* A11 := A11 + U11 B11^H */
      status = aux_ApUBH(&U11.matrix, &B11.matrix, &A11.matrix);
      if (status)
        return status;

      /* A21 := A21 + U22 B12^H */
      status = aux_ApUBH(&U22.matrix, &B12.matrix, &A21.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t M1 = M / 2;
      const size_t M2 = M - M1;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_complex_view A11 = gsl_matrix_complex_submatrix(A, 0, 0, M1, N1);
      gsl_matrix_complex_view A12 = gsl_matrix_complex_submatrix(A, 0, N1, M1, N2);
      gsl_matrix_complex_view A21 = gsl_matrix_complex_submatrix(A, M1, 0, M2, N1);
      gsl_matrix_complex_view A22 = gsl_matrix_complex_submatrix(A, M1, N1, M2, N2);

      gsl_matrix_complex_const_view U11 = gsl_matrix_complex_const_submatrix(U, 0, 0, M1, M1);
      gsl_matrix_complex_const_view U12 = gsl_matrix_complex_const_submatrix(U, 0, M1, M1, M2);
      gsl_matrix_complex_const_view U22 = gsl_matrix_complex_const_submatrix(U, M1, M1, M2, M2);

      gsl_matrix_complex_const_view B11 = gsl_matrix_complex_const_submatrix(B, 0, 0, N1, M1);
      gsl_matrix_complex_const_view B12 = gsl_matrix_complex_const_submatrix(B, 0, M1, N1, M2);
      gsl_matrix_complex_const_view B21 = gsl_matrix_complex_const_submatrix(B, N1, 0, N2, M1);
      gsl_matrix_complex_const_view B22 = gsl_matrix_complex_const_submatrix(B, N1, M1, N2, M2);

      /* A11 := A11 + U11 B11^H */
      status = aux_ApUBH(&U11.matrix, &B11.matrix, &A11.matrix);
      if (status)
        return status;

      /* A11 := A11 + U12 B12^H */
      gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, &U12.matrix, &B12.matrix, GSL_COMPLEX_ONE, &A11.matrix);

      /* A12 := A12 + U11 B21^H */
      status = aux_ApUBH(&U11.matrix, &B21.matrix, &A12.matrix);
      if (status)
        return status;

      /* A12 := A12 + U12 B22^H */
      gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, &U12.matrix, &B22.matrix, GSL_COMPLEX_ONE, &A12.matrix);

      /* A21 := A21 + U22 B12^H */
      status = aux_ApUBH(&U22.matrix, &B12.matrix, &A21.matrix);
      if (status)
        return status;

      /* A22 := A22 + U22 B22^H */
      status = aux_ApUBH(&U22.matrix, &B22.matrix, &A22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}
