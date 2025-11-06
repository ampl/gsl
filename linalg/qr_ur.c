/* linalg/qr_ur.c
 * 
 * Copyright (C) 2019, 2020 Patrick Alken
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

/*
 * this module contains routines for the QR factorization of a matrix
 * using the recursive Level 3 BLAS algorithm of Elmroth and Gustavson with
 * additional modifications courtesy of Julien Langou.
 */

static double qrtr_householder_transform (double *v0, gsl_vector * v);

/*
gsl_linalg_QR_UR_decomp()
  Compute the QR decomposition of the "triangle on top of rectangle" matrix

  [ S ] = Q [ R ]
  [ A ]     [ 0 ]

where S is N-by-N upper triangular and A is M-by-N dense.

Inputs: S   - on input, upper triangular N-by-N matrix
              on output, R factor in upper triangle
        A   - on input, dense M-by-N matrix
              on output, Householder matrix V
        T   - (output) block reflector matrix, N-by-N

Notes:
1) Based on the Elmroth/Gustavson algorithm, taking into account the
sparse structure of the S matrix

2) The Householder matrix V has the special form:

      N
V = [ I  ] N
    [ V~ ] M

The matrix V~ is stored in A on output; the identity is not stored

3) The orthogonal matrix is

Q = I - V T V^T
*/

int
gsl_linalg_QR_UR_decomp (gsl_matrix * S, gsl_matrix * A, gsl_matrix * T)
{
  const size_t M = A->size1;
  const size_t N = S->size1;

  if (N != S->size2)
    {
      GSL_ERROR ("S matrix must be square", GSL_ENOTSQR);
    }
  else if (N != A->size2)
    {
      GSL_ERROR ("S and A have different number of columns", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* base case, compute Householder transform for single column matrix */
      double * T00 = gsl_matrix_ptr(T, 0, 0);
      double * S00 = gsl_matrix_ptr(S, 0, 0);
      gsl_vector_view v = gsl_matrix_column(A, 0);
      *T00 = qrtr_householder_transform(S00, &v.vector);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * partition matrices:
       *
       *       N1  N2              N1  N2
       * N1 [ S11 S12 ] and  N1 [ T11 T12 ]
       * N2 [  0  S22 ]      N2 [  0  T22 ]
       * M  [  A1  A2 ]
       */
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view S11 = gsl_matrix_submatrix(S, 0, 0, N1, N1);
      gsl_matrix_view S12 = gsl_matrix_submatrix(S, 0, N1, N1, N2);
      gsl_matrix_view S22 = gsl_matrix_submatrix(S, N1, N1, N2, N2);

      gsl_matrix_view A1 = gsl_matrix_submatrix(A, 0, 0, M, N1);
      gsl_matrix_view A2 = gsl_matrix_submatrix(A, 0, N1, M, N2);

      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      /*
       * Eq. 2: recursively factor
       *
       *       N1           N1
       * N1 [ S11 ] = Q1 [ R11 ] N1
       * N2 [  0  ]      [  0  ] N2
       * M  [  A1 ]      [  0  ] M
       */
      status = gsl_linalg_QR_UR_decomp(&S11.matrix, &A1.matrix, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       *      N2              N2            N2
       * N1 [ R12  ] = Q1^T [ S12 ] = [   S12 - W  ] N1
       * N2 [ S22~ ]        [ S22 ]   [     S22    ] N2
       * M  [  A2~ ]        [  A2 ]   [ A2 - V1~ W ] M
       *
       * where W = T11^T ( S12 + V1~^T A2 ), using T12 as temporary storage, and
       *
       *        N1
       * V1 = [  I  ] N1
       *      [  0  ] N2
       *      [ V1~ ] M
       */
      gsl_matrix_memcpy(&T12.matrix, &S12.matrix);                                                    /* W := S12 */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);        /* W := S12 + V1~^T A2 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix); /* W := T11^T W */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A1.matrix, &T12.matrix, 1.0, &A2.matrix);     /* A2 := A2 - V1~ W */
      gsl_matrix_sub(&S12.matrix, &T12.matrix);                                                       /* R12 := S12 - W */

      /*
       * Eq. 4: recursively factor
       *
       * [ S22~ ] = Q2~ [ R22 ]
       * [  A2~ ]       [  0  ]
       */
      status = gsl_linalg_QR_UR_decomp(&S22.matrix, &A2.matrix, &T22.matrix);
      if (status)
        return status;

      /*
       * Eq. 13: update T12 := -T11 * V1^T * V2 * T22
       *
       * where:
       *
       *        N1                N2
       * V1 = [  I  ] N1   V2 = [  0  ] N1
       *      [  0  ] N2        [  I  ] N2
       *      [ V1~ ] M         [ V2~ ] M
       *
       * Note: V1^T V2 = V1~^T V2~
       */

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 0.0, &T12.matrix);           /* T12 := V1~^T * V2~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 := -T11 * T12 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   [ U ] x = b
 *   [ A ]
 *  
 * using the QR factorization [ U; A ] = Q R. 
 *
 * Inputs: R    - upper triangular R matrix, N-by-N
 *         Y    - dense Y matrix, M-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         b    - right hand side, size N+M
 *         x    - (output) solution, size N+M
 *                x(1:N) = least squares solution vector
 *                x(N+1:N+M) = vector whose norm equals ||b - [U; A] x||
 *         work - workspace, size N
 */

int
gsl_linalg_QR_UR_lssolve (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T,
                          const gsl_vector * b, gsl_vector * x, gsl_vector * work)
{
  const size_t N = R->size1;
  const size_t M = Y->size1;

  if (R->size2 != N)
    {
      GSL_ERROR ("R matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size2 != N)
    {
      GSL_ERROR ("Y matrix must have N columns", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (N+M != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (N+M != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (N != work->size)
    {
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      gsl_vector_memcpy(x, b);
      status = gsl_linalg_QR_UR_lssvx(R, Y, T, x, work);

      return status;
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   [ U ] x = b
 *   [ A ]
 *  
 * using the QR factorization [ U; A ] = Q R. 
 *
 * Inputs: R    - upper triangular R matrix, N-by-N
 *         Y    - dense Y matrix, M-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         x    - (input/output) solution, size N+M
 *                on input, right hand side vector b, length N+M
 *                on output,
 *                  x(1:N) = least squares solution vector
 *                  x(N+1:N+M) = vector whose norm equals ||b - [U; A] x||
 *         work - workspace, size N
 */

int
gsl_linalg_QR_UR_lssvx (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T,
                        gsl_vector * x, gsl_vector * work)
{
  const size_t N = R->size1;
  const size_t M = Y->size1;

  if (R->size2 != N)
    {
      GSL_ERROR ("R matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size2 != N)
    {
      GSL_ERROR ("Y matrix must have N columns", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (N+M != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (N != work->size)
    {
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_view x1 = gsl_vector_subvector(x, 0, N);

      /* compute x = Q^T b */
      gsl_linalg_QR_UR_QTvec (Y, T, x, work);

      /* Solve R x = Q^T b */
      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, R, &x1.vector);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_QR_UR_QTvec()
  Apply 2N-by-2N Q^T to the 2N-by-1 vector b

Inputs: Y    - upper triangular Y matrix encoded by gsl_linalg_QR_UR_decomp, N-by-N
        T    - block reflector matrix, N-by-N
        b    - 2N-by-1 vector replaced by Q^T b on output
        work - workspace, length N

Notes:
1) Q^T b = (I - V T^T V^T) b
         = b - V T^T [ I Y^T ] [ b1 ]
                                   [ b2 ]
         = b - V T^T [ b1 + Y^T b2 ]
         = [ b1 ] - [  w  ]
           [ b2 ]   [ Y w ]

where w = T^T ( b1 + Y^T b2 )
*/

int
gsl_linalg_QR_UR_QTvec(const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * b, gsl_vector * work)
{
  const size_t M = Y->size1;
  const size_t N = Y->size2;

  if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (b->size != M + N)
    {
      GSL_ERROR ("b vector must have length M+N", GSL_EBADLEN);
    }
  else if (work->size != N)
    {
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_view b1 = gsl_vector_subvector(b, 0, N);
      gsl_vector_view b2 = gsl_vector_subvector(b, N, M);

      /* work := YT b2 */
      gsl_blas_dgemv(CblasTrans, 1.0, Y, &b2.vector, 0.0, work);

      /* work = b1 + YT b2 */
      gsl_vector_add(work, &b1.vector);

      /* work = T^T * work */
      gsl_blas_dtrmv(CblasUpper, CblasTrans, CblasNonUnit, T, work);

      /* b1 := b1 - work */
      gsl_vector_sub(&b1.vector, work);

      /* b2 := b2 - Y w */
      gsl_blas_dgemv(CblasNoTrans, -1.0, Y, work, 1.0, &b2.vector);

      return GSL_SUCCESS;
    }
}

/*
qrtr_householder_transform()
  This routine is an optimized version of
gsl_linalg_householder_transform(), designed for the QR
decomposition of M-by-N matrices of the form:

B = [ S ]
    [ A ]

where S is N-by-N upper triangular, and A is M-by-N dense.
This routine computes a householder transformation (tau,v) of a 
x so that P x = [ I - tau*v*v' ] x annihilates x(1:n-1). x will
be a subcolumn of the matrix B, and so its structure will be:

x = [ x0 ] <- 1 nonzero value for the diagonal element of S
    [ 0  ] <- N - j - 1 zeros, where j is column of matrix in [0,N-1]
    [ x  ] <- M nonzero values for the dense part A

Inputs: v0 - pointer to diagonal element of S
             on input, v0 = x0;
        v  - on input, x vector
             on output, householder vector v
*/

static double
qrtr_householder_transform (double *v0, gsl_vector * v)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  double alpha, beta, tau ;
  
  /* compute xnorm = || [ 0 ; v ] ||, ignoring zero part of vector */
  double xnorm = gsl_blas_dnrm2(v);

  if (xnorm == 0) 
    {
      return 0.0; /* tau = 0 */
    }

  alpha = *v0;
  beta = - GSL_SIGN(alpha) * hypot(alpha, xnorm) ;
  tau = (beta - alpha) / beta ;
  
  {
    double s = (alpha - beta);
    
    if (fabs(s) > GSL_DBL_MIN) 
      {
        gsl_blas_dscal (1.0 / s, v);
        *v0 = beta;
      }
    else
      {
        gsl_blas_dscal (GSL_DBL_EPSILON / s, v);
        gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, v);
        *v0 = beta;
      }
  }
  
  return tau;
}
