/* linalg/qr_ud.c
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

/*
 * this module contains routines for the QR factorization of a matrix
 * using the recursive Level 3 BLAS algorithm of Elmroth and Gustavson with
 * additional modifications courtesy of Julien Langou.
 */

static int aux_QR_TRD_decomp (gsl_matrix * U, gsl_matrix * A, const gsl_vector * D, gsl_matrix * T);
static double qrtd_householder_transform (double *v0, double *v1);
static double qrtrd_householder_transform (double *v0, gsl_vector * v, double *d);

/*
gsl_linalg_QR_UD_decomp()
  Compute the QR decomposition of the "triangle on top of diagonal" matrix

  [ U ] = Q [ R ]
  [ D ]     [ 0 ]

where U is N-by-N upper triangular, D is N-by-N diagonal

Inputs: U - on input, upper triangular N-by-N matrix
            on output, R factor in upper triangle
        D - diagonal N-by-N matrix
        Y - (output) upper triangular N-by-N matrix (lower part of V)
        T - (output) block reflector matrix, N-by-N

Notes:
1) Based on the Elmroth/Gustavson algorithm, taking into account the
sparse structure of the U,S matrices

2) The Householder matrix V has the special form:

      N
V = [ I ] N
    [ Y ] N

with Y upper triangular

3) The orthogonal matrix is

Q = I - V T V^T
*/

int
gsl_linalg_QR_UD_decomp (gsl_matrix * U, const gsl_vector * D, gsl_matrix * Y, gsl_matrix * T)
{
  const size_t N = U->size1;

  if (N != U->size2)
    {
      GSL_ERROR ("U matrix must be square", GSL_ENOTSQR);
    }
  else if (D->size != N)
    {
      GSL_ERROR ("D matrix does not match dimensions of U", GSL_EBADLEN);
    }
  else if (Y->size1 != Y->size2)
    {
      GSL_ERROR ("Y matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size1 != N)
    {
      GSL_ERROR ("Y matrix does not match dimensions of U", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* base case, compute Householder transform for single column matrix */
      double * T00 = gsl_matrix_ptr(T, 0, 0);
      double * U00 = gsl_matrix_ptr(U, 0, 0);
      double * Y00 = gsl_matrix_ptr(Y, 0, 0);
      *Y00 = gsl_vector_get(D, 0);
      *T00 = qrtd_householder_transform(U00, Y00);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * partition matrices:
       *
       *       N1  N2              N1  N2
       * N1 [ U11 U12 ] and  N1 [ T11 T12 ]
       * N2 [  0  U22 ]      N2 [  0  T22 ]
       * N1 [ D11 D12 ]
       * N2 [  0  D22 ]
       */
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view U11 = gsl_matrix_submatrix(U, 0, 0, N1, N1);
      gsl_matrix_view U12 = gsl_matrix_submatrix(U, 0, N1, N1, N2);
      gsl_matrix_view U22 = gsl_matrix_submatrix(U, N1, N1, N2, N2);

      gsl_vector_const_view D1 = gsl_vector_const_subvector(D, 0, N1);
      gsl_vector_const_view D2 = gsl_vector_const_subvector(D, N1, N2);

      gsl_matrix_view Y11 = gsl_matrix_submatrix(Y, 0, 0, N1, N1);
      gsl_matrix_view Y12 = gsl_matrix_submatrix(Y, 0, N1, N1, N2);

      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      gsl_matrix_view m;

      /*
       * Eq. 2: recursively factor
       *
       *       N1           N1
       * N1 [ U11 ] = Q1 [ R11 ] N1
       * N2 [  0  ]      [  0  ] N2
       * N1 [ D11 ]      [  0  ] N1
       * N2 [  0  ]      [  0  ] N2
       */
      status = gsl_linalg_QR_UD_decomp(&U11.matrix, &D1.vector, &Y11.matrix, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       *      N2              N2            N2
       * N1 [ R12  ] = Q1^T [ U12 ] = [   U12 - W   ] N1
       * N2 [ U22~ ]        [ U22 ]   [     U22     ] N2
       * N1 [ D12~ ]        [  0  ]   [  - V31 W    ] N1
       * N2 [ D22~ ]        [ D22 ]   [     D22     ] N2
       *
       * where W = T11^T U12, using T12 as temporary storage, and
       *
       *        N1
       * V1 = [  I  ] N1
       *      [  0  ] N2
       *      [ V31 ] N1
       *      [  0  ] N2
       */
      gsl_matrix_memcpy(&T12.matrix, &U12.matrix);                                                    /* W := U12 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix); /* W := T11^T W */
      gsl_matrix_sub(&U12.matrix, &T12.matrix);                                                       /* R12 := U12 - W */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &Y11.matrix, &T12.matrix); /* W := -V31 W */
      gsl_matrix_memcpy(&Y12.matrix, &T12.matrix);                                                    /* Y12 := -V31 W */

      /*
       * Eq. 4: factor the "triangle on top of rectangle on top of diagonal" matrix
       *
       * N2 [ U22 ] = Q2~ [ R22 ] N2
       * N1 [ Y12 ]       [  0  ] N1
       * N2 [ D22 ]       [  0  ] N2
       */
      m = gsl_matrix_submatrix(Y, 0, N1, N, N2);
      status = aux_QR_TRD_decomp(&U22.matrix, &m.matrix, &D2.vector, &T22.matrix);
      if (status)
        return status;

      /*
       * Eq. 13: update T12 := -T11 * V1^T * V2 * T22
       *
       * where:
       *
       *        N1                N2
       * V1 = [  I   ] N1   V2 = [  0   ] N1
       *      [  0   ] N2        [  I   ] N2
       *      [ V31~ ] N1        [ V32~ ] N1
       *      [  0   ] N2        [ V42~ ] N2
       *
       * Note: V1^T V2 = V31~^T V32~
       */

      gsl_matrix_memcpy(&T12.matrix, &Y12.matrix);                                                    /* T12 := V32~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &Y11.matrix, &T12.matrix); /* T12 := V31~^T V32~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 := -T11 * T12 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   [ U ] x = b
 *   [ D ]
 *  
 * using the QR factorization [ U; D ] = Q R. 
 *
 * Inputs: R    - upper triangular R matrix, N-by-N
 *         Y    - upper triangular Y matrix, N-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         b    - right hand side, size 2*N
 *         x    - (output) solution, size 2*N
 *                x(1:N) = least squares solution vector
 *                x(N+1:2*N) = vector whose norm equals ||b - Ax||
 *         work - workspace, size N
 */

int
gsl_linalg_QR_UD_lssolve (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T,
                          const gsl_vector * b, gsl_vector * x, gsl_vector * work)
{
  const size_t N = R->size1;
  const size_t M = 2 * N;

  if (R->size2 != N)
    {
      GSL_ERROR ("R matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size1 != Y->size2)
    {
      GSL_ERROR ("Y matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size1 != N)
    {
      GSL_ERROR ("Y and R must have same dimensions", GSL_EBADLEN);
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
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      return gsl_linalg_QR_UU_lssolve(R, Y, T, b, x, work);
    }
}

/* Find the least squares solution to the overdetermined system 
 *
 *   [ U ] x = b
 *   [ D ]
 *  
 * using the QR factorization [ U; D ] = Q R. 
 *
 * Inputs: R    - upper triangular R matrix, N-by-N
 *         Y    - upper triangular Y matrix, N-by-N
 *         T    - upper triangular block reflector, N-by-N
 *         x    - (input/output) solution, size 2*N
 *                on input, right hand side vector b, length 2*N;
 *                on output,
 *                  x(1:N) = least squares solution vector
 *                  x(N+1:2*N) = vector whose norm equals ||b - Ax||
 *         work - workspace, size N
 */

int
gsl_linalg_QR_UD_lssvx (const gsl_matrix * R, const gsl_matrix * Y, const gsl_matrix * T,
                        gsl_vector * x, gsl_vector * work)
{
  const size_t N = R->size1;
  const size_t M = 2 * N;

  if (R->size2 != N)
    {
      GSL_ERROR ("R matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size1 != Y->size2)
    {
      GSL_ERROR ("Y matrix must be square", GSL_ENOTSQR);
    }
  else if (Y->size1 != N)
    {
      GSL_ERROR ("Y and R must have same dimensions", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (M != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else if (N != work->size)
    {
      GSL_ERROR ("workspace must be length N", GSL_EBADLEN);
    }
  else
    {
      return gsl_linalg_QR_UU_lssvx(R, Y, T, x, work);
    }
}

/*
gsl_linalg_QR_UD_QTvec()
  Apply 2N-by-2N Q^T to the 2N-by-1 vector b

Inputs: Y    - upper triangular Y matrix encoded by gsl_linalg_QR_UD_decomp, N-by-N
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
gsl_linalg_QR_UD_QTvec(const gsl_matrix * Y, const gsl_matrix * T, gsl_vector * b, gsl_vector * work)
{
  const size_t N = Y->size1;
  const size_t M = 2 * N;

  if (Y->size2 != N)
    {
      GSL_ERROR ("Y matrix must be square", GSL_ENOTSQR);
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
      return gsl_linalg_QR_UU_QTvec(Y, T, b, work);
    }
}

/*
aux_QR_TRD_decomp()
  Compute the QR decomposition of the "triangle on top of rectangle on top of diagonal" matrix

  [ U ] = Q [ R ]
  [ A ]     [ 0 ]
  [ D ]     [ 0 ]

where U is N-by-N upper triangular, A is M-by-N dense, D is N-by-N diagonal

Inputs: U - on input, upper triangular N-by-N matrix
            on output, R factor in upper triangle
        A - on input, A is stored in A(1:M,1:N)
            on output, Y is stored in A(1:M+N,1:N)
            size (M+N)-by-N
        D - diagonal N-by-N matrix
        T - (output) block reflector matrix, N-by-N

Notes:
1) Based on the Elmroth/Gustavson algorithm, taking into account the
sparse structure of the U,S matrices

2) The Householder matrix V has the special form:

      N
V = [ I ] N
    [ Y ] M+N

with Y upper trapezoidal

3) The orthogonal matrix is

Q = I - V T V^T
*/

static int
aux_QR_TRD_decomp (gsl_matrix * U, gsl_matrix * A, const gsl_vector * D, gsl_matrix * T)
{
  const size_t N = U->size1;

  if (N != U->size2)
    {
      GSL_ERROR ("U matrix must be square", GSL_ENOTSQR);
    }
  else if (A->size1 <= N)
    {
      GSL_ERROR ("A matrix must have > N rows", GSL_EBADLEN);
    }
  else if (D->size != N)
    {
      GSL_ERROR ("D matrix does not match dimensions of U", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* base case, compute Householder transform for single column matrix */
      const size_t M = A->size1 - N;
      double * T00 = gsl_matrix_ptr(T, 0, 0);
      double * U00 = gsl_matrix_ptr(U, 0, 0);
      gsl_vector_view v = gsl_matrix_subcolumn(A, 0, 0, M);
      double * D00 = gsl_matrix_ptr(A, M, 0);
      *D00 = gsl_vector_get(D, 0);
      *T00 = qrtrd_householder_transform(U00, &v.vector, D00);
      return GSL_SUCCESS;
    }
  else
    {
      /*
       * partition matrices:
       *
       *       N1  N2              N1  N2
       * N1 [ U11 U12 ] and  N1 [ T11 T12 ]
       * N2 [  0  U22 ]      N2 [  0  T22 ]
       * N1 [ D11 D12 ]
       * N2 [  0  D22 ]
       */
      int status;
      const size_t M = A->size1 - N;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view U11 = gsl_matrix_submatrix(U, 0, 0, N1, N1);
      gsl_matrix_view U12 = gsl_matrix_submatrix(U, 0, N1, N1, N2);
      gsl_matrix_view U22 = gsl_matrix_submatrix(U, N1, N1, N2, N2);

      gsl_matrix_view A1 = gsl_matrix_submatrix(A, 0, 0, M, N1);
      gsl_matrix_view A2 = gsl_matrix_submatrix(A, 0, N1, M, N2);

      gsl_vector_const_view D1 = gsl_vector_const_subvector(D, 0, N1);
      gsl_vector_const_view D2 = gsl_vector_const_subvector(D, N1, N2);

      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      gsl_matrix_view Y41 = gsl_matrix_submatrix(A, M, 0, N1, N1);
      gsl_matrix_view Y42 = gsl_matrix_submatrix(A, M, N1, N1, N2);

      gsl_matrix_view m;

      /*
       * Eq. 2: recursively factor
       *
       *       N1           N1
       * N1 [ U11 ] = Q1 [ R11 ] N1
       * N2 [  0  ]      [  0  ] N2
       * M  [ A1  ]      [  0  ] M
       * N1 [ D11 ]      [  0  ] N1
       * N2 [  0  ]      [  0  ] N2
       */
      m = gsl_matrix_submatrix(A, 0, 0, M + N1, N1);
      status = aux_QR_TRD_decomp(&U11.matrix, &m.matrix, &D1.vector, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       *      N2              N2            N2
       * N1 [ R12  ] = Q1^T [ U12 ] = [   U12 - W   ] N1
       * N2 [ U22~ ]        [ U22 ]   [     U22     ] N2
       * M  [ A2~  ]        [ A2  ]   [ A2 - V31 W  ] M
       * N1 [ D12~ ]        [  0  ]   [  - V41 W    ] N1
       * N2 [ D22~ ]        [ D22 ]   [     D22     ] N2
       *
       * where W = T11^T ( U12 + V31^T A2 ), using T12 as temporary storage, and
       *
       *        N1
       * V1 = [  I  ] N1
       *      [  0  ] N2
       *      [ V31 ] M
       *      [ V41 ] N1
       *      [  0  ] N2
       */
      gsl_matrix_memcpy(&T12.matrix, &U12.matrix);                                                    /* W := U12 */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);        /* W := U12 + V31^T A2 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix); /* W := T11^T W */
      gsl_matrix_sub(&U12.matrix, &T12.matrix);                                                       /* R12 := U12 - W */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A1.matrix, &T12.matrix, 1.0, &A2.matrix);     /* A2 := A2 - V31 W */

      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &Y41.matrix, &T12.matrix); /* W := -V41 W */
      gsl_matrix_memcpy(&Y42.matrix, &T12.matrix);                                                       /* V42 := -431 W */

      /*
       * Eq. 4: factor the "triangle on top of rectangle on top of diagonal" matrix
       *
       * N2 [ U22 ] = Q2~ [ R22 ] N2
       * N1 [ Y12 ]       [  0  ] N1
       * N2 [ D22 ]       [  0  ] N2
       */
      m = gsl_matrix_submatrix(A, 0, N1, M + N, N2);
      status = aux_QR_TRD_decomp(&U22.matrix, &m.matrix, &D2.vector, &T22.matrix);
      if (status)
        return status;

      /*
       * Eq. 13: update T12 := -T11 * V1^T * V2 * T22
       *
       * where:
       *
       *        N1                N2
       * V1 = [  I   ] N1   V2 = [  0   ] N1
       *      [  0   ] N2        [  I   ] N2
       *      [ V31~ ] M         [ V32~ ] M 
       *      [ V41~ ] N1        [ V42~ ] N1
       *      [  0   ] N2        [ V52~ ] N2
       *
       * Note: V1^T V2 = V31~^T V32~ + V41~^T V42~
       */

      gsl_matrix_memcpy(&T12.matrix, &Y42.matrix);                                                    /* T12 := V42~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &Y41.matrix, &T12.matrix); /* T12 := V41~^T V42~ */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);        /* T12 := V31~^T V32~^T + V41~^T V42~ */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 := -T11 * T12 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
}

/*
qrtd_householder_transform()
  This routine is an optimized version of
gsl_linalg_householder_transform(), designed for the QR
decomposition of M-by-N matrices of the form:

B = [ U ]
    [ S ]

where U,S are N-by-N upper triangular.
This routine computes a householder transformation (tau,v) of a 
x so that P x = [ I - tau*v*v' ] x annihilates x(1:n-1). x will
be a subcolumn of the matrix B, and so its structure will be:

x = [ x0 ] <- 1 nonzero value for the diagonal element of U
    [ y0 ] <- 1 nonzero value for the diagonal element of S

Inputs: v0 - pointer to diagonal element of U
             on input, v0 = x0;
        v1 - on input, pointer to diagonal element of S
             on output, householder vector v
*/

static double
qrtd_householder_transform (double *v0, double *v1)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  double alpha, beta, tau ;
  
  /* compute xnorm = || [ 0 ; v ] ||, ignoring zero part of vector */
  double xnorm = fabs(*v1);

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
        *v1 /= s;
        *v0 = beta;
      }
    else
      {
        *v1 *= GSL_DBL_EPSILON / s;
        *v1 /= GSL_DBL_EPSILON;
        *v0 = beta;
      }
  }
  
  return tau;
}

/*
qrtrd_householder_transform()
  This routine is an optimized version of
gsl_linalg_householder_transform(), designed for the QR
decomposition of M-by-N matrices of the form:

B = [ U ]
    [ A ]
    [ D ]

where U is N-by-N upper triangular A is M-by-N dense and D is N-by-N diagonal.
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
        d  - on input, diagonal element of D
             on output, householder element
*/

static double
qrtrd_householder_transform (double *v0, gsl_vector * v, double *d)
{
  /* replace v[0:M-1] with a householder vector (v[0:M-1]) and
     coefficient tau that annihilate v[1:M-1] */

  double alpha, beta, tau ;
  
  /* compute xnorm = || [ 0 ; v ; 0 ; d] ||, ignoring zero part of vector */
  double xnorm;
  
  xnorm = gsl_blas_dnrm2(v);
  xnorm = gsl_hypot(xnorm, *d);

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
        *d /= s;
        *v0 = beta;
      }
    else
      {
        gsl_blas_dscal (GSL_DBL_EPSILON / s, v);
        gsl_blas_dscal (1.0 / GSL_DBL_EPSILON, v);
        *d *= GSL_DBL_EPSILON / s;
        *d /= GSL_DBL_EPSILON;
        *v0 = beta;
      }
  }
  
  return tau;
}
