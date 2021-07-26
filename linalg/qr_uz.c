/* linalg/qr_uz.c
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

static double qrtz_householder_transform (double *v0, gsl_vector * v);

/*
gsl_linalg_QR_UZ_decomp()
  Compute the QR decomposition of the "triangle on top of trapezoidal" matrix

  [ S ] = Q [ R ]
  [ A ]     [ 0 ]

where S is N-by-N upper triangular and A is M-by-N upper trapezoidal

Inputs: S   - on input, upper triangular N-by-N matrix
              on output, R factor in upper triangle
        A   - on input, M-by-N upper trapezoidal matrix
              on output, upper trapezoidal of Householder matrix V
        T   - (output) block reflector matrix, N-by-N

Notes:
1) Based on the Elmroth/Gustavson algorithm, taking into account the
sparse structure of the S,A matrices

2) The Householder matrix V has the special form:

      N
V = [ I  ] N
    [ V~ ] M

with V~ upper trapezoidal: The matrix V~ is stored in A on output;
the identity is not stored

3) The orthogonal matrix is

Q = I - V T V^T
*/

int
gsl_linalg_QR_UZ_decomp (gsl_matrix * S, gsl_matrix * A, gsl_matrix * T)
{
  const size_t M = A->size1;
  const size_t N = S->size1;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (N != S->size2)
    {
      GSL_ERROR ("S matrix must be square", GSL_ENOTSQR);
    }
  else if (N != A->size2)
    {
      GSL_ERROR ("S and A must have same number of columns", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (M == N)
    {
      /* triangle on top of triangle */
      return gsl_linalg_QR_UU_decomp(S, A, T);
    }
  else if (N == 1)
    {
      /* base case, compute Householder transform for single column matrix */
      double * T00 = gsl_matrix_ptr(T, 0, 0);
      double * S00 = gsl_matrix_ptr(S, 0, 0);
      gsl_vector_view v = gsl_matrix_column(A, 0);
      *T00 = qrtz_householder_transform(S00, &v.vector);
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
       *  M [  A1  A2 ]
       * N1 [ U11 U12 ]
       * N2 [  0  U22 ]
       */
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view S11 = gsl_matrix_submatrix(S, 0, 0, N1, N1);
      gsl_matrix_view S12 = gsl_matrix_submatrix(S, 0, N1, N1, N2);
      gsl_matrix_view S22 = gsl_matrix_submatrix(S, N1, N1, N2, N2);

      gsl_matrix_view A1 = gsl_matrix_submatrix(A, 0, 0, M - N, N1);
      gsl_matrix_view A2 = gsl_matrix_submatrix(A, 0, N1, M - N, N2);

      gsl_matrix_view U11 = gsl_matrix_submatrix(A, M - N, 0, N1, N1);
      gsl_matrix_view U12 = gsl_matrix_submatrix(A, M - N, N1, N1, N2);

      gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

      gsl_matrix_view m;

      /*
       * Eq. 2: recursively factor
       *
       *       N1           N1
       * N1  [ S11 ] = Q1 [ R11 ] N1
       * N2  [  0  ]      [  0  ] N2
       * M-N [  A1 ]      [  0  ] M-N 
       * N1  [ U11 ]      [  0  ] N1
       * N2  [  0  ]      [  0  ] N2
       */
      m = gsl_matrix_submatrix(A, 0, 0, M - N2, N1);
      status = gsl_linalg_QR_UZ_decomp(&S11.matrix, &m.matrix, &T11.matrix);
      if (status)
        return status;

      /*
       * Eq. 3:
       *
       *       N2              N2            N2
       * N1  [ R12  ] = Q1^T [ S12 ] = [   S12 - W   ] N1
       * N2  [ S22~ ]        [ S22 ]   [     S22     ] N2
       * M-N [  A2~ ]        [  A2 ]   [  A2 - V11 W ] M-N
       * N1  [ U12~ ]        [ U12 ]   [ U12 - V21 W ] N1
       * N2  [ U22~ ]        [ U22 ]   [     U22     ] N2
       *
       * where W = T11^T ( S12 + V11^T A2 + V21^T U12 ), using T12 as temporary storage, and
       *
       *        N1
       * V1 = [  I  ] N1
       *      [  0  ] N2
       *      [ V11 ] M-N
       *      [ V21 ] N1
       *      [  0  ] N2
       */
      gsl_matrix_memcpy(&T12.matrix, &U12.matrix);                                                    /* W := U12 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &U11.matrix, &T12.matrix); /* W := V21^T U12 */
      gsl_matrix_add(&T12.matrix, &S12.matrix); /* W := S12 + V21^T U12 */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);        /* W := S12 + V21^T U12 + V11^T A2 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix); /* W := T11^T W */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A1.matrix, &T12.matrix, 1.0, &A2.matrix);     /* A2 := A2 - V11 W */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &U11.matrix, &T12.matrix, 1.0, &U12.matrix);   /* U12 := U12 - V21 W */
      gsl_matrix_sub(&S12.matrix, &T12.matrix);                                                       /* R12 := S12 - W */

      /*
       * Eq. 4: recursively factor
       *
       * [ S22~ ] = Q2~ [ R22 ]
       * [  A2~ ]       [  0  ]
       * [ U12~ ]       [  0  ]
       * [ U22~ ]       [  0  ]
       */
      m = gsl_matrix_submatrix(A, 0, N1, M, N2);
      status = gsl_linalg_QR_UZ_decomp(&S22.matrix, &m.matrix, &T22.matrix);
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
       *      [ V11 ] M-N       [ V12 ] M-N
       *      [ V21 ] N1        [ V22 ] N1
       *      [  0  ] N2        [ V32 ] N2
       *
       * Note: V1^T V2 = V11^T V12 + V21^T V22
       */

      gsl_matrix_memcpy(&T12.matrix, &U12.matrix); /* T12 := V22 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &U11.matrix, &T12.matrix);    /* T12 := V21^T V22 */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A1.matrix, &A2.matrix, 1.0, &T12.matrix);           /* T12 := T12 + V11^T * V12 */
      gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 := -T11 * T12 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 := T12 * T22 */

      return GSL_SUCCESS;
    }
}

/*
qrtz_householder_transform()
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
qrtz_householder_transform (double *v0, gsl_vector * v)
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
