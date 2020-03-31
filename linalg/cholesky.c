/* Cholesky Decomposition
 *
 * Copyright (C) 2000 Thomas Walter
 * Copyright (C) 2000, 2001, 2002, 2003, 2005, 2007 Brian Gough, Gerard Jungman
 * Copyright (C) 2016, 2019 Patrick Alken
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * 03 May 2000: Modified for GSL by Brian Gough
 * 29 Jul 2005: Additions by Gerard Jungman
 * 04 Mar 2016: Change Cholesky algorithm to gaxpy version by Patrick Alken
 * 26 May 2019: implement recursive Cholesky with Level 3 BLAS by Patrick Alken
 */

/*
 * Cholesky decomposition of a symmetric positive definite matrix.
 *
 * This algorithm does:
 *   A = L * L'
 * with
 *   L  := lower left triangle matrix
 *   L' := the transposed form of L.
 *
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "recurse.h"

static double cholesky_norm1(const gsl_matrix * LLT, gsl_vector * work);
static int cholesky_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);
static int cholesky_decomp_L2 (gsl_matrix * A);
static int cholesky_decomp_L3 (gsl_matrix * A);

/*
In GSL 2.2, we decided to modify the behavior of the Cholesky decomposition
to store the Cholesky factor in the lower triangle, and store the original
matrix in the upper triangle. Previous versions stored the Cholesky factor in
both places. The routine gsl_linalg_cholesky_decomp1 was added for the new
behavior, and gsl_linalg_cholesky_decomp is maintained for backward compatibility.
It will be removed in a future release.
*/

int
gsl_linalg_cholesky_decomp (gsl_matrix * A)
{
  int status;

  status = gsl_linalg_cholesky_decomp1(A);
  if (status == GSL_SUCCESS)
    {
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, A, A);
    }

  return status;
}

/*
gsl_linalg_cholesky_decomp1()
  Perform Cholesky decomposition of a symmetric positive
definite matrix using lower triangle using Level 3 BLAS algorithm.

Inputs: A - (input) symmetric, positive definite matrix
            (output) lower triangle contains Cholesky factor

Return: success/error

Notes:
1) original matrix is saved in upper triangle on output
*/

int
gsl_linalg_cholesky_decomp1 (gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("Cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      /* save original matrix in upper triangle for later rcond calculation */
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, A, A);

      return cholesky_decomp_L3(A);
    }
}

int
gsl_linalg_cholesky_solve (const gsl_matrix * LLT,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_cholesky_svx(LLT, x);

      return status;
    }
}

int
gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                         gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* solve for c using forward-substitution, L c = b */
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* perform back-substitution, L^T x = c */
      gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit, LLT, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_solve_mat (const gsl_matrix * LLT,
                               const gsl_matrix * B,
                               gsl_matrix * X)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != B->size1)
    {
      GSL_ERROR ("matrix size must match B size", GSL_EBADLEN);
    }
  else if (LLT->size2 != X->size1)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy X <- B */
      gsl_matrix_memcpy (X, B);

      status = gsl_linalg_cholesky_svx_mat(LLT, X);

      return status;
    }
}

int
gsl_linalg_cholesky_svx_mat (const gsl_matrix * LLT,
                             gsl_matrix * X)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != X->size1)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* solve for C using forward-substitution, L C = B */
      gsl_blas_dtrsm (CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0,
                      LLT, X);

      /* perform back-substitution, L^T X = C */
      gsl_blas_dtrsm (CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0,
                      LLT, X);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_cholesky_invert()
  Compute the inverse of a symmetric positive definite matrix in
Cholesky form.

Inputs: LLT - matrix in cholesky form on input
              A^{-1} = L^{-t} L^{-1} on output

Return: success or error
*/

int
gsl_linalg_cholesky_invert(gsl_matrix * LLT)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      int status;

      /* invert the lower triangle of LLT */
      status = gsl_linalg_tri_invert(CblasLower, CblasNonUnit, LLT);
      if (status)
        return status;

      /* compute A^{-1} = L^{-T} L^{-1} */
      status = gsl_linalg_tri_LTL(LLT);
      if (status)
        return status;

      /* copy lower triangle to upper */
      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, LLT, LLT);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_decomp_unit(gsl_matrix * A, gsl_vector * D)
{
  const size_t N = A->size1;
  size_t i, j;

  /* initial Cholesky */
  int stat_chol = gsl_linalg_cholesky_decomp1(A);

  if(stat_chol == GSL_SUCCESS)
  {
    /* calculate D from diagonal part of initial Cholesky */
    for(i = 0; i < N; ++i)
    {
      const double C_ii = gsl_matrix_get(A, i, i);
      gsl_vector_set(D, i, C_ii*C_ii);
    }

    /* multiply initial Cholesky by 1/sqrt(D) on the right */
    for(i = 0; i < N; ++i)
    {
      for(j = 0; j < N; ++j)
      {
        gsl_matrix_set(A, i, j, gsl_matrix_get(A, i, j) / sqrt(gsl_vector_get(D, j)));
      }
    }

    /* Because the initial Cholesky contained both L and transpose(L),
       the result of the multiplication is not symmetric anymore;
       but the lower triangle _is_ correct. Therefore we reflect
       it to the upper triangle and declare victory.
       */
    for(i = 0; i < N; ++i)
      for(j = i + 1; j < N; ++j)
        gsl_matrix_set(A, i, j, gsl_matrix_get(A, j, i));
  }

  return stat_chol;
}

/*
gsl_linalg_cholesky_scale()
  This function computes scale factors diag(S), such that

diag(S) A diag(S)

has a condition number within a factor N of the matrix
with the smallest condition number over all possible
diagonal scalings. See Corollary 7.6 of:

N. J. Higham, Accuracy and Stability of Numerical Algorithms (2nd Edition),
SIAM, 2002.

Inputs: A - symmetric positive definite matrix
        S - (output) scale factors, S_i = 1 / sqrt(A_ii)
*/

int
gsl_linalg_cholesky_scale(const gsl_matrix * A, gsl_vector * S)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("A is not a square matrix", GSL_ENOTSQR);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      /* compute S_i = 1/sqrt(A_{ii}) */
      for (i = 0; i < N; ++i)
        {
          double Aii = gsl_matrix_get(A, i, i);

          if (Aii <= 0.0)
            gsl_vector_set(S, i, 1.0); /* matrix not positive definite */
          else
            gsl_vector_set(S, i, 1.0 / sqrt(Aii));
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_cholesky_scale_apply()
  This function applies scale transformation to A:

A <- diag(S) A diag(S)

Inputs: A     - (input/output)
                on input, symmetric positive definite matrix
                on output, diag(S) * A * diag(S) in lower triangle
        S     - (input) scale factors
*/

int
gsl_linalg_cholesky_scale_apply(gsl_matrix * A, const gsl_vector * S)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("A is not a square matrix", GSL_ENOTSQR);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      size_t i, j;

      /* compute: A <- diag(S) A diag(S) using lower triangle */
      for (j = 0; j < N; ++j)
        {
          double sj = gsl_vector_get(S, j);

          for (i = j; i < N; ++i)
            {
              double si = gsl_vector_get(S, i);
              double *Aij = gsl_matrix_ptr(A, i, j);
              *Aij *= si * sj;
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_decomp2(gsl_matrix * A, gsl_vector * S)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (N != S->size)
    {
      GSL_ERROR("S must have length N", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* compute scaling factors to reduce cond(A) */
      status = gsl_linalg_cholesky_scale(A, S);
      if (status)
        return status;

      /* apply scaling factors */
      status = gsl_linalg_cholesky_scale_apply(A, S);
      if (status)
        return status;

      /* compute Cholesky decomposition of diag(S) A diag(S) */
      status = gsl_linalg_cholesky_decomp1(A);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_svx2 (const gsl_matrix * LLT,
                          const gsl_vector * S,
                          gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != S->size)
    {
      GSL_ERROR ("matrix size must match S", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* b~ = diag(S) b */
      gsl_vector_mul(x, S);

      /* Solve for c using forward-substitution, L c = b~ */
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, L^T x~ = c */
      gsl_blas_dtrsv (CblasLower, CblasTrans, CblasNonUnit, LLT, x);

      /* compute original solution vector x = S x~ */
      gsl_vector_mul(x, S);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_solve2 (const gsl_matrix * LLT,
                            const gsl_vector * S,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != S->size)
    {
      GSL_ERROR ("matrix size must match S size", GSL_EBADLEN);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* Copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_cholesky_svx2(LLT, S, x);

      return status;
    }
}

int
gsl_linalg_cholesky_rcond (const gsl_matrix * LLT, double * rcond,
                           gsl_vector * work)
{
  const size_t M = LLT->size1;
  const size_t N = LLT->size2;

  if (M != N)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      double Anorm = cholesky_norm1(LLT, work); /* ||A||_1 */
      double Ainvnorm;                          /* ||A^{-1}||_1 */

      *rcond = 0.0;

      /* don't continue if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      /* estimate ||A^{-1}||_1 */
      status = gsl_linalg_invnorm1(N, cholesky_Ainv, (void *) LLT, &Ainvnorm, work);

      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/* compute 1-norm of original matrix, stored in upper triangle of LLT;
 * diagonal entries have to be reconstructed */
static double
cholesky_norm1(const gsl_matrix * LLT, gsl_vector * work)
{
  const size_t N = LLT->size1;
  double max = 0.0;
  size_t i, j;

  for (j = 0; j < N; ++j)
    {
      double sum = 0.0;
      gsl_vector_const_view lj = gsl_matrix_const_subrow(LLT, j, 0, j + 1);
      double Ajj;

      /* compute diagonal (j,j) entry of A */
      gsl_blas_ddot(&lj.vector, &lj.vector, &Ajj);

      for (i = 0; i < j; ++i)
        {
          double *wi = gsl_vector_ptr(work, i);
          double Aij = gsl_matrix_get(LLT, i, j);
          double absAij = fabs(Aij);

          sum += absAij;
          *wi += absAij;
        }

      gsl_vector_set(work, j, sum + fabs(Ajj));
    }

  for (i = 0; i < N; ++i)
    {
      double wi = gsl_vector_get(work, i);
      max = GSL_MAX(max, wi);
    }

  return max;
}

/* x := A^{-1} x = A^{-t} x, A = L L^T */
static int
cholesky_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  int status;
  gsl_matrix * A = (gsl_matrix * ) params;

  (void) TransA; /* unused parameter warning */

  /* compute L^{-1} x */
  status = gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, A, x);
  if (status)
    return status;

  /* compute L^{-t} x */
  status = gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, A, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
cholesky_decomp_L2()
  Perform Cholesky decomposition of a symmetric positive
definite matrix using lower triangle

Inputs: A - (input) symmetric, positive definite matrix
            (output) lower triangle contains Cholesky factor

Return: success/error

Notes:
1) Based on algorithm 4.2.1 (Gaxpy Cholesky) of Golub and
Van Loan, Matrix Computations (4th ed), using Level 2 BLAS.
*/

static int
cholesky_decomp_L2 (gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("Cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t j;

      for (j = 0; j < N; ++j)
        {
          double ajj;
          gsl_vector_view v = gsl_matrix_subcolumn(A, j, j, N - j); /* A(j:n,j) */

          if (j > 0)
            {
              gsl_vector_view w = gsl_matrix_subrow(A, j, 0, j);           /* A(j,1:j-1)^T */
              gsl_matrix_view m = gsl_matrix_submatrix(A, j, 0, N - j, j); /* A(j:n,1:j-1) */

              gsl_blas_dgemv(CblasNoTrans, -1.0, &m.matrix, &w.vector, 1.0, &v.vector);
            }

          ajj = gsl_matrix_get(A, j, j);

          if (ajj <= 0.0)
            {
              GSL_ERROR("matrix is not positive definite", GSL_EDOM);
            }

          ajj = sqrt(ajj);
          gsl_vector_scale(&v.vector, 1.0 / ajj);
        }

      return GSL_SUCCESS;
    }
}

/*
cholesky_decomp_L3()
  Perform Cholesky decomposition of a symmetric positive
definite matrix using Level 3 BLAS.

Inputs: A - (input) symmetric, positive definite matrix in lower triangle
            (output) lower triangle contains Cholesky factor

Return: success/error

Notes:
1) Based on ReLAPACK recursive block Cholesky algorithm using Level 3 BLAS

2) 28 May 2019: performed several benchmark tests of this recursive variant
against the right-looking block variant from LAPACK. This recursive variant
performed faster in all cases, so it is now the default algorithm.
*/

static int
cholesky_decomp_L3 (gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("Cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_CHOLESKY)
    {
      /* use unblocked Level 2 algorithm */
      return cholesky_decomp_L2(A);
    }
  else
    {
      /*
       * partition matrix:
       *
       * A11 A12
       * A21 A22
       *
       * where A11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT(N);
      const size_t N2 = N - N1;
      gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, N2, N1);
      gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, N2, N2);

      /* recursion on A11 */
      status = cholesky_decomp_L3(&A11.matrix);
      if (status)
        return status;

      /* A21 = A21 * L11^{-T} */
      gsl_blas_dtrsm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, &A11.matrix, &A21.matrix);

      /* A22 -= L21 L21^T */
      gsl_blas_dsyrk(CblasLower, CblasNoTrans, -1.0, &A21.matrix, 1.0, &A22.matrix);

      /* recursion on A22 */
      status = cholesky_decomp_L3(&A22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

