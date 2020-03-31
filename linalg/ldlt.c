/* linalg/ldlt.c
 * 
 * Copyright (C) 2018 Patrick Alken
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

/* L D L^T decomposition of a symmetric positive semi-definite matrix */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static double ldlt_norm1(const gsl_matrix * A);
static int ldlt_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params);

/*
gsl_linalg_ldlt_decomp()
  Perform L D L^T decomposition of a symmetric positive
semi-definite matrix using lower triangle

Inputs: A - (input) symmetric, positive semi-definite matrix
            (output) lower triangle contains L factor;
                     diagonal contains D

Return: success/error

Notes:
1) Based on algorithm 4.1.1 of Golub and Van Loan, Matrix Computations (4th ed).
2) The first subrow A(1, 2:end) is used as temporary workspace
3) The 1-norm ||A||_1 of the original matrix is stored in the upper right corner on output
*/

int
gsl_linalg_ldlt_decomp (gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR ("LDLT decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t i, j;
      double a00, anorm;
      gsl_vector_view work, v;

      /* check for quick return */
      if (N == 1)
        return GSL_SUCCESS;

      /* compute ||A||_1 */
      anorm = ldlt_norm1(A);

      /* special case first column */
      a00 = gsl_matrix_get(A, 0, 0);
      if (a00 == 0.0)
        {
          GSL_ERROR ("matrix is singular", GSL_EDOM);
        }

      v = gsl_matrix_subcolumn(A, 0, 1, N - 1);
      gsl_vector_scale(&v.vector, 1.0 / a00);

      /* use first subrow A(1, 2:end) as temporary workspace */
      work = gsl_matrix_subrow(A, 0, 1, N - 1);

      for (j = 1; j < N; ++j)
        {
          gsl_vector_view w = gsl_vector_subvector(&work.vector, 0, j);
          double ajj = gsl_matrix_get(A, j, j);
          double dval;

          for (i = 0; i < j; ++i)
            {
              double aii = gsl_matrix_get(A, i, i);
              double aji = gsl_matrix_get(A, j, i);
              gsl_vector_set(&w.vector, i, aji * aii);
            }

          v = gsl_matrix_subrow(A, j, 0, j); /* A(j,1:j-1) */
          gsl_blas_ddot(&v.vector, &w.vector, &dval);
          ajj -= dval;

          if (ajj == 0.0)
            {
              GSL_ERROR ("matrix is singular", GSL_EDOM);
            }

          gsl_matrix_set(A, j, j, ajj);

          if (j < N - 1)
            {
              double ajjinv = 1.0 / ajj;
              gsl_matrix_view m = gsl_matrix_submatrix(A, j + 1, 0, N - j - 1, j); /* A(j+1:n, 1:j-1) */
              v = gsl_matrix_subcolumn(A, j, j + 1, N - j - 1);                    /* A(j+1:n, j) */
              gsl_blas_dgemv(CblasNoTrans, -ajjinv, &m.matrix, &w.vector, ajjinv, &v.vector);
            }
        }

      /* save ||A||_1 in upper right corner */
      gsl_matrix_set(A, 0, N - 1, anorm);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_ldlt_solve (const gsl_matrix * LDLT,
                       const gsl_vector * b,
                       gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_ldlt_svx(LDLT, x);

      return status;
    }
}

int
gsl_linalg_ldlt_svx (const gsl_matrix * LDLT,
                     gsl_vector * x)
{
  if (LDLT->size1 != LDLT->size2)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (LDLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_const_view diag = gsl_matrix_const_diagonal(LDLT);

      /* solve for z using forward-substitution, L z = b */
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasUnit, LDLT, x);

      /* solve for y, D y = z */
      gsl_vector_div(x, &diag.vector);

      /* perform back-substitution, L^T x = y */
      gsl_blas_dtrsv (CblasLower, CblasTrans, CblasUnit, LDLT, x);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_ldlt_rcond (const gsl_matrix * LDLT, double * rcond,
                       gsl_vector * work)
{
  const size_t M = LDLT->size1;
  const size_t N = LDLT->size2;

  if (M != N)
    {
      GSL_ERROR ("LDLT matrix must be square", GSL_ENOTSQR);
    }
  else if (work->size != 3 * N)
    {
      GSL_ERROR ("work vector must have length 3*N", GSL_EBADLEN);
    }
  else
    {
      int status;
      double Anorm;    /* ||A||_1 */
      double Ainvnorm; /* ||A^{-1}||_1 */

      if (N == 1)
        Anorm = fabs(gsl_matrix_get(LDLT, 0, 0));
      else
        Anorm = gsl_matrix_get(LDLT, 0, N - 1);

      *rcond = 0.0;

      /* don't continue if matrix is singular */
      if (Anorm == 0.0)
        return GSL_SUCCESS;

      /* estimate ||A^{-1}||_1 */
      status = gsl_linalg_invnorm1(N, ldlt_Ainv, (void *) LDLT, &Ainvnorm, work);

      if (status)
        return status;

      if (Ainvnorm != 0.0)
        *rcond = (1.0 / Anorm) / Ainvnorm;

      return GSL_SUCCESS;
    }
}

/* compute 1-norm of symmetric matrix stored in lower triangle */
static double
ldlt_norm1(const gsl_matrix * A)
{
  const size_t N = A->size1;
  double max = 0.0;
  size_t i, j;

  for (j = 0; j < N; ++j)
    {
      gsl_vector_const_view v = gsl_matrix_const_subcolumn(A, j, j, N - j);
      double sum = gsl_blas_dasum(&v.vector);

      /* now add symmetric elements from above diagonal */
      for (i = 0; i < j; ++i)
        {
          double Aij = gsl_matrix_get(A, i, j);
          sum += fabs(Aij);
        }

      if (sum > max)
        max = sum;
    }

  return max;
}

/* x := A^{-1} x = A^{-t} x, A = L D L^T */
static int
ldlt_Ainv(CBLAS_TRANSPOSE_t TransA, gsl_vector * x, void * params)
{
  int status;
  gsl_matrix * A = (gsl_matrix * ) params;
  gsl_vector_const_view diag = gsl_matrix_const_diagonal(A);

  (void) TransA; /* unused parameter warning */

  /* compute L^{-1} x */
  status = gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasUnit, A, x);
  if (status)
    return status;

  /* compute D^{-1} x */
  gsl_vector_div(x, &diag.vector);

  /* compute L^{-t} x */
  status = gsl_blas_dtrsv(CblasLower, CblasTrans, CblasUnit, A, x);
  if (status)
    return status;

  return GSL_SUCCESS;
}
